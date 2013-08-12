#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <values.h>
#include <stack>
#include <signal.h>
#include <unistd.h>

#include <typeinfo>
#include <math.h>
#include <fstream>
#include <iostream>
#include <string>
#include "mydll.h"


normal::normal(int __dim) : dim(__dim){
  assert(dim > 0);
  internalsUpToDate = false;

  _mean.destructiveLazyResize(1, dim);
  _cov.destructiveLazyResize(dim, dim);
  _invcov.destructiveLazyResize(dim, dim);
  _eigval.destructiveLazyResize(1, dim);
  _eigvec.destructiveLazyResize(dim, dim);

  _mean = 0;
  _cov.makeIdentity();
}

normal::~normal(){

}


// **************************************************
// *** direct element access
// **************************************************
//


double &normal::operator[](int i){
  assert(i >= 0 && i < dim);
  return _mean.rvalue(i);
}



double &normal::mean(int i){
  assert(i >= 0 && i < dim);
  internalsUpToDate = false;	// user might modify _mean value...
  return _mean.rvalue(i);
}



double &normal::cov(int i, int j){
  assert(i >= 0 && i < dim && j >= 0 && j < dim);
  internalsUpToDate = false;	// user might modify matrix...
  return _cov.value(i, j);
}



double normal::corr(int i, int j){
  assert(i >= 0 && i < dim);
  assert(j >= 0 && j < dim);

  assert(_cov.dvalue(i) > 0.0);
  assert(_cov.dvalue(j) > 0.0);
  
  return _cov.value(i, j) / (sqrt(_cov.dvalue(i)) * sqrt(_cov.dvalue(j)));
}


cmatrix &normal::mean(){
  internalsUpToDate = false;	// user might modify matrix...
  return _mean;
}


cmatrix &normal::cov(){
  internalsUpToDate = false;	// user might modify matrix...
  return _cov;
}




/*
// **************************************************
// *** project to lower dimal normal (good for display purposes)
// **************************************************
//

void normal::marginalize(normal &theNormal, int *projections){
  // --------------------------------------------------
  //
  // check various  stuff
  //

  assert(theNormal.dim > 0);
  for (int i = 0; i < theNormal.dim; i++)
    assert(projections[i] >= 0 && projections[i] < dim);

  // the following might be avoided by dynamically allocating more memory
  assert(memoryAllocated);
  assert(theNormal.memoryAllocated);

  for (int i = 0; i < theNormal.dim; i++){
    theNormal._mean(i) = _mean(projections[i]);
    for (int j = 0; j < theNormal.dim; j++)
      theNormal.cov(i, j) = cov(projections[i], projections[j]);
  }

  theNormal.internalsUpToDate = false;
  // calculateInternals(true);
}
*/

// **************************************************
// *** calculate likelihood of a data item
// **************************************************
//

double normal::unNormalizedLogLikelihood(const cmatrix &dat){
  assert(dim > 0);
  cmatrix diff(dat , CMAT_COPY);
  cmatrix bufa, bufb;

  diff -= _mean;;
  bufa.gemm(_invcov, diff.trans());
  bufb.gemm(diff, bufa);

  return - 0.5 * bufb.rvalue(0);
}




double normal::unNormalizedLikelihood(const cmatrix &dat){
  return exp(unNormalizedLogLikelihood(dat));
}

double normal::logLikelihood(const cmatrix &dat){
  return logNormalizer() + unNormalizedLogLikelihood(dat);
}

double normal::likelihood(const cmatrix &dat){
  double logLik = logLikelihood(dat);
  if (logLik > 80.0) logLik = 80.0;
  return exp(logLik);
}

double normal::logNormalizer(){
  calculateInternals();

  return -0.5 * (dim * log(2.0 * M_PI) + log(fabs(_det)));
}

double normal::normalizer(){
  return exp(logNormalizer());
}


void normal::calculateInternals(bool force){
  if (!force && internalsUpToDate)
    return;

  _invcov.inverse(_cov, _det);
  _cov.symmEig_nr(_eigval, _eigvec);

  internalsUpToDate = true;
}






// **************************************************
// *** generate normally distributed samples
// **************************************************
//

double __dgasdev(){return gasdev();}

void normal::generateSample(cmatrix &targetData, int numSamples){
  assert(numSamples > 0);
  assert(dim == targetData.cols());

  calculateInternals();

  cmatrix evecT(_eigvec.trans(), CMAT_VIEW);
  cmatrix evalt(_eigval, CMAT_COPY);
  cmatrix randvec(1, dim);
  cmatrix bufa;

  evalt.op(sqrt);
  for (int i = 0; i < numSamples; i++){	// loop through samples
    //cout<<"xi"<<endl<<targetData.row(i)<<endl;
    randvec.makeRand(__dgasdev);
    randvec.dotTimes(evalt);
    targetData.row(i).add(_mean, bufa.gemm(randvec, evecT));
  }
}


/*


// **************************************************
// *** fit data (calculate new Gaussian)
// **************************************************
//



void normal::fitData(cmatrix &targetData, double *weights ){
  assert(targetData.rows() >= dim);
  assert(dim == targetData.cols());

  calculateInternals();

  // --------------------------------------------------
  //
  // this is a check on our ability to fit a Gaussian
  //

  bool doDummyCheck(true);


  double xBefore   = 0.0;
  if (doDummyCheck){
    for (int k = 0; k < targetData.numDataItems; k++)
      if (weights)
	xBefore += (weights[k] * logLikelihood(&(targetData.theData[k][0])));
      else
	xBefore += logLikelihood(&(targetData.theData[k][0]));
  }
  
  // --------------------------------------------------
  //
  // extract basic statistics
  //

  double sum[dim];
  double crossSum[dim][dim];
  double totalFact = 0.0;
  
  for (int i = 0; i < dim; i++){
    sum[i] = 0.0;
    for (int j = 0; j < dim; j++)
      crossSum[i][j] = 0.0;
  }

  for (int k = 0; k < targetData.numDataItems; k++){

    double fact(1.0);
    if (weights)
      fact = weights[k];

    for (int i = 0; i < dim; i++){
      sum[i] += fact * targetData.theData[k][i];

      for (int j = 0; j < dim; j++)
	crossSum[i][j] += fact * (targetData.theData[k][i] *
				  targetData.theData[k][j]);
    }

    totalFact += fact;
  }

  assert(totalFact > 0.0);

  // --------------------------------------------------
  //
  // compute new Gaussian
  //

  double multiplier = 1.0 / totalFact;

  for (int i = 0; i < dim; i++)
    _mean.value(i) = multiplier * sum[i];

  for (int i = 0; i < dim; i++)
    for (int j = 0; j < dim; j++)
      _cov.value(i, j) = 
	(multiplier * crossSum[i][j]) -
	(_mean.value(i) * _mean.value(j));
  
  // --------------------------------------------------
  //
  // update statistics
  //
  
  calculateInternals(true);	// forced calculation


  if (doDummyCheck){
    double xAfter   = 0.0;
    for (int k = 0; k < targetData.numDataItems; k++)
      if (weights)
	xAfter += (weights[k] * logLikelihood(&(targetData.theData[k][0])));
      else
	xAfter += logLikelihood(&(targetData.theData[k][0]));
    if (xBefore > xAfter + 0.00001){
      cerr << endl << "### DummyCheck failed: " << xBefore << " -> " 
	   << xAfter << " ###" << endl << endl;
      putc(7, stderr);
    }
  }
}





// **************************************************
// *** compute weighted error over data set
// **************************************************
//

double
normal::calculateLogLikelihood(cmatrix &targetData, double *weights)
{

  assert(0);			// we never tested this routine...

  // --------------------------------------------------
  //
  // check
  //

  assert(targetData.numDataItems >= dim);
  assert(dim == targetData.dim);
  
  
  // --------------------------------------------------
  //
  // iterate through data and calculate log likelihood
  //


  double logLik = 0.0;

  for (int k = 0; k < targetData.numDataItems; k++){

    double fact(1.0);
    if (weights)
      fact = weights[k];

    logLik += fact * logLikelihood(targetData.theData[k]);
  }

  // --------------------------------------------------
  //
  // return (weighted) log likelihood
  //

  return logLik;

}
  

// **************************************************
// *** test of factorization (needs improvement)
// **************************************************
//


bool
normal::factorTest(data &targetData, double *weights )
{
  double LogLikCond[dim];
  double LogLikIndep[dim];
  double ScoreCond = 0.0, ScoreIndep = 0.0;
  
  for (int i = 0; i < dim; i++){
    ScoreCond += (LogLikCond[i]  = 
		  calculateConditionalLogLikelihood(targetData, i, weights));
    ScoreIndep += (LogLikIndep[i] = 
		   calculateIndependentLogLikelihood(targetData, i, weights));
  }

  doFactor = (ScoreCond < ScoreIndep);

  cerr << "*** factor test? ***    --->";
  if (doFactor)
    cerr << " FACTORIZE!";
  else
    cerr << " don't factorize";
  cerr << endl;
  cerr << "original:";
  for (int i = 0; i < dim; i++)
    cerr << " " << LogLikCond[i];
  cerr << " -> " << ScoreIndep << endl;
  cerr << "factored:";
  for (int i = 0; i < dim; i++)
    cerr << " " << LogLikIndep[i];
  cerr << " -> " << ScoreCond << endl;
  cerr << endl;

  return doFactor;
}


// **************************************************
// *** compute the error of data under Gaussian
// **************************************************
//


double normal::calculateConditionalLogLikelihood(cmatrix &targetData, 
						 int d, 
						 cmatrix &weights ){
  assert(targetData.rows() > 0); // not necessary, but for precaution
  assert(targetData.cols() == dim);


  double sum(0.0), norm(0.0);

  if (0)
    for (int k = 0; k < targetData.numDataItems; k++){
      sum  += (weights[k] * unNormalizedLogLikelihood(targetData.theData[k]));
      norm += weights[k];
    }
  else
    for (int k = 0; k < targetData.rows(); k++){
      sum  += unNormalizedLogLikelihood(targetData.theData[k]);
      norm += 1.0;
    }

  assert(norm > 0.0);
  sum /= norm;
  return sum;
}




double normal::calculateIndependentLogLikelihood(cmatrix &targetData, 
						 int d,
						 cmatrix &weights ){
  assert(targetData.rows() > 0); // not necessary, but for precaution
  assert(targetData.cols() == dim);

  double sum(0.0), norm(0.0);
  double testDat[dim];

  for (int i = 0; i < dim; i++)
    testDat[i] = _mean.rvalue(i);

  for (int k = 0; k < targetData.rows(); k++){
    testDat[d] = targetData.value(k, d);
    if (0){
      sum  += (weights[k] * unNormalizedLogLikelihood(testDat));
      norm += weights[k];
    }
    else{
      sum += unNormalizedLogLikelihood(testDat);
      norm += 1.0;
    }
  }

  assert(norm > 0.0);
  sum /= norm;
  return sum;  
}
*/
