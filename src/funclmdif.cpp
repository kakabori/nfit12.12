/*  11/27/06  This version uses only backsig for estimated errors
    Comparison of data - fits indicates that errors are independent
    of the size of the signal, so that feature is set to zero.
    See line 143 */
/*  9/14/2012 This is an updated version with revisions by Jonathan
    Michel. This implementation now supports multithreaded calciulation
    of chi squared to take advantage of multiple cores to speed up
    the calculation. */
#include <sstream>
#include <tcl.h>
#include <tk.h>
#include <pthread.h>

#include "nfit.h"
#include "funclmdif.h"
#include "globalVariables.h"

extern void updatelinks(double xisquare, char *chain);
extern double refine(double theor, double exper);
extern Para g_ParaStruct;
extern Tcl_Interp *NKinterp;

//////////////////////////////////////////////////////////////
// compute I_e - I_{model} to be used for least square optimization
// lmdif which is linked upon compilation - see makefile
//////////////////////////////////////////////////////////////

void FuncLmdif::setPara(Para *p)
{
  int i;
  para = p;
  for(i = 0; i < Para::s_numParams; i ++) {
    bestParams[i] = p->getValue(i);
  }
}

/****************************************************************************************
This function logs the combination of parameters that yields the optimal
chi squared value, along with the corresponding chi squared value.
****************************************************************************************/
void FuncLmdif::logBest(double chisq, char *chain)
{
  int i;
  bestChisq = chisq;
  for(i = 0; i < Para::s_numParams; i ++) {
    bestParams[i] = para->getValue(i);
  }
  sprintf(bestChain, "Best values: %s", chain);
}

/****************************************************************************************
This function retirieves the combination of parameters that yielded the
optimal chi squared value, along with the corresponding chi squared value.
****************************************************************************************/
double FuncLmdif::recoverBestParams(Para *p)
{
  int i;
  for(i = 0; i < Para::s_numParams; i ++){
    p->setValue(bestParams[i], i);
  }
  return bestChisq;
}

/****************************************************************************************
This function returns a string reporting the optimal value for each
parameter.
****************************************************************************************/
char *FuncLmdif::getBestChain()
{
  return bestChain;
}


/**************************************************************************************************
typedef int (*LmdifFunc)(int, int, double *, double *, void *, void *);
Lmdif object gets linked to this function at the beginning of a fit

n < m
m: an input variable set to the number of functions (= number of data points)
n: an input variable set to the number of variables (= number of free parameters)
par: an array of length n. On input par must contain an initial estimate of the
solution vector. On output par contains the final estimate of the solution vector.
fvec: an output array of length m which contains the functions evaluated at the output x.

**************************************************************************************************/
int FuncLmdif::WrapperFunclmdif(int m, int n, double *par, double *fvec, void *client, void* ctrl)
{
  return ((FuncLmdif*)client)->funclmdif(m, n, par, fvec, ctrl);
}

/**************************************************************************************************
This function calculates the sum of squares

n: number of parameters to be fitted
par: an array that stores initial values of the free parameters
fvec: an array to store the residuals (= weighted difference between the model and a data point)
ctrl: ?
**************************************************************************************************/
int FuncLmdif::funclmdif(int m, int n, double *par, double *fvec, void* ctrl)
{

  if(stopflag) return 0;
  printf("funclmdif ");
  char chain[256]={0};
  //A set of structures to pass information to each chi squared worker thread.
  MCThreadPars **threadPars;
  //pthread_t structures to hold information about each worker thread
  pthread_t *threads;
  double sum = 0;
  double chisq;
  size_t nslice = data->qs.size();
  int i;
  size_t minIter, maxIter, iterChunk, sliceWidth;
  ModelCalculator *myMC;

  // loop through n free parameters and update
  // corresponding ModelCalculator parameters
  for (int k = 0; k < n; k++) {
    //double tmp;
    //tmp = par[k];

    // for variables other than bc2b, don't accept a negative value
		if (k != Var_bc2b) {
		  par[k] = fabs(par[k]);
		  printf("%g ", par[k]);
		  sprintf(chain, "%s %g ", chain, par[k]);
		  *(para->xp[k]) = par[k];
		  mc->setpara(par[k], para->idx[k]);
		// bc2b can be negative when the beam is below the detector edge
		} else {
			printf("%g ", par[k]);
			sprintf(chain, "%s %g ", chain, par[k]);
			*(para->xp[k]) = par[k];
			mc->setpara(par[k], para->idx[k]);
		}
  }

  fflush(stdout);

  /*Set up the multithreaded calculation of chi squared. If the number of
    slices, then the number of threads should simply be set to the number of
    slice. Prepare the threads and structures containing additional tools
    needed by the worker threads.
  */
  if(nslice < (unsigned int) nthreads) nthreads = nslice;

  threadPars = (MCThreadPars **) malloc(nthreads*sizeof(MCThreadPars *));
  threads = (pthread_t *) malloc(nthreads * sizeof(pthread_t));

  /*Determine the number of slices given to the first n - 1 worker threads,
    given n total worker threads. This policy is designed to guarantee that
    the final thread does not take on a signficantly different number of
    slices, thus preventing overcrowding in some threads and under-utilization
    of others.
  */
  iterChunk = nslice / nthreads;
  if((nslice % nthreads) > (nthreads / 2)) iterChunk ++;
  //Get the slice width so it can be known how to divide up the fvec array
  //amongst the threads.
  sliceWidth = data->qs[0].sdata.size();

  for(i = 0; i < nthreads; i++){
    threadPars[i] = (MCThreadPars *) malloc(sizeof(MCThreadPars));
  }

  /*Prepare the data needed for each thread. Each thread must have its own
    ModelCalculator object because in calculating theoretical values, certain
    variables used by a ModelCalculator object are set for a specific qz value.
    Further, tell the thread the first and last qz slices it should consider.
    Also give each thread a reference to the fitting parameter information and
    the experimental values stored within the Data data structure, which is
    helpfully named data.
  */
  for(i = 0; i < nthreads; i ++){
    minIter = i * iterChunk;
    maxIter = (i < nthreads - 1) ? minIter + iterChunk - 1 : nslice - 1;
    threadPars[i]->minIter = minIter;
    threadPars[i]->maxIter = maxIter;
    threadPars[i]->fvec = fvec + i * iterChunk * sliceWidth;
    threadPars[i]->m = m;
    myMC = (i) ? new ModelCalculator(*mc) : mc;
    threadPars[i]->myMC = myMC;
    threadPars[i]->data = data;
    threadPars[i]->para = para;
    pthread_create(&threads[i], NULL, &FuncLmdif::modelCalcThread, (void *) threadPars[i]);
  }
  //Join each worker thread to the parent to ensure that the parent waits for
  //all worker threads to finish before continuing.
  for(i = 0; i < nthreads; i ++){
    pthread_join(threads[i], NULL);
  }
  //Reap chi squared values calculated by each thread, then free all unneeded
  //data structures.
  for(i = 0; i < nthreads; i ++){
    sum += threadPars[i]->mySum;
    if(threadPars[i]->m == -1) m = -1;
    if(i){
      myMC = threadPars[i]->myMC;
      myMC->cleanup();
      free(myMC);
    }
    free(threadPars[i]);
  }

  free(threadPars);
  free(threads);

  //Finalize the chi squared caluclation and report current parameter values.
  chisq = sum/m;
  printf("Xr: %g  / %d = %g  /  %g \n",sum,m, chisq,backgroundSigmaSquare); fflush(stdout);
  sprintf(chain, "%s Xr: %g  / %d = %g ", chain, sum, m, chisq);
  updatelinks(chisq, chain);
  if(chisq < bestChisq) logBest(chisq, chain);
  // update associated image with 'data'
  data->writeimg();
  data->writefrm((char *) "frm.dat", para);

  /*
  // update paramArray, which is a TCL global variable, whose elements are
  // linked to the fitting panel entry fields
  Tcl_Interp *interp = NKinterp;
  std::ostringstream strs;
  for (int i = 0; i < Para::s_numParams; i++) {
    strs << g_ParaStuct.getValue(i);
    std::string newValue = strs.str();
    Tcl_SetVar2(interp, "paramArray", Varname[i], newValue.c_str(), TCL_GLOBAL_ONLY);
  }
  */

  return 0;
}

/****************************************************************************************
This function is run within a worker thread to calculate a chi sqared value
for a subset of all slices being fit.
****************************************************************************************/
void* FuncLmdif::modelCalcThread(void *args){
  MCThreadPars *myPars = (MCThreadPars *) args;
  ModelCalculator *myMC = myPars->myMC;
  Data *data = myPars->data;
  Para *para = myPars->para;
  size_t minIter = myPars->minIter;
  size_t maxIter = myPars->maxIter;
  double *fvec = myPars->fvec;
  int m = myPars->m;
  size_t i;
  double sum1=0,sum2=0, sum3=0, sum4 =0, sum5= 0, sum=0,scale, bias;
  vector<DataPoint>::iterator it;
  vector<DataPoint>::iterator itend;

  for (i = minIter; i <= maxIter; i++) {
    /* iterate through all the slices */

    // tell mc which slice to calculation for this iteration
    myMC->setSliceParameter( para->setup.getqz(data->qs[i].qz) );

    // specify the qx range for model calculation
    it = data->qs[i].sdata.begin();
    itend = data->qs[i].sdata.end();
    myMC->QxSlice( para->setup.getqr(it->qx), para->setup.getqr((itend-1)->qx) );

    /* Obtain linear parameters c(z) background adjustment and
	the overall scale s(z) that appear in Eq. 5.1 and Eq. 2.1,
	given fixed values of all the other parameters that determine Sccd.
	Compute the summations for determining linear parameters,*/
    sum1=sum2=sum3=sum4=sum5=0;
    for(; it != itend; it++){
      double vcal;
      double sigma2;

      /*Notation in Numerical Recipes, chapter 15.2, Eq. (15.2.4)
        vcal = x_i
        sigma2 = sigma^2
        sum1 = S_{xy}
        sum2 = S_{xx}
        sum3 = S_x
        sum4 = S_y
        sum5 = S*/
      vcal = myMC->getCCDStrFct( para->setup.getqr(it->qx) );
      sigma2 = it->sigma * it->sigma;
      sum1 += it->inte * vcal / sigma2;
      sum2 += vcal * vcal / sigma2;
      sum3 += vcal / sigma2;
      sum4 += it->inte / sigma2;
      sum5 += 1 / sigma2;
    }

    // obtain 'scale' and 'bias' using sum*
    if( nocz ){
      scale = sum1 / sum2;
      bias=0;
    } else {
        double Delta_NR = sum5*sum2 - sum3*sum3;        // -Delta defined in Numerical Recipes
        scale = (sum1*sum5 - sum4*sum3) / Delta_NR;     // s(z) (or phi(z))
        bias = -(sum2*sum4 - sum3*sum1) / Delta_NR;     // c(z) parameter
        double sigma_scale = sqrt(sum5 / Delta_NR);     // uncertainty in s(z)
        double sigma_bias = sqrt(sum2 / Delta_NR);      // uncertainty in c(z)

        // store uncertainties in Data object
    	data->qs[i].sigma_scale = sigma_scale;
    	data->qs[i].sigma_bias = sigma_bias;
    }
    if( noscale ) {
      scale = 1;
      bias = 0;
    }

    // store 'scale' and 'bias' in Data object
    data->qs[i].scale = scale;
    data->qs[i].bias  = bias;



    // compute (I_e - I_{model})/(\sigma_p^2) and write them to fvec
    // also update 'data' with model values
    it = data->qs[i].sdata.begin();
    for(; it != itend; it++){
      double tmp;
      it->cal = tmp = myMC->getCCDStrFct( para->setup.getqr(it->qx) ) * scale - bias;
  //do not include outstanding points (refinement)
      if (NKfilter != 0 && refine(tmp, it->inte) > NKfilter) {tmp=0; m-=1;}
      else {tmp = ( tmp - it->inte ) / it->sigma;}  //looks like the command used.
      /*      else {tmp = ( tmp - it->inte ) / sqrt(backgroundSigmaSquare);}
	      gives much smaller chi2 b/c backgroundSigmaSquare is fixed at 10.
            it->sigma comes from dataset.cxx  */
      *(fvec++) = tmp;
      sum += tmp * tmp;
    }
  }
  myPars->mySum = sum;
  return NULL;
}
