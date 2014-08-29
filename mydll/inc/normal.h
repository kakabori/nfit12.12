#ifndef __NORMAL_HH
#define __NORMAL_HH



// ***********************************************************************
// ***********************************************************************
// ***********************************************************************
// 
// normal - class for normal distribution
// 
// ***********************************************************************
// ***********************************************************************
// ***********************************************************************



class normal{

public:

  // constructor(s)
  normal(int __dim);

  // destructor(s)
  ~normal();

  // copy operator
  //normal(const normal& x); 
  //normal operator=(normal inNormal);

  // direct access to the elements of the normal distribution
  double &operator[](int i);
  double &mean(int i);
  double &cov(int i, int j);
  double corr(int i, int j);
  cmatrix &mean();
  cmatrix &cov();

  // calculate inverse, determinnant, eigen values and eigen vectors
  void calculateInternals(bool force = false);

  // compute likelihood
  double unNormalizedLogLikelihood(const cmatrix &dat);
  double unNormalizedLikelihood(const cmatrix &dat);
  double logLikelihood(const cmatrix &dat);
  double likelihood(const cmatrix &dat);
  double logNormalizer();
  double normalizer();

  // print operator
  //friend ostream& operator<<(ostream& out, const normal& n);

  // getDim()
  int getDim() { return dim; }

  // project to lower dimal mixture (good for display purposes)
  // void marginalize(normal *theNormal, int *projections);

  // sample
  void generateSample(cmatrix &targetData, int numSamples = 1);

  // fit Gaussian
  //void fitData(data &targetData, double *weights = NULL);

  // evaluating the error
  //double calculateLogLikelihood(data &targetData, double *weights = NULL);

  // factorization test
  //  bool factorTest(data &targetData, double *weights = NULL);


private:
  // internal data structures
  int dim;

  cmatrix _mean, _cov, _invcov, _eigval, _eigvec;
  double _det;

  bool internalsUpToDate;

  /*
  double calculateConditionalLogLikelihood(data &targetData, int d, 
					   double *weights = NULL);
  double calculateIndependentLogLikelihood(data &targetData, int d, 
					   double *weights = NULL);
  */
};

#endif
