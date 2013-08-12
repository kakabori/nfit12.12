#define PI 3.1415926535897932384626433832795

double xlow = 0;
double xhigh = 2048;
int nocz = 1;
int noscale = 1;
int kiyomask = 0;
double NKfilter = 0.0;
double dupe = 1;
double backgroundSigmaSquare = 10;
double aFactor = 0.2;
int updateSigma = 1;
int startflag = 0;
int stopflag = 0;
int nthreads = 2;

/*FunSupport parameters*/
double g_spRotatedAbserr = 1e-20;
double g_spRotatedRelerr = 5e-4;
double g_spRotatedMindx = 0.001;
double g_spRotatedMaxdx = 10;
double g_spStrFctAbserr = 1e-20;
double g_spStrFctRelerr = 1e-3;
double g_spStrFctMindx = 0.001;
double g_spStrFctMaxdx = 10;
double g_spsumnAbserr = 1e-20;
double g_spsumnRelerr = 1e-4;
double g_spsumnMindx = 0.001;
double g_spsumnMaxdx = 10;
double g_spHrAbserr = 1e-20;
double g_spHrRelerr = 1e-4;
double g_spHrMindx = 0.01;
double g_spHrMaxdx = 100000;

//For qag routine
//int g_workspaceSize = 10000; //workspace size
double g_epsabs = 0; //absolute error tolerance
double g_epsrel = 1e-8; //relative error tolerance
//int g_key = 6; //key
