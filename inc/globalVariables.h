#ifndef GUARD_GLOBALVARIABLES_H
#define GUARD_GLOBALVARIABLES_H

extern double xlow;
extern double xhigh;
extern int nocz;
extern int noscale;
extern int kiyomask;
extern double NKfilter;
extern double dupe;
extern double backgroundSigmaSquare;
extern double aFactor;
extern int updateSigma;
extern int startflag;
extern int stopflag;
extern int nthreads;

/*FunSupport parameters*/
//extern double g_spRotatedAbserr;
//extern double g_spRotatedRelerr;
//extern double g_spRotatedMindx;
//extern double g_spRotatedMaxdx;
extern double g_spStrFctAbserr;
extern double g_spStrFctRelerr;
extern double g_spStrFctMindx;
extern double g_spStrFctMaxdx;
extern double g_spsumnAbserr;
extern double g_spsumnRelerr;
extern double g_spsumnMindx;
extern double g_spsumnMaxdx;
extern double g_spHrAbserr;
extern double g_spHrRelerr;
extern double g_spHrMindx;
extern double g_spHrMaxdx;


//For qag routine
//extern int g_workspaceSize; //workspace size
extern double g_epsabs; //absolute error tolerance
extern double g_epsrel; //relative error tolerance
//extern int g_key; //key
extern double g_omega;

#endif
