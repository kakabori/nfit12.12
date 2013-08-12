#ifndef _MYDLL_HEADER_
#define _MYDLL_HEADER_

#ifndef EXPORT
#define EXPORT
#endif

#define MYDLL_API extern EXPORT

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <iostream>
#include <assert.h>
using namespace std;

typedef double(*Integrand)(void *, double);
typedef double(*CHIFunc)(double *);
typedef double(*E2eOp)(double *, double *);

#define YFPI 3.1415926535897932384626433832795


class HankelIntegrator
{
  //  double ydb[800];
  //  double lastb;
  //  double *cutoff;
 public:
  //	HankelIntegral(){lastb=-1;cutoff=NULL;}
  double integrate(Integrand func, void *instance, double b, double cutoff);
  double integrate(double (*func)(double), double b, double cutoff);
};

class Qromb
{
private:
  double* hForExterp;
  double* sForExterp;
  double sum;
  double eps;
  int iteration;
  int jmax;
  int initNumOfIter;
  int numForExterp;
  void *instance;
  Integrand func;
  double trapzd(double a, double b, int n);
  double trapzd(Integrand f, double a, double b, int n);
  void inith(double);
  void default_setup();
public:
  Qromb(){default_setup();}
  Qromb(Integrand f, void *p){default_setup(); func=f; instance=p;}
  double calc(double a, double b);
  double calc(Integrand f, double a, double b);
  void setIntegrand(Integrand f){func=f;}
  void setInstance(void *p){instance=p;}
  void custom_setup(double a);
};

//This creates the type LmdifFunc, for "pointer to function (of these arguments) returning int," which can be used like LmdifFunc pfcn;
typedef int (*LmdifFunc)(int,int,double *,double *,void *,void *);

class Lmdif {
public:
  Lmdif();
  ~Lmdif(){cleanup();}
  void setup(LmdifFunc,int, int, int, double,	double, double, double, double*, void *);
  int fit();
  void init(double *,int);
  void init(double *b[],int);
  void covar_cminpack();
private:
  LmdifFunc pfcn;
  int m,n,ldfjac,maxfev;
  int *ipvt;
  double *x, *fvec, **fjac, *qtf, *diag, *wa1, *wa2, *wa3, *wa4;
	double *fjac_cminpack;    //to use cminpack covar function, **fjac stored as one dimensional array, added by KA, 9/20/2012
	double ftol, xtol, gtol, factor, *epsfcn;
  void *client;
  void cleanup();
};



typedef double (*YLMinimizationFunc)(double *x,void *,void *);
typedef void (*YLMinimizationDFunc)(double *x, double *g, void *,void *);

enum MinMethod {DFPVM, SUBPLEX};
struct DfpminCtrl{
  double tolx;
  double gtol;
  double eps;
  double alf;
  double stpmx;

  int maxiter;
  int nvar;
  YLMinimizationFunc func;
  YLMinimizationDFunc dfunc;

  DfpminCtrl();
};

struct YLSUBPLEXctrl {
  double alpha, /* reflection coefficient       alpha > 0 */
    beta,       /* contraction coefficient   0 < beta < 1 */
    gamma,      /* expansion coefficient         gamma > 1 */
    delta,      /* shrinkage (massive contraction) coefficient 0 < delta < 1 */
    psi,        /* simplex reduction coefficient   0 < psi < 1 */
    omega;      /* step reduction coefficient      0 < omega < 1 */

  int nsmin,
    nsmax,
    irepl,
    ifxsw;
  double bonus, fstop;
  int nfstop, nfxe;
  double fxstat[4], ftest;
  bool  initx, newx;
  /* addtional */
  double tol;
  double maxnfe;
  int nvar;
  void *client;
  YLMinimizationFunc func;

  YLSUBPLEXctrl();
};

class YLMinimization{
 public:
  double *x;
  void *client;
  int iter;
  int n;
  double val;
  void *ctrl;
  int method;
  int fit();
  int setup(void*,int);
  YLMinimization(){client=ctrl=NULL;x=NULL;}
};









MYDLL_API double d_sign(double a, double b);

MYDLL_API double bessj0(double x);
MYDLL_API double polint3(double z,double *ya);
MYDLL_API double polint3_log(double xmxa,float *ya);
MYDLL_API double polint30_log(double xmxa,float *ya);
MYDLL_API double polint4t(double xmxa,double *ya);
MYDLL_API double nr_erf(double x);
MYDLL_API double expint(double x);
MYDLL_API void amoeba(double **p, double *y,int ndim, double ftol, double (*funk)(double *), int *nfunk);
MYDLL_API double poidev(double xm, int *idum);
MYDLL_API double YFfindmin(double *x,double *step, int *xpin,int varnum, CHIFunc chifunc, int iter, void *client);
MYDLL_API int mydllNearNeighbour(double *xa, int incx, int n, double q);



MYDLL_API int yl_subplex_subplx(YLMinimizationFunc f,
				 int n,
				 double tol,
				 int maxnfe,
				 int mode,
				 double *scale,
				 double *x,
				 double *fx,
				 int *nfe,
				 double *work,
				 int *iwork,
				 int *iflag,
				 void *ctrl);




template<class T> inline T yfmax(T a, T b){return a>b? a : b; }
template<class T> inline T yfmin(T a, T b){return a<b? a : b; }

template<class T>
struct YFDuple{
  T x,y;
  YFDuple(){}
  YFDuple(T a, T b){set(a,b);}
  void set(T a, T b){x=a; y=b;}
  YFDuple<T> operator+(YFDuple<T>& a){YFDuple<T> b(x+a.x, y+a.y); return b;}
  YFDuple<T>& operator+=(YFDuple<T>& a){x+=a.x; y+=a.y; return *this;}
  YFDuple<T> operator-(YFDuple<T>& a){YFDuple<T> b(x-a.x, y-a.y); return b;}
  YFDuple<T>& operator-=(YFDuple<T>& a){x-=a.x; y-=a.y; return *this;}
  T norm2(){return x*x+y*y;}
};

template<class T>
class YFRect{
  YFDuple<T> dll;
  YFDuple<T> dur;
 public:
  T& llx(){return dll.x;}
  T& urx(){return dur.x;}
  T& lly(){return dll.y;}
  T& ury(){return dur.y;}
  T width(){return dur.x-dll.x;}
  T height(){return dur.y-dll.y;}
  YFDuple<T> ll(){return dll;}
  YFDuple<T> ur(){return dur;}
  YFDuple<T> lr(){YFDuple<T> a(urx(),lly()); return a;}
  YFDuple<T> ul(){YFDuple<T> a(llx(),ury()); return a;}
  T midx(){return (llx()+urx())/2;}
  T midy(){return (lly()+ury())/2;}
  void setll(YFDuple<T> a){dll=a;}
  void setur(YFDuple<T> a){dur=a;}
  void setll(T a, T b){dll.x=a, dll.y=b;}
  void setur(T a, T b){dur.x=a, dur.y=b;}
  void set(T a, T b, T c, T d);
  void shift(YFDuple<T>);
  bool overlap(YFRect<T>& a, YFRect<T>& b);
  YFDuple<T> overlap_ll(YFRect &);
  YFDuple<T> overlap_ur(YFRect &);
  YFRect<T> overlap(YFRect<T>& rect_b);
  T area(){return width()*height();}
};

template<class T>
void YFRect<T>::set(T a, T b, T c, T d){
  if(a<c){ dll.x=a; dur.x=c;}
  else{    dll.x=c; dur.x=a;}
  if(b<d){ dll.y=b; dur.y=d;}
  else{    dll.y=d; dur.y=b;}
}

template<class T>
bool YFRect<T>::overlap(YFRect<T>& a, YFRect<T>& b){
  T x1, y1, x2, y2;
  x1 = yfmax( llx(), a.llx() );
  y1 = yfmax( lly(), a.lly() );
  x2 = yfmin( urx(), a.urx() );
  y2 = yfmin( ury(), a.ury() );
  b.set(x1, y1, x2, y2);
  return ( x1 < x2  && y1 < y2 );
}

template<class T>
YFDuple<T> YFRect<T>::overlap_ll(YFRect<T>& rect_b){
  YFDuple<T> point;
  point.x=yfmax(llx(), rect_b.llx());
  point.y=yfmax(lly(), rect_b.lly());
  return point;
}

template<class T>
YFDuple<T> YFRect<T>::overlap_ur(YFRect<T>& rect_b){
  YFDuple<T> point;
  point.x=yfmin(urx(), rect_b.urx());
  point.y=yfmin(ury(), rect_b.ury());
  return point;
}

template<class T>
YFRect<T> YFRect<T>::overlap(YFRect<T>& rect_b){
  YFRect<T> overlap;

  setll(this->overlap_ll(rect_b));
  setur(this->overlap_ur(rect_b));
  return *this;
}

template<class T>
void YFRect<T>::shift(YFDuple<T> a){
  dll -= a; dur -= a;
}

MYDLL_API void fct2d(double f[], int nrows, int ncols);
MYDLL_API void ifct2d(double f[], int nrows, int ncols);

template<class T> inline T yl_max(T a, T b){return a>b?a:b ;}
template<class T> inline T yl_min(T a, T b){return a<b?a:b ;}
template<class T> inline T sqr(T a){ return a*a; }

template<class T> void tvec(T** a, int n){ *a = (T*)calloc(n, sizeof(T)); }

template<class T> void tvecr(T** a, int n){ *a = (T*)realloc(*a, n*sizeof(T)); }


template<class T> void tmat(T*** a, int n, int m){
  *a = (T**)malloc(n* sizeof(T*));
  (*a)[0] = (T*)calloc(m*n, sizeof(T));
  for(int i=1; i<n; i++){
    (*a)[i] = (*a)[i-1] + m;
  }
}
template<class T> void freetmat(T** a){
  free(a[0]); free(a);
}

template<class T> inline void tswap(T&a, T&b){T t=a; a=b; b=t;}

template<class T> inline T d_sign(T a, T b){
  T x = (a>=0) ? a : -a;
  return (b>=0) ? x : -x;
}

template<class A, class B> void tconvert(int n, A *a, int inca, B *b, int incb){
  B *bend = b + n * incb ;
  while(b != bend){
    *b = (B)*a;
    b += incb;
    a += inca;
  }
}


void locate(float *xx, int n, float x, int *j);


struct dduple {
  double x, y;
};
#define degToRad(x) ((x) / 180.0 * M_PI)
#define radToDeg(x) ((x) * 180.0 / M_PI)
double regularTruncate(double x, double min, double max);
double cyclicTruncate(double x, double min, double max);
double L2Dist(double x1, double y1, double x2, double y2);
bool boundCheck(double value, double min, double max);
double normalSample(double mean, double variance);

void svdfit(int ndata, double *a, int ma, double *u, double *v, double *w, double *b);


float gasdev();
int bisearch(float *xx, int n, float x);
float urand0();
double pythag(double a, double b);

#include <../mydll/inc/yl_blas.h>
#include <../mydll/inc/cmatrix.h>
#include <../mydll/inc/uniform.h>
#include <../mydll/inc/normal.h>
#include <../mydll/inc/matcher.h>
#endif
