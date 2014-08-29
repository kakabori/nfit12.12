#include <math.h>
#include <stdlib.h>
#include "mydll.h"

#define TOL 1.0e-10

void linfit_pre(double *x, double *y, double *sig, int ndata, int ma, int xn, double *u, double *b,void (*funcs)(double *,double *,int)){
  int i,j;
  double *afunc;
  double tmp;
  afunc=(double *)malloc(sizeof(double)*ma);
  for (i=0;i<ndata;i++) {
    (*funcs)(x+i*xn,afunc,ma);
    tmp=1.0/sig[i];
    for (j=0;j<ma;j++) (u+j*ndata)[i]=afunc[j]*tmp;
    b[i]=y[i]*tmp;
  }
  free(afunc);
}

double linfit_post(double *x, double *y, double *sig, int ndata, int ma, int xn, double *a, void (*funcs)(double *,double *,int)){
  double chisq=0.0, sum, *afunc, tmp;
  int i,j;
  afunc=(double *)malloc(sizeof(double)*ma);
  for (i=0;i<ndata;i++) {
    (*funcs)(x+i*xn,afunc,ma);
    for (sum=0.0,j=0;j<ma;j++) sum += a[j]*afunc[j];
    chisq += (tmp=(y[i]-sum)/sig[i],tmp*tmp);
  }
  free(afunc);
  return chisq;
}

void svdfit(int ndata,double *a,int ma,double *u,double *v,double *w, double *b){
  int j;
  double wmax,thresh;
  yl_cblas_nr_svdcmp(u,ndata,ma,w,v);
  wmax=0.0;
  for (j=0;j<ma;j++)
    if (w[j] > wmax) wmax=w[j];
  thresh=TOL*wmax;
  for (j=0;j<ma;j++)
    if (w[j] < thresh) w[j]=0.0;
  yl_cblas_nr_svbksb(u,w,v,ndata,ma,b,a);
}
#undef TOL



