#include "mydll.h"
#include <cmath>
#include <cstdlib>
#include <cstdio>


#define NMAX 5000
#define ALPHA 1.0
#define BETA 0.5
#define GAMMA 2.0
void free_vector(double *v,int nl,int nh);

double *nrvector(int nl,int nh){
  double *v;
  v=(double *)malloc((unsigned) (nh-nl+1)*sizeof(double));
  if (!v) printf("allocation failure in nrvector ()");
  return v-nl;
}


#define GET_PSUM for (j=1;j<=ndim;j++) { for (i=1,sum=0.0;i<=mpts;i++) sum += p[i][j]; psum[j]=sum;}

double amotry(double **p,double *y,double *psum,int ndim, double (*funk)(double *),int ihi, int *nfunk, double fac)
{
  int j;
  double fac1,fac2,ytry,*ptry;

  ptry=nrvector(1,ndim);
  fac1=(1.0-fac)/ndim;
  fac2=fac1-fac;
  for (j=1;j<=ndim;j++) ptry[j]=psum[j]*fac1-p[ihi][j]*fac2;
  ytry=(*funk)(ptry);
  ++(*nfunk);
  if (ytry < y[ihi]) {
    y[ihi]=ytry;
    for (j=1;j<=ndim;j++) {
      psum[j] += ptry[j]-p[ihi][j];
      p[ihi][j]=ptry[j];
    }
  }
  free_vector(ptry,1,ndim);
  return ytry;
}


void amoeba(double **p, double *y,int ndim, double ftol, double (*funk)(double *), int *nfunk)
{
  int i,j,ilo,ihi,inhi,mpts=ndim+1;
  double ytry,ysave,sum,rtol,*psum;

  psum=nrvector(1,ndim);
  *nfunk=0;
  GET_PSUM
    for (;;) {
      ilo=1;
      ihi = y[1]>y[2] ? (inhi=2,1) : (inhi=1,2);
      for (i=1;i<=mpts;i++) {
	if (y[i] < y[ilo]) ilo=i;
	if (y[i] > y[ihi]) {
	  inhi=ihi;
	  ihi=i;
	} else if (y[i] > y[inhi])
	  if (i != ihi) inhi=i;
      }
      rtol=2.0*fabs(y[ihi]-y[ilo])/(fabs(y[ihi])+fabs(y[ilo]));
      if (rtol < ftol) break;
      if (*nfunk >= NMAX) printf("Too many iterations in AMOEBA");
      ytry=amotry(p,y,psum,ndim,funk,ihi,nfunk,-ALPHA);
      if (ytry <= y[ilo])
	ytry=amotry(p,y,psum,ndim,funk,ihi,nfunk,GAMMA);
      else if (ytry >= y[inhi]) {
	ysave=y[ihi];
	ytry=amotry(p,y,psum,ndim,funk,ihi,nfunk,BETA);
	if (ytry >= ysave) {
	  for (i=1;i<=mpts;i++) {
	    if (i != ilo) {
	      for (j=1;j<=ndim;j++) {
		psum[j]=0.5*(p[i][j]+p[ilo][j]);
		p[i][j]=psum[j];
	      }
	      y[i]=(*funk)(psum);
	    }
	  }
	  *nfunk += ndim;
	  GET_PSUM
	    }
      }
    }
  free_vector(psum,1,ndim);
}

#undef ALPHA
#undef BETA
#undef GAMMA
#undef NMAX
