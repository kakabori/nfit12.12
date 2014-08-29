#include "mydll.h"
#include <math.h>
#include <stdio.h>
#include <assert.h>

#define a_ptr(a_ld, a_inc) (a + lda * (a_ld) + inca * (a_inc))
#define b_ptr(b_inc) ( b + (b_inc) * incb )
#define TINY 1e-30

void yl_cblas_nr_ludcmp(double *a, int lda, int inca,
			int n, int * indx, double *d__){
  int i, imax=0, j;
  double big, dum;
  int d = 1;
  double *vv = (double *)malloc( sizeof(double) * n ) ;
  
  for( i=0; i<n; i++){
    yl_blas_idamax2( n, a+i*lda, inca, &big );
    assert( big != 0.0);
    vv[i] = 1 / big;
  }
    
  for( j=0; j<n; j++ ){
    double *ajj = a + (lda + inca) * j;

    for(i=1; i<j; i++)
      *a_ptr(i,j) -= yl_blas_ddot(i, a+i*lda, inca, a+j*inca, lda);
    big = 0;
    for( i=j; i<n; i++ ){
      dum = ( vv[i] * 
	      fabs (*a_ptr(i,j) -= yl_blas_ddot(j, a+i*lda, inca, a+j*inca, lda) ) );
      if( dum >= big ){
	big = dum;
	imax = i;
      }
    }
    if( j != imax ) {
      yl_blas_dswap(n, a+imax*lda, inca, a+j*lda, inca);
      d = -d;
      vv[imax] = vv[j];
    }
    indx[j] = imax;
    if ( *ajj == 0.0 ) *ajj = TINY;
    if( j != n-1 )
      yl_blas_dscal(n-j-1, 1.0 / *ajj, ajj+lda, lda );
  }
  *d__ = d;
  free(vv);
}


void yl_cblas_nr_lubksb( double *a, int lda, int inca,
		     int n, int *indx, double *b, int incb){
  int i, ii=-1;
  double sum;
  for(i=0; i<n; i++){
    double *bip = b_ptr( indx[i] );
    sum = *bip;
    *bip = *b_ptr(i);
    if(ii!=-1) sum -= yl_blas_ddot( i - ii, a_ptr(i, ii), inca, b_ptr(ii), incb );
    else if( sum != 0 ) ii = i;
    *b_ptr(i) = sum;
  }
  for( i=n-1; i>=0; i--){
    double *tb = b_ptr(i);
    *tb -= yl_blas_ddot(n-i-1, a_ptr(i, i+1), inca, tb+incb, incb);
    *tb = *b_ptr(i) / *a_ptr(i,i);
  }
}


void yl_cblas_nr_minverse(double *a, int lda, int inca,
			  double *c, int ldc, int incc,
			  int n){
  int *indx = (int *)malloc(sizeof(int) * n);
  double d;
  int j;
  yl_cblas_nr_ludcmp(a, lda, inca, n, indx, &d);
  for(j=0; j<n; j++){
    yl_blas_dset(n, 0, c+j*incc, ldc);
    *(c + j * (ldc + incc) ) = 1;
    yl_cblas_nr_lubksb(a, lda, inca, n, indx, c+j*incc, ldc);
  }
  free(indx);
}


double yl_blas_dpro(int n, double *dx, int incx){
  double *xend;
  double dtemp;
  if(incx == 0)	  return pow(*dx, n);

  xend = dx + n*incx;
  if(n % 2) {dtemp = *dx; dx += incx;}
  else dtemp = 1;
  while (dx != xend) {
    dtemp *= *dx;	dx += incx;
    dtemp *= *dx;	dx += incx;
  }
  return dtemp;
}

double yl_blas_dlogpro(int n, double *dx, int incx){
  double *xend;
  double dtemp;
  if(incx == 0)	  return log(fabs(*dx))*n;

  xend = dx + n*incx;
  if(n % 2) {dtemp = log(fabs(*dx)); dx += incx;}
  else dtemp = 0;
  while (dx != xend) {
    dtemp += log(fabs(*dx));	dx += incx;
    dtemp += log(fabs(*dx));	dx += incx;
  }
  return dtemp;
}

double yl_cblas_nr_mdeterminant(int n, double *a, int lda, int inca){
  int *indx = (int *)malloc(sizeof(int) * n);
  double d;
  yl_cblas_nr_ludcmp(a, lda, inca, n, indx, &d);
  d *= yl_blas_dpro(n, a, lda+inca);
  free(indx);
  return d;
}

double yl_cblas_nr_mlogdet(int n, double *a, int lda, int inca){
  int *indx = (int *)malloc(sizeof(int) * n);
  double d;
  yl_cblas_nr_ludcmp(a, lda, inca, n, indx, &d);
  d = d_sign( yl_blas_dlogpro(n, a, lda+inca), d);
  free(indx);
  return d;
}

#undef a_ptr
#undef TINY
