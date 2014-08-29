#include "mydll.h"
#include <math.h>
#include <stdio.h>

int yl_blas_idamax(int n, double* x, int incx){
  double dmax = fabs( *x );
  double *xmax = x;
  double *xend = x + n * incx;
  
  while( (x+=incx) != xend ){
    double tmp;
    if( ( tmp = fabs(*x) ) > dmax ){
      dmax = tmp;
      xmax = x;
    }
  }
  return n - (xend - xmax) / incx;
}

int yl_blas_idamax2(int n, double* x, int incx, double *max__){
  double dmax = fabs( *x );
  double *xmax = x;
  double *xend = x + n * incx;
  
  while( (x+=incx) != xend ){
    double tmp;
    if( ( tmp = fabs(*x) ) > dmax ){
      dmax = tmp;
      xmax = x;
    }
  }
  *max__ = dmax;
  return n - (xend - xmax) / incx;
}

// L, U, U has unit value on diag
// fblas based
int yl_blas_dgetf2(int m, int n,
		   double *a, int lda, int inca,
		   int *ipiv){
  int info = 0;
  int minmn = yfmin(m,n);


  for(int j=0; j< minmn; j++){
    int jp1 = j + 1;
    double *ajj = a + (lda + inca) * j;
    int jpmj = yl_blas_idamax(m-j, ajj , inca);

    //    yl_blas_mprint(m,n,a,lda,inca,0);

    ipiv[j] = jpmj + j;
    if( *(ajj + jpmj * inca) != 0 ){
      if( jpmj != 0 )
	yl_blas_dswap(n, a+j*inca, lda, a+(jpmj+j)*inca, lda);
      if( jp1 < m )
	yl_blas_dscal(m-jp1, 1 / (*ajj), ajj+inca, inca);
    } else if( info == 0 ){
      info = j;
    }

    if( jp1 < minmn ){
      yl_blas_dger( m-jp1 , n-jp1, -1,
		    ajj+inca,  inca,
		    ajj+lda, lda,
		    ajj+lda+inca, lda, inca);
    }
  }
  return info;
}
