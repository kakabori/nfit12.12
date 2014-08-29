#include "mydll.h"
#include <math.h>
#include <stdio.h>


/* y += alpha * x */
void yl_blas_dadd1(int n, double alpha, double *x, int incx, double *y, int incy){
  double *yend = y + n;
  if(n%2) { *y += alpha * *x; y += incy; x += incx; }
  while(y != yend) {
    *y += alpha * *x; y += incy; x += incx;
    *y += alpha * *x; y += incy; x += incx;
  }
}


/*     A := alpha*x*y' + A, */
// fblas: m row #, n col #
void yl_blas_dger(int m, int n, double alpha,
		  double *x, int incx, double *y, int incy,
		  double *a, int lda, int inca)
{

  if (alpha == 0.)  return;

  if( inca > lda ){
    tswap( lda, inca );
    tswap( x, y);
    tswap( incx, incy);
    tswap( m, n );
  }

  double *ta1, *ta2;
  int minca = m * inca;
  ta2 = a + n * lda;
  lda -= minca;

  while( a != ta2 ){
    double temp;
    ta1 = a + minca;
    temp = alpha * *y;
    while( a != ta1 ){
      *a += *x * temp;
      a += inca;
      x += incx;
    }
    a += lda;
    x -= incx*m;
    y += incy;
  }

}

