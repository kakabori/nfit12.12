#include "mydll.h"
#include <math.h>
#include <stdio.h>
#include "yl_simpletools.h"

double  yl_blas_ddot(int n,
		   double *dx,
		   int incx,
		   double *dy,
		   int incy)
{
  double *yend;
  double sum=0;

  if( incy == 0 ){
    if( incx == 0 ) return n * *dx * *dy;
    else {
      yl_swap(incx, incy);
      yl_swap(dx, dy);
    }
  }

  yend = dy + n*incy;
  if(n%2) { sum += *dy * *dx; dy += incy; dx += incx;}
  while (dy != yend) {
    sum += *dy * *dx; dy += incy; dx += incx;
    sum += *dy * *dx; dy += incy; dx += incx;
  }
  return sum;
}

double  yl_blas_dnrm22(int n,
		   double *dy,
		   int incy)
{
  double *yend;
  double sum=0;

  if( incy == 0 ) return n * *dy * *dy;

  yend = dy + n*incy;
  if(n%2) { sum += *dy * *dy; dy += incy; }
  while (dy != yend) {
    sum += *dy * *dy; dy += incy;
    sum += *dy * *dy; dy += incy;
  }
  return sum;
}

