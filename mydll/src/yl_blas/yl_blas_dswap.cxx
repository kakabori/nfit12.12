#include "mydll.h"
#include <math.h>
#include <stdio.h>

void yl_blas_dswap(int n, double *dx, int incx, double *dy, int incy)
{
  double *yend;
  double tmp;
  /*
  if ( incy == 0 ){
	  *dy = *dx;
	  return;
  }
  */
  yend = dy + n*incy;

  if(n%2) { tmp =*dy;  *dy = *dx; *dx = tmp; dy+=incy; dx+=incx;}
  while(dy != yend){
    tmp =*dy;  *dy = *dx; *dx = tmp; dy+=incy; dx+=incx;
    tmp =*dy;  *dy = *dx; *dx = tmp; dy+=incy; dx+=incx;
  }
}

