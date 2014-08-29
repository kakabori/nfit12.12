#include "mydll.h"
#include <math.h>
#include <stdio.h>
/*    y: da*x + y  */

void yl_blas_daxpy(int n,
		   double da,
		   double *dx,
		   int incx,
		   double *dy,
		   int incy)
{
  double *yend;
  if (da == 0.) return;

//#define YLFORWARD
#ifdef YLFORWARD
  if( incy < 0 ){
	  dy += (n-1)*incy;
	  incy = -incy;
	  dx += (n-1)*incx;
	  incx = -incx;
  }
#endif //YLFORWARD

//#define YLSPECIAL // this does not help
#ifdef YLSPECIAL
  if ( incx == 1 && incy == 1) {
	  yend = dy+n;
	  if (n % 2) *(dy++) += da * *(dx++);
	  while (dy != yend) {
		*(dy++) += da * *(dx++);
		*(dy++) += da * *(dx++);
	  }
	  return;
  }
  /* code for unequal increments or equal increments not equal to 1 */
#endif //YLSPECIAL

  if ( incy == 0 ){
	  *dy += da * *(dx+=incx);
	  return;
  }

  yend = dy + n*incy;
  if(n%2) { *dy += da * *dx; dy += incy; dx += incx;}
  while (dy != yend) {
    *dy += da * *dx; dy += incy; dx += incx;
    *dy += da * *dx; dy += incy; dx += incx;
  }
  return;
}


void yl_blas_saxpy_(int n,
		   float da,
		   float *dx,
		   int incx,
		   float *dy,
		   int incy)
{
  float *yend;
  if (da == 0.) return;

//#define YLFORWARD
#ifdef YLFORWARD
  if( incy < 0 ){
	  dy += (n-1)*incy;
	  incy = -incy;
	  dx += (n-1)*incx;
	  incx = -incx;
  }
#endif //YLFORWARD

//#define YLSPECIAL // this does not help
#ifdef YLSPECIAL
  if ( incx == 1 && incy == 1) {
	  yend = dy+n;
	  if (n % 2) *(dy++) += da * *(dx++);
	  while (dy != yend) {
		*(dy++) += da * *(dx++);
		*(dy++) += da * *(dx++);
	  }
	  return;
  }
  /* code for unequal increments or equal increments not equal to 1 */
#endif //YLSPECIAL

  if ( incy == 0 ){
	  *dy += da * *(dx+=incx);
	  return;
  }

  yend = dy + n*incy;
  if(n%2) { *dy += da * *dx; dy += incy; dx += incx;}
  while (dy != yend) {
    *dy += da * *dx; dy += incy; dx += incx;
    *dy += da * *dx; dy += incy; dx += incx;
  }
  return;
}
