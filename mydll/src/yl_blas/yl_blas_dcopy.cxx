#include "mydll.h"
#include <math.h>
#include <stdio.h>
/*********************************************/
/*     copies a vector, x, to a vector, y.   */
/*********************************************/
void yl_blas_dcopy(int n, double *dx, int incx, double *dy, int incy)
{
  double *yend;

//#define YLFORWARD
#ifdef YLFORWARD
  if( incy < 0 ){
	  dy += (n-1)*incy;
	  incy = -incy;
	  dx += (n-1)*incx;
	  incx = -incx;
  }
#endif //YLFORWARD

//#define YLSPECIAL
#ifdef YLSPECIAL
  if (incx == 1 && incy == 1){
	  yend = dy + n;
	  if(n%2) *(dy++) = *(dx++);
	  while(dy != yend){
		  *(dy++) = *(dx++);
		  *(dy++) = *(dx++);
	  }
	  //memcpy(dy, dx, sizeof(double)*n);
	  // this is not fast as I expected, overhead or bad written?
      return;
  }
#endif //YLSPECIAL
  /*
  if ( incy == 0 ){
	  *dy = *dx;
	  return;
  }
  */
  yend = dy + n*incy;

  if(n%2) { *dy = *dx; dy+=incy; dx+=incx;}
  if ( incx == 0 ){
	  while(dy != yend){
		  *dy = *dx; dy+=incy;
		  *dy = *dx; dy+=incy;
	  }
  } else {
	  while (dy != yend){
		  *dy = *dx; dy+=incy; dx+=incx;
		  *dy = *dx; dy+=incy; dx+=incx;
	  }
  }
}


void yl_blas_scopy_(int n, float *dx, int incx, float *dy, int incy)
{
  float *yend;

//#define YLFORWARD
#ifdef YLFORWARD
  if( incy < 0 ){
	  dy += (n-1)*incy;
	  incy = -incy;
	  dx += (n-1)*incx;
	  incx = -incx;
  }
#endif //YLFORWARD

//#define YLSPECIAL
#ifdef YLSPECIAL
  if (incx == 1 && incy == 1){
	  yend = dy + n;
	  if(n%2) *(dy++) = *(dx++);
	  while(dy != yend){
		  *(dy++) = *(dx++);
		  *(dy++) = *(dx++);
	  }
	  //memcpy(dy, dx, sizeof(double)*n);
	  // this is not fast as I expected, overhead or bad written?
      return;
  }
#endif //YLSPECIAL
  /*
  if ( incy == 0 ){
	  *dy = *dx;
	  return;
  }
  */

  yend = dy + n*incy;

  if(n%2) { *dy = *dx; dy+=incy; dx+=incx;}
  if ( incx == 0 ){
	  while(dy != yend){
		  *dy = *dx; dy+=incy;
		  *dy = *dx; dy+=incy;
	  }
  } else {
	  while (dy != yend){
		  *dy = *dx; dy+=incy; dx+=incx;
		  *dy = *dx; dy+=incy; dx+=incx;
	  }
  }
}



void yl_blas_dset(int n, double da, double *dx, int incx)
{
  double *xend;
  xend = dx + n*incx;
  if(n%2) {*dx = da; dx += incx;}
  while (dx != xend) {
    *dx = da; dx += incx;
    *dx = da; dx += incx;
  }
}
