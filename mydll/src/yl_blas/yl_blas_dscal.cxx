#include "mydll.h"
#include <math.h>
#include <stdio.h>
/*     scales a vector by a constant.  dx *= da */
void yl_blas_dscal(int n, double da, double *dx, int incx)
{
  double *xend;

//#define YLFORWARD
#ifdef YLFORWARD
  if ( incx < 0 ){
	dx += (n-1)*incx;
	incx = -incx;
  }
#endif //YLFORWARD

//#define YLSPECIAL
#ifdef YLSPECIAL
  if( da == 1. ) return;
  if ( incx == 1 ) {
	  int i, m;
	  xend = dx+n;
	  m = n % 4;
	  for (i = 0; i < m; ++i) *(dx++) *= da;
	  while (dx < xend) {
		*(dx++) *= da;
		*(dx++) *= da;
		*(dx++) *= da;
		*(dx++) *= da;
	  }
	  //printf("h");
	  return;
  }
#endif //YLSPECIAL

  if( incx == 0 ) {
	  *dx *= da;
	  return;
  }
  xend = dx + n*incx;
  if(n%2) {*dx *= da; dx += incx;}
  while (dx != xend) {
    *dx *= da; dx += incx;
    *dx *= da; dx += incx;
  }
}
void yl_blas_sscal_(int n, float da, float *dx, int incx)
{
  float *xend;

//#define YLFORWARD
#ifdef YLFORWARD
  if ( incx < 0 ){
	dx += (n-1)*incx;
	incx = -incx;
  }
#endif //YLFORWARD

//#define YLSPECIAL
#ifdef YLSPECIAL
  if( da == 1. ) return;
  if ( incx == 1 ) {
	  int i, m;
	  xend = dx+n;
	  m = n % 4;
	  for (i = 0; i < m; ++i) *(dx++) *= da;
	  while (dx < xend) {
		*(dx++) *= da;
		*(dx++) *= da;
		*(dx++) *= da;
		*(dx++) *= da;
	  }
	  //printf("h");
	  return;
  }
#endif //YLSPECIAL

  if( incx == 0 ) {
	  *dx *= da;
	  return;
  }
  xend = dx + n*incx;
  if(n%2) {*dx *= da; dx += incx;}
  while (dx != xend) {
    *dx *= da; dx += incx;
    *dx *= da; dx += incx;
  }
}
