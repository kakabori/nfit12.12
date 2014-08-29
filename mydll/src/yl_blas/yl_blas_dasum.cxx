#include "mydll.h"
#include <math.h>

/* takes the sum of the absolute values. */

double yl_blas_dasum(int n, double *dx, int incx)
{
  double *xend;
  double dtemp = 0;

//#define YLFORWARD
#ifdef YLFORWARD
  if (incx < 0) {
	dx += (n-1)*incx;
	incx = -incx;
  }
#endif //YLFORWARD

//#define YLSPECIAL
#ifdef YLSPECIAL
  if ( incx == 1 ) {
	  xend = dx+n;
				 // was  n % 6 (by YL) mysterously, 2 is better than 6 on the P3 400MHz, cache size?
				 // I should try it on another type of OS or machine
				 // pointer is better than []
	  if(n % 2) dtemp += fabs(*(dx++));
	  while (dx != xend) {
		dtemp += fabs(*(dx++));
		dtemp += fabs(*(dx++));
	  }
	  return dtemp;
  }
#endif //YLSPECIAL

  if(incx == 0)	  return n*fabs(*dx);

  xend = dx + n*incx;
  if(n % 2) {dtemp += fabs(*dx); dx += incx;}
  while (dx != xend) {
    dtemp += fabs(*dx);	dx += incx;
    dtemp += fabs(*dx);	dx += incx;
  }
  return dtemp;
}

float yl_blas_sasum_(int n, float *dx, int incx)
{
  float *xend;
  float dtemp = 0;

//#define YLFORWARD
#ifdef YLFORWARD
  if (incx < 0) {
	dx += (n-1)*incx;
	incx = -incx;
  }
#endif //YLFORWARD

//#define YLSPECIAL
#ifdef YLSPECIAL
  if ( incx == 1 ) {
	  xend = dx+n;
				 // was  n % 6 (by YL) mysterously, 2 is better than 6 on the P3 400MHz, cache size?
				 // I should try it on another type of OS or machine
				 // pointer is better than []
	  if(n % 2) dtemp += fabs(*(dx++));
	  while (dx != xend) {
		dtemp += fabs(*(dx++));
		dtemp += fabs(*(dx++));
	  }
	  return dtemp;
  }
#endif //YLSPECIAL

  if(incx == 0)	  return n*fabs(*dx);

  xend = dx + n*incx;
  if(n % 2) {dtemp += fabs(*dx); dx += incx;}
  while (dx != xend) {
    dtemp += fabs(*dx);	dx += incx;
    dtemp += fabs(*dx);	dx += incx;
  }
  return dtemp;
}



double yl_blas_dsum(int n, double *dx, int incx)
{
  double *xend;
  double dtemp = 0;

  xend = dx + n*incx;
  if(n % 2) { dtemp += *dx; dx += incx;}
  while (dx != xend) {
    dtemp += *dx ;	dx += incx;
    dtemp += *dx;	dx += incx;
  }
  return dtemp;
}
