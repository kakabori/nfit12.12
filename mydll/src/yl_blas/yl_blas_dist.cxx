#include "mydll.h"
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include "yl_simpletools.h"
/* dist calculates the distance between the points x,y. */

double yl_blas_dist(int n, double *dx, int incx, double *dy, int incy)
{
    double d1, scale=0, absxmy, sum=1;
    double *yend;

	if( incy == 0 ){
		if( incx == 0 ) return sqrt((double)n)*fabs(*dx - *dy);
		else {
			yl_swap(incx, incy);
			yl_swap(dx, dy);
		}
	}
//#define YLFORWARD
#ifdef YLFORWARD
	if( incy < 0 ){
		dy += (n-1)*incy;
		incy = -incy;
		dx += (n-1)*incx;
		incx = -incx;
	}
#endif //YLFORWARD

	yend = dy + n*incy;
    while ( dy != yend ) {
		if( ( absxmy = fabs(*dx - *dy) ) != 0. ){
			if (absxmy > scale) {
				d1 = scale / absxmy;
				sum = sum * (d1 * d1) + 1.;
				scale = absxmy;
			} else {
				d1 = absxmy / scale;
				sum += d1 * d1;
			}
		}
		dx += incx;
		dy += incy;
    }
    return scale * sqrt(sum);
}
/*
the following sounds good but not in pratice, why?

	bool ifswap = false;

	if( abs(incy) > abs(incx) ) ifswap = (incx == 0);
	else if( incy == 0 ){
		if( incx == 0 ) return sqrt(n)*fabs(*dx - *dy);
		else ifswap = true;
	}
	if( ifswap ){
		yl_swap(incx, incy);
		yl_swap(dx, dy);
	}

*/

double yl_blas_sist(int n, float *dx, int incx, float *dy, int incy)
{
    float d1, scale=0, absxmy, sum=1;
    float *yend;

	if( incy == 0 ){
		if( incx == 0 ) return sqrt((double)n)*fabs(*dx - *dy);
		else {
			yl_swap(incx, incy);
			yl_swap(dx, dy);
		}
	}
//#define YLFORWARD
#ifdef YLFORWARD
	if( incy < 0 ){
		dy += (n-1)*incy;
		incy = -incy;
		dx += (n-1)*incx;
		incx = -incx;
	}
#endif //YLFORWARD

	yend = dy + n*incy;
    while ( dy != yend ) {
		if( ( absxmy = fabs(*dx - *dy) ) != 0. ){
			if (absxmy > scale) {
				d1 = scale / absxmy;
				sum = sum * (d1 * d1) + 1.;
				scale = absxmy;
			} else {
				d1 = absxmy / scale;
				sum += d1 * d1;
			}
		}
		dx += incx;
		dy += incy;
    }
    return scale * sqrt(sum);
}
