#include "mydll.h"
#include <math.h>
/*  DNRM2 := sqrt( x'*x )  */ 

double yl_blas_dnrm2(int n, double *dx, int incx)
{
    double d1, scale=0, absx, sum=1;
    double *xend;

	if( incx == 0 ) return fabs(*dx);

//#define YLFORWARD
#ifdef YLFORWARD
	if( incx < 0 ){
		dx += (n-1)*incx;
		incx = -incx;
	}
#endif //YLFORWARD

	xend = dx + n*incx;
    while ( dx != xend ) {
		if( ( absx = fabs(*dx) ) != 0. ){
			if (absx > scale) {
				d1 = scale / absx;
				sum = sum * (d1 * d1) + 1.;
				scale = absx;
			} else {
				d1 = absx / scale;
				sum += d1 * d1;
			}
		}
		dx += incx;
    }
    return scale * sqrt(sum);
}

float yl_blas_snrm2(int n, float *dx, int incx)
{
    float d1, scale=0, absx, sum=1;
    float *xend;

	if( incx == 0 ) return fabs(*dx);

//#define YLFORWARD
#ifdef YLFORWARD
	if( incx < 0 ){
		dx += (n-1)*incx;
		incx = -incx;
	}
#endif //YLFORWARD

	xend = dx + n*incx;
    while ( dx != xend ) {
		if( ( absx = fabs(*dx) ) != 0. ){
			if (absx > scale) {
				d1 = scale / absx;
				sum = sum * (d1 * d1) + 1.;
				scale = absx;
			} else {
				d1 = absx / scale;
				sum += d1 * d1;
			}
		}
		dx += incx;
    }
    return scale * sqrt(sum);
}

