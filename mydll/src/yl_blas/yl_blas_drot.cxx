#include "mydll.h"
#include <math.h>
#include <stdio.h>

void yl_blas_drot(int n,
		  double *dx, int incx,
		  double *dy, int incy,
		  double c, double s)
{
  if( incy == 0 ) {
    if( incx == 0 ) {
      double dtemp;
      dtemp = c * *dx + s * *dy;
      *dy   = c * *dy - s * *dx;
      *dx   = dtemp;
    }
    tswap(dx, dy);
    tswap(incx, incy);
  }
  double *yend = dy + n * incy;
  /*     applies a plane rotation.   */
  while( dy != yend ){
    double dtemp;
      dtemp = c * *dx + s * *dy;
      *dy   = c * *dy - s * *dx;
      *dx   = dtemp;
      dx += incx;
      dy += incy;
  }
}


int yl_blas_drotg(double *_da, double *_db, double *c, double *s)
{
  double da = *_da;
  double db = *_db;
  double d1, d2;

  double r, scale, z, roe;
  double fabsda = fabs(da);
  double fabsdb = fabs(db);

  roe = yfmax( fabsdb, fabsdb );
  scale = fabsda + fabsdb;

  if (scale != 0.){
    d1 = da / scale;
    d2 = db / scale;
    r = scale * sqrt(d1 * d1 + d2 * d2);
    //??? r = ( roe >= 0) ? r : -r;
    *c = da / r;
    *s = db / r;
    z = 1.;
    if (fabsda > fabsdb)    z = *s;
    else if (*c != 0.) z = 1. / *c;
  } else {
    *c = 1.;
    *s = 0.;
    r = 0.;
    z = 0.;
  }
  *_da = r;
  *_db = z;
  return 0;
}
