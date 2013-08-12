inline double yl_blas_inline_ddot(int n, double *dx, int incx, double *dy, int incy){
  double *yend;
  double sum=0;
  yend = dy + n*incy;
  if(n%2) { sum += *dy * *dx; dy += incy; dx += incx;}
  while (dy != yend) {
    sum += *dy * *dx; dy += incy; dx += incx;
    sum += *dy * *dx; dy += incy; dx += incx;
  }
  return sum;
}

inline double yl_blas_inline_ddot_s(int n, double *dx, double *dy){
  double *yend;
  double sum=0;
  yend = dy + n;
  if(n%2) { sum += *dy * *dx; dy++; dx++;}
  while (dy != yend) {
    sum += *dy * *dx; dy++; dx++;
    sum += *dy * *dx; dy++; dx++;
  }
  return sum;
}

inline void yl_blas_inline_daxpy_(int n, double da, double *dx, int incx, double *dy, int incy){
  double *yend;
  if (da == 0.) return;
  yend = dy + n*incy;
  if(n%2) { *dy += da * *dx; dy += incy; dx += incx;}
  while (dy != yend) {
    *dy += da * *dx; dy += incy; dx += incx;
    *dy += da * *dx; dy += incy; dx += incx;
  }
  return;
}

