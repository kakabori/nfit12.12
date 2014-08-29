#include "mydll.h"
#include "yl_cblas_.h"


// a : the value the matrix will be set to
// n : number of leading dimension
// m : the secondary dimension
// mincx : the secondary dimension * incx = memsize of one leading sequence
void yl_blas_dmset(double a,
		   double* x, int ldx, int incx,
		   int n, int mincx){
  double *x1, *x2;
  x2 = x + n * ldx;
  n = ldx - mincx;
  while(x != x2){
    x1 = x + mincx;
    while(x != x1){
      *x = a;
      x += incx;
    }
    x += n;
  }
}

void yl_blas_dmscal(double a,
		    double* x, int ldx, int incx,
		    int n, int mincx){
  double *x1, *x2;
  x2 = x + n * ldx;
  n = ldx - mincx;
  while(x != x2){
    x1 = x + mincx;
    while(x != x1){
      *x *= a;
      x += incx;
    }
    x += n;
  }
}


// n is the #  of leading  dimension
void yl_blas_dmcopy( double* x, int ldx, int incx,
		      double* y, int ldy, int incy,
		      int n, int m){
  double *y1, *y2;
  int mincy = m*incy;
  
  y2 = y + n * ldy;
  n = ldx - m * incx;
  m = ldy - mincy;

  while(y != y2){
    y1 = y + mincy;
    while(y != y1){
      *y = *x;
      x += incx;   y += incy;
    }
    x += n; y += m;
  }
}

//**********************************
// n is the #  of leading  dimension
void yl_blas_dmxpytc( double* x, int ldx, int incx,
		      double* y, int ldy, int incy,
		      double* c, int ldc, int incc,
		      int n, int m ){
  double *y1, *y2 = y + n * ldy;
  ldc -= m * incc;
  ldx -= m * incx;
  ldy -= (m *= incy);

  while(y != y2){
    y1 = y + m;
    while(y != y1){
      *c = *x + *y;
      x += incx;
      y += incy;
      c += incc;
    }
    x += ldx; y += ldy; c += ldc;
  }
}



void yl_blas_dmxmytc( double* x, int ldx, int incx,
		       double* y, int ldy, int incy,
		       double* c, int ldc, int incc,
		       int n, int m){
  double *y1, *y2 = y + n * ldy;
  ldc -= m * incc;
  ldx -= m * incx;
  ldy -= (m *= incy);

  while(y != y2){
    y1 = y + m;
    while(y != y1){
      *c = *x - *y;
      x += incx;
      y += incy;
      c += incc;
    }
    x += ldx; y += ldy; c += ldc;
  }
}




// c = alpha * x + beta * y
void yl_blas_dmaxpbytc( double* x, int ldx, int incx,
		       double* y, int ldy, int incy,
		       double* c, int ldc, int incc,
		       int n, int m, double alpha, double beta){
  double *y1, *y2;
  
  y2 = y + n * ldy;
  ldc -= m * incc;
  ldx -= m * incx;
  m *= incy;
  ldy -= m;

  while(y != y2){
    y1 = y + m;
    while(y != y1){
      *c = alpha * *x + beta * *y;
      x += incx;
      y += incy;
      c += incc;
    }
    x += ldx; y += ldy; c += ldc;
  }
}
// c = x + beta * y
void yl_blas_dmxpbytc( double* x, int ldx, int incx,
		       double* y, int ldy, int incy,
		       double* c, int ldc, int incc,
		       int n, int m, double beta){
  double *y1, *y2;
  
  y2 = y + n * ldy;
  ldc -= m * incc;
  ldx -= m * incx;
  m *= incy;
  ldy -= m;

  while(y != y2){
    y1 = y + m;
    while(y != y1){
      *c = *x + *y * beta;
      x += incx;
      y += incy;
      c += incc;
    }
    x += ldx; y += ldy; c += ldc;
  }
}
// y +=  alpha * x
void yl_blas_dmaxpy( double* x, int ldx, int incx,
		     double* y, int ldy, int incy,
		     int n, int m, double alpha){
  double *y1, *y2;
  
  y2 = y + n * ldy;
  ldx -= m * incx;
  m *= incy;
  ldy -= m;

  while(y != y2){
    y1 = y + m;
    while(y != y1){
      *y += alpha * *x;
      x += incx;
      y += incy;
    }
    x += ldx; y += ldy;
  }
}


// c =  x * y
void yl_blas_dmxytc( double* x, int ldx, int incx,
		     double* y, int ldy, int incy,
		     double* c, int ldc, int incc,
		     int n, int m ){
  double *y1, *y2;
  
  y2 = y + n * ldy;
  ldc -= m * incc;
  ldx -= m * incx;
  m *= incy;
  ldy -= m;

  while(y != y2){
    y1 = y + m;
    while(y != y1){
      *c = *x * *y;
      x += incx;
      y += incy;
      c += incc;
    }
    x += ldx; y += ldy; c += ldc;
  }
}



void yl_blas_dme2e( double* x, int ldx, int incx,
		    double* y, int ldy, int incy,
		    int n, int m, E2eOp op){
  double *y1, *y2;
  int mincy = m*incy;
  
  y2 = y + n * ldy;
  n = ldx - m * incx;
  m = ldy - mincy;

  while(y != y2){
    y1 = y + mincy;
    while(y != y1){
      op(x , y);
      x += incx;   y += incy;
    }
    x += n; y += m;
  }
}

