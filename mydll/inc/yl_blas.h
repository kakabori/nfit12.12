#ifndef _YL_BLAS_HEADER_
#define _YL_BLAS_HEADER_


MYDLL_API double yl_blas_dasum(int n, double *dx, int incx);
MYDLL_API double yl_blas_dsum(int n, double *dx, int incx);
MYDLL_API float  yl_blas_sasum(int n, float  *dx, int incx);

MYDLL_API void yl_blas_dcopy(int n, double *dx, int incx, double *dy, int incy);
MYDLL_API void yl_blas_scopy(int n, float  *dx, int incx, float  *dy, int incy);

MYDLL_API double yl_blas_dist(int n, double *dx, int incx, double *dy, int incy);
MYDLL_API double yl_blas_sist(int n, float  *dx, int incx, float  *dy, int incy);

MYDLL_API double yl_blas_dnrm2(int n, double *dx, int incx);
MYDLL_API double yl_blas_dnrm22(int n, double *dx, int incx);
MYDLL_API float  yl_blas_snrm2(int n, float  *dx, int incx);

MYDLL_API void yl_blas_dscal(int n, double da, double *dx, int incx);
MYDLL_API void yl_blas_sscal(int n, float  da, float  *dx, int incx);

MYDLL_API void yl_blas_daxpy(int n, double da, double *dx, int incx, double *dy, int incy);
MYDLL_API void yl_blas_saxpy(int n, float  da, float  *dx, int incx, float  *dy, int incy);

MYDLL_API void yl_blas_daxty(int n, double da, double *dx, int incx, double *dy, int incy);
MYDLL_API void yl_blas_saxty(int n, float  da, float  *dx, int incx, float  *dy, int incy);

MYDLL_API double yl_blas_ddot(int n, double *dx, int incx, double *dy, int incy);

MYDLL_API void yl_blas_dset(int n, double da, double *dx, int incx);

MYDLL_API void yl_blas_dmset(double a, double* x, int ldx, int incx, int n, int mincx);

MYDLL_API void yl_blas_dmscal(double a, double* x, int ldx, int incx, int n, int mincx);

MYDLL_API void yl_blas_dmcopy(double* x, int ldx, int incx, double* y, int ldy, int incy, int n, int m);
MYDLL_API void yl_blas_dger(int m, int n, double alpha, double *x, int incx, double *y, int incy, double *a, int lda, int inca);


MYDLL_API void yl_blas_dswap(int n, double *dx, int incx, double *dy, int incy);

MYDLL_API void yl_blas_dme2e( double* x, int ldx, int incx, double* y, int ldy, int incy, int n, int m, E2eOp op);
MYDLL_API void yl_blas_dmxpytc( double* x, int ldx, int incx, double* y, int ldy, int incy, double* c, int ldc, int incc, int n, int m );
MYDLL_API void yl_blas_dmxmytc( double* x, int ldx, int incx, double* y, int ldy, int incy, double* c, int ldc, int incc, int n, int m );
MYDLL_API void yl_blas_dmaxpbytc( double* x, int ldx, int incx, double* y, int ldy, int incy, double* c, int ldc, int incc, int n, int m, double alpha, double beta);
MYDLL_API void yl_blas_dmxpbytc( double* x, int ldx, int incx, double* y, int ldy, int incy, double* c, int ldc, int incc, int n, int m, double beta);
MYDLL_API void yl_blas_dmxytc( double* x, int ldx, int incx, double* y, int ldy, int incy, double* c, int ldc, int incc, int n, int m);
MYDLL_API void yl_blas_dmaxpy( double* x, int ldx, int incx, double* y, int ldy, int incy, int n, int m, double alpha);

MYDLL_API void yl_blas_mprint(ostream& out, int m, int n, double *a, int lda, int inca, int trans);
MYDLL_API void yl_blas_vprint(ostream& out, int n, double *a, int inca);

MYDLL_API int yl_blas_idamax(int n, double* x, int incx);
MYDLL_API int yl_blas_idamax2(int n, double* x, int incx, double *max__);

MYDLL_API int yl_blas_dgetf2(int m, int n, double *a, int lda, int inca, int *ipiv);
MYDLL_API void yl_cblas_nr_ludcmp(double *a, int lda, int inca, int n, int * indx, double *d__);
MYDLL_API void yl_cblas_nr_lubksb( double *a, int lad, int inca, int n, int *indx, double *b, int incb);
MYDLL_API void yl_cblas_nr_minverse(double *a, int lda, int inca, double *c, int ldc, int incc, int n);
MYDLL_API void yl_cblas_nr_svbksb(double *u,double *w,double *v,int m,int n,double *b,double *x);
MYDLL_API void yl_cblas_nr_svdcmp(double *a,int m,int n,double *w,double *v);
MYDLL_API void yl_cblas_nr_svdcmp( int m, int n,
			 double *a, int lda, int inca,
			 double *w, int incw,
			 double *v, int ldv, int incv);
MYDLL_API void yl_cblas_nr_svbksb(int m,int n,
			double *u, int ldu, int incu,
			double *w, int incw,
			double *v, int ldv, int incv,
			double *b, int incb,
			double *x, int incx);


MYDLL_API void yl_cblas_nr_jacobi(int n, double *a, int lda, int inca, double *v, int ldv, int incv, double *d, int *nrot);
MYDLL_API double yl_cblas_nr_mdeterminant(int n, double *a, int lda, int inca);
MYDLL_API double yl_cblas_nr_mlogdet(int n, double *a, int lda, int inca);

MYDLL_API void yl_blas_nr_tred2(int n, double *a, int lda, int inca, double *d, double *e);
MYDLL_API void yl_blas_nr_tqli(int n, double *z, int ld, int inc, double *d, double *e );
MYDLL_API void yl_blas_nr_elmhes(int n, double *a, int lda, int inca);
MYDLL_API void yl_blas_nr_hqr(int n, double *a, int lda, int inca, double *wr, double *wi);


MYDLL_API void yl_blas_dgemm_n_ijk(double *a, int lda, int inca, double *b, int ldb, int incb, double *c, int ldc, int incc, int ni, int nk, int nj, char);

MYDLL_API void yl_blas_dgemm_np_jik(double *a, int lda, int inca,double *b, int ldb, int incb, double *c, int ldc, int incc, int ni, int nk, int nj);

MYDLL_API void yl_blas_dgemm_np_jik_opt(double *a, int lda, int inca, double *b, int ldb, int incb, double *c, int ldc, int incc, int ni, int nk, int nj, char);




MYDLL_API void yl_cblas_dgemm2(double *a, int lda, int inca, double *b, int ldb, int incb, double *c, int ldc, int incc, int ni, int nk, int nj);

MYDLL_API void yl_cblas_dgemm(double *a, int lda, int inca, double *b, int ldb, int incb, double *c, int ldc, int incc, int ni, int nk, int nj, double, double);

MYDLL_API void yl_cblas_dgemm_opt(double *a, int lda, int inca, double *b, int ldb, int incb, double *c, int ldc, int incc, int ni, int nk, int nj, double, double, char);

MYDLL_API void yl_blas_dgemm_auto(double *a, int lda, int inca, double *b, int ldb, int incb, double *c, int ldc, int incc, int ni, int nk, int nj, double, double, char);

MYDLL_API void yl_cblas_dgemm_np_ijk(double *a, int lda, int inca,double *b, int ldb, int incb,	double *c, int ldc, int incc, int ni, int nk, int nj);

MYDLL_API void yl_cblas_dgemm_np_kij(double *a, int lda, int inca, double *b, int ldb, int incb, double *c, int ldc, int incc, int ni, int nk, int nj);

MYDLL_API void yl_cblas_dgemm_n_ijk(double *a, int lda, int inca, double *b, int ldb, int incb,	double *c, int ldc, int incc, int ni, int nk, int nj, char);

MYDLL_API void yl_blas_dgemv_auto(double *a, int lda, int inca, double *b, int incb, double *c, int incc, int ni, int nk, double, double, char);

typedef void (*YL_CBLAS_DMSCAL)(double a,
				double* x, int ldx, int incx,
				int n, int mincx);

typedef void (*YL_CBLAS_DSCAL)(int n, double a, double* x, int incx);

typedef void (*YL_CBLAS_DGEMM_NP)(double *a, int lda, int inca,
				  double *b, int ldb, int incb,
				  double *c, int ldc, int incc,
				  int ni, int nk, int nj);
typedef void (*YL_CBLAS_DGEMM_NPA)(double *a, int lda, int inca,
				  double *b, int ldb, int incb,
				  double *c, int ldc, int incc,
				  int ni, int nk, int nj, double);

typedef void (*GEVFUNC)();

#endif
