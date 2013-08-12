CXX = g++
CXXFLAGS = -Wall -fPIC -I. -I../inc -O3 -ffast-math

MISCOBJ =	qromb.o \
		tools.o \
		tvLmdif.o \
		simplex.o \
		polyinterp.o \
		hankelintegrator.o \
		fastdct.o \
		aux.o \
		covar.o \
		linfit.o \
		matcher.o

BLASOBJ = 	yl_blas_dscal.o \
		yl_blas_dnrm2.o \
		yl_blas_dcopy.o \
		yl_blas_dasum.o \
		yl_blas_daxpy.o \
		yl_blas_daxty.o \
		yl_blas_dist.o \
		yl_blas_ddot.o \
		yl_blas_dswap.o \
		yl_cblas_dgemm.o \
		yl_blas_dgemm.o \
		yl_blas_dgemv.o \
		yl_blas_drot.o \
		yl_blas_dger.o \
		yl_blas_dgetf2.o \
		yl_blas_print.o \
		yl_cblas_nr_ludcmp.o \
		yl_cblas_nr_svd.o \
		yl_blas_dm2eop.o \
		yl_cblas_nr_jacobi.o \
		yl_blas_nr_eig.o \
		cmatrix.o


SUBPLEXOBJ =	yl_subplex_calcc.o \
		yl_subplex_simplx.o \
		yl_subplex_subplx.o \
		yl_subplex_partx.o

MINIMIZATIONOBJ =	dfpmin.o \
		yl_minimization.o 

STATOBJ =	uniform.o \
		normal.o \
		statmisc.o


OBJ = $(MISCOBJ) $(BLASOBJ) $(SUBPLEXOBJ) $(MINIMIZATIONOBJ) $(STATOBJ)
