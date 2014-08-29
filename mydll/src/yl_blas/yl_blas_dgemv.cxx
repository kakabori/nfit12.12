#include "mydll.h"
#include "yl_cblas_.h"

void yl_blas_dgemv_np_ik(double *a, int lda, int inca,
			  double *b, int incb,
			  double *c, int incc,
			  int ni, int nk){
  double *tc1;
  tc1 = c + ni * incc;
  while( c != tc1 ){
    *c += yl_blas_inline_ddot (nk, a, lda, b, incb );
    c += incc;
    a += inca;
  }
}

void yl_blas_dgemv_npa_ik(double *a, int lda, int inca,
			  double *b, int incb,
			  double *c, int incc,
			  int ni, int nk, double alpha){
  double *tc1;
  tc1 = c + ni * incc;
  while( c != tc1 ){
    *c += alpha * yl_blas_inline_ddot (nk, a, lda, b, incb );
    c += incc;
    a += inca;
  }
}

#include <iostream>
#define  BIGBZ 1024
#define  SMLBZ 10

void yl_blas_dgemv_auto(double *a, int lda, int inca,
			double *b, int incb,
			double *c, int incc,
			int ni, int nk,
			double alpha, double beta, char opt){
  int bsz, bszk;
  YL_CBLAS_DSCAL dmsf = (beta==0.) ? yl_blas_dset : yl_blas_dscal;

  dmsf(ni, beta, c, incc);
  if(alpha==0.)    return;
  
  if(opt & 2) {
    bszk = BIGBZ;
    bsz  = SMLBZ;
    tswap(lda, inca);
  } else {
    bsz =  BIGBZ;
    bszk = SMLBZ;
  }

  int ii, kk;
  int bszc;
  for(ii=0; ii<ni; ii+=bsz){
    bszc = yfmin ( bsz , ni-ii );
    for(kk=0; kk<nk; kk+=bszk){
      if(alpha==1.) yl_blas_dgemv_np_ik(a + ii * inca + kk * lda, lda, inca, 
					b + kk * incb, incb, 
					c + ii * incc, incc, 
					bszc, yfmin ( bszk , nk-kk ));
      else         yl_blas_dgemv_npa_ik(a + ii * inca + kk * lda, lda, inca, 
					b + kk * incb, incb, 
					c + ii * incc, incc, 
					bszc, yfmin ( bszk , nk-kk ), alpha);
    }
  }
    /*   for(int kk=0; kk<nk; kk+=bszk){
    if(alpha==1.) yl_blas_dgemv_np_ik(a + kk * lda, lda, inca, 
    b + kk * incb, incb, 
    c, incc, 
    ni, yfmin ( bszk , nk-kk ));
    else         yl_blas_dgemv_npa_ik(a + kk * lda, lda, inca, 
    b + kk * incb, incb, 
    c, incc, 
    ni, yfmin ( bszk , nk-kk ), alpha);
    }*/
}

