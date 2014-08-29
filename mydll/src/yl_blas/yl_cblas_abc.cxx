#include "mydll.h"


#include <iostream>
#define  BIGBZ 1024
#define  SMLBZ 10
void yl_blas_dgemm_auto(double *a, int lda, int inca,
			 double *b, int ldb, int incb,
			 double *c, int ldc, int incc,
			 int ni, int nk, int nj,
			 double alpha, double beta, char opt){
  int bsz, bszk;
  YL_CBLAS_DMSCAL dmsf = (beta==0.) ? yl_blas_dmset : yl_blas_dmscal;
  YL_CBLAS_DGEMM_NP dgemmnp;

  if(alpha==0.){
    dmsf(0, c, ldc, incc, nj, ni*incc);
    return;
  }
  
  bsz = bszk = SMLBZ;
  if(opt & 2)   {bszk = BIGBZ; bsz = SMLBZ; tswap(lda, inca); dgemmnp = yl_blas_dgemm_np_jik; }
  else          {bsz = bszk = SMLBZ; dgemmnp = yl_blas_dgemm_np_ijk; }
  if(opt & 1)   {bszk = SMLBZ; tswap(ldb, incb);}
  cout<<(int)opt<<' '<<bsz<<' '<<bszk<<endl;


  if(! (opt & 2) ){  
    int jj, kk;
    int bszc;
    for(jj=0; jj<nj; jj+=bsz){
      bszc = yfmin ( bsz , nj-jj );
      dmsf(beta, c + jj * ldc, ldc, incc, bszc, ni*incc);
      for(kk=0; kk<nk; kk+=bszk){
	dgemmnp(a + kk * lda, lda, inca,
		b + kk * incb + jj * ldb,  ldb,  incb,
		c + jj * ldc, ldc, incc,
		ni, yfmin ( bszk , nk-kk ), bszc);
      }
    }
  } else if(0){
    int jj, kk;
    int bszc;
    for(kk=0; kk<nk; kk+=bszk){
      for(jj=0; jj<nj; jj+=bsz){
	bszc = yfmin ( bsz , nj-jj );
	if(kk==0) dmsf(beta, c + jj * ldc, ldc, incc, bszc, ni*incc);
	dgemmnp(a + kk * lda, lda, inca,
		b + kk * incb + jj * ldb,  ldb,  incb,
		c + jj * ldc, ldc, incc,
		ni, yfmin ( bszk , nk-kk ), bszc);
      }
    }
  }
  else if(1) {
    int ii, kk;
    int bszc;
    for(ii=0; ii<ni; ii+=bsz){
      bszc = yfmin ( bsz , ni-ii );
      dmsf(beta, c + ii * incc, ldc, incc, nj, bszc*incc);
      for(kk=0; kk<nk; kk+=bszk){
	dgemmnp(a + ii * inca + kk * lda, lda, inca, 
		b + kk * incb, ldb, incb, 
		c + ii * incc, ldc, incc, 
		bszc, yfmin ( bszk , nk-kk ), nj);
      }
    }
  } else if(0){
    int ii, kk;
    int bszc;
    for(kk=0; kk<nk; kk+=bszk){
      for(ii=0; ii<ni; ii+=bsz){
	bszc = yfmin ( bsz , ni-ii );
	if(kk==0) dmsf(beta, c + ii * incc, ldc, incc, nj, bszc*incc);
	dgemmnp(a + ii * inca + kk * lda, lda, inca, 
		b + kk * incb, ldb, incb, 
		c + ii * incc, ldc, incc, 
		bszc, yfmin ( bszk , nk-kk ), nj);
      }
    }
  } else {
    int ii, jj, kk;
    int bszc, bszj=SMLBZ, bsza;
    bszk = 10;
    for(jj=0; jj<nj; jj+=bszj){
      bszc = yfmin ( bszj , nj-jj );
      dmsf(beta, c + jj * ldc, ldc, incc, bszc, ni*incc);
      for(ii=0; ii<ni; ii+=bsz){
	bsza = yfmin ( bsz , ni-ii );
	for(kk=0; kk<nk; kk+=bszk){
	  dgemmnp(a + ii * inca + kk * lda, lda, inca, 
		  b + kk * incb + jj * ldb, ldb, incb, 
		  c + ii * incc + jj * ldc, ldc, incc, 
		  bsza, yfmin ( bszk , nk-kk ), bszc);
	}
      }
    }
  }

}
