#include "mydll.h"
#include "yl_cblas_.h"

/*C += A * B (ijk)  (_ | .) , interchange ldx vs. incx & a vs. b & ni vs. nj -> ijk */
void yl_blas_dgemm_np_jik(double *a, int lda, int inca,
			  double *b, int ldb, int incb,
			  double *c, int ldc, int incc,
			  int ni, int nk, int nj){
  double *tc1, *tc2;
  
  tc2 = c + nj * ldc;
  ldc -= ni * incc;
  nj = ni * inca;
  ni = ni * incc;

  while( c != tc2 ){
    tc1 = c + ni;
    while( c != tc1 ){
      *c += yl_blas_inline_ddot (nk, a, lda, b, incb );
      c += incc;
      a += inca;
    }
    c += ldc;
    a -= nj;
    b += ldb;
  }
}

void yl_blas_dgemm_np_ijk(double *a, int lda, int inca,
			  double *b, int ldb, int incb,
			  double *c, int ldc, int incc,
			  int ni, int nk, int nj){
  double *tc1, *tc2;
  
  tc2 = c + ni * incc;
  incc -= nj * ldc;
  ni = nj * ldb;
  nj = nj * ldc;

  while( c != tc2 ){
    tc1 = c + nj;
    while( c != tc1 ){
      *c += yl_blas_inline_ddot (nk, b, incb, a, lda );
      c += ldc;
      b += ldb;
    }
    c += incc;
    b -= ni;
    a += inca;
  }
}

void yl_blas_dgemm_npa_jik(double *a, int lda, int inca,
			  double *b, int ldb, int incb,
			  double *c, int ldc, int incc,
			  int ni, int nk, int nj, double alpha){
  double *tc1, *tc2;
  
  tc2 = c + nj * ldc;
  ldc -= ni * incc;
  nj = ni * inca;
  ni = ni * incc;

  while( c != tc2 ){
    tc1 = c + ni;
    while( c != tc1 ){
      *c += alpha * yl_blas_inline_ddot (nk, a, lda, b, incb );
      c += incc;
      a += inca;
    }
    c += ldc;
    a -= nj;
    b += ldb;
  }
}

void yl_blas_dgemm_npa_ijk(double *a, int lda, int inca,
			  double *b, int ldb, int incb,
			  double *c, int ldc, int incc,
			  int ni, int nk, int nj, double alpha){
  double *tc1, *tc2;
  
  tc2 = c + ni * incc;
  incc -= nj * ldc;
  ni = nj * ldb;
  nj = nj * ldc;

  while( c != tc2 ){
    tc1 = c + nj;
    while( c != tc1 ){
      *c += alpha * yl_blas_inline_ddot (nk, b, incb, a, lda );
      c += ldc;
      b += ldb;
    }
    c += incc;
    b -= ni;
    a += inca;
  }
}

void yl_blas_dgemm_np_jik_opt(double *a, int lda, int inca,
			      double *b, int ldb, int incb,
			      double *c, int ldc, int incc,
			      int ni, int nk, int nj, char opt){
  
  if(opt & 2) tswap(lda, inca);
  if(opt & 1) tswap(ldb, incb); 
  yl_blas_dgemm_np_jik(a, lda, inca,
		       b, ldb, incb,
		       c, ldc, incc,
		       ni, nk, nj);
}


/*
void yl_blas_dgemm_auto(double *a, int lda, int inca,
			 double *b, int ldb, int incb,
			 double *c, int ldc, int incc,
			 int ni, int nk, int nj,
			 double alpha, double beta, YL_CBLAS_DGEMM_NP dgemmnp){
 
  int jj, kk;
  int bszc;

  
  for(jj=0; jj<nj; jj+=bszj){
    bszc = yfmin ( bszj , nj-jj );
    dmsf(beta, c + jj * ldc, ldc, incc, bszc, ni*incc);
    for(kk=0; kk<nk; kk+=bszk){
      dgemmnp(a + kk * lda, lda, inca,
	      b + kk * incb + jj * ldb,  ldb,  incb,
	      c + jj * ldc, ldc, incc,
	      ni, yfmin ( bszk , nk-kk ), bszc);
    }
  }

*/
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

  if(alpha==0.){
    dmsf(beta, c, ldc, incc, nj, ni*incc);
    return;
  }
  

  if(opt & 1)   tswap(ldb, incb);

  if(opt & 2) {
    bszk = BIGBZ;
    bsz  = SMLBZ;
    tswap(lda, inca);
    if(opt & 1)   bszk = SMLBZ;
  } else {
    bsz =  BIGBZ;
    bszk = SMLBZ;
  }

  //  cout<<(int)opt<<' '<<bsz<<' '<<bszk<<endl;

  int ii, kk;
  int bszc;
  for(ii=0; ii<ni; ii+=bsz){
    bszc = yfmin ( bsz , ni-ii );
    if(beta != 1.) dmsf(beta, c + ii * incc, ldc, incc, nj, bszc*incc);
    for(kk=0; kk<nk; kk+=bszk){
      if(alpha==1.) yl_blas_dgemm_np_jik(a + ii * inca + kk * lda, lda, inca, 
					 b + kk * incb, ldb, incb, 
					 c + ii * incc, ldc, incc, 
					 bszc, yfmin ( bszk , nk-kk ), nj);
      else         yl_blas_dgemm_npa_jik(a + ii * inca + kk * lda, lda, inca, 
					 b + kk * incb, ldb, incb, 
					 c + ii * incc, ldc, incc, 
					 bszc, yfmin ( bszk , nk-kk ), nj, alpha);
    }
  }
}


void yl_blas_dgemm_n_ijk(double *a, int lda, int inca,
			   double *b, int ldb, int incb,
			   double *c, int ldc, int incc,
			   int ni, int nk, int nj, char opt){
  if(opt == 0){
    for(int i=0; i<ni; i++){
      for(int j=0; j<nj; j++){
	double sum = 0;
	for(int k=0; k<nk; k++){
	  sum += a[k*lda+i*inca] * b[j*ldb+k*incb];
	}
	c[j*ldc+i*incc] += sum;
      }
    }
  } else if(opt==1){
    for(int i=0; i<ni; i++){
      for(int j=0; j<nj; j++){
	double sum = 0;
	for(int k=0; k<nk; k++){
	  sum += a[k*lda+i*inca] * b[k*ldb+j*incb];
	}
	c[j*ldc+i*incc] += sum;
      }
    }
  }
  else if(opt==3){
    for(int i=0; i<ni; i++){
      for(int j=0; j<nj; j++){
	double sum = 0;
	for(int k=0; k<nk; k++){
	  sum += a[i*lda+k*inca] * b[k*ldb+j*incb];
	}
	c[j*ldc+i*incc] += sum;
      }
    }
  } else {
    for(int i=0; i<ni; i++){
      for(int j=0; j<nj; j++){
	double sum = 0;
	for(int k=0; k<nk; k++){
	  sum += a[i*lda+k*inca] * b[j*ldb+k*incb];
	}
	c[j*ldc+i*incc] += sum;
      }
    }
  }
}
