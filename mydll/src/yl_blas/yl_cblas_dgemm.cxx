#include "mydll.h"
#include "yl_cblas_.h"
#include <iostream>

using namespace std;

/*C += A * B (ijk)  (_ | .) , interchange ldx vs. incx & a vs. b -> jik */
void yl_cblas_dgemm_np_ijk(double *a, int lda, int inca,
			    double *b, int ldb, int incb,
			    double *c, int ldc, int incc,
			    int ni, int nk, int nj){
  double *tc1, *tc2;
  
  tc2 = c + ni * ldc;
  ldc -= nj * incc;
  ni = nj * incb;
  nj = nj * incc;

  while( c != tc2 ){
    tc1 = c + nj;
    while( c != tc1 ){
      *c += yl_blas_inline_ddot (nk, a, inca, b, ldb );
      c += incc;
      b += incb;
    }
    c += ldc;
    a += lda;
    b -= ni;
  }
}

void yl_cblas_dgemm_npa_ijk(double *a, int lda, int inca,
			    double *b, int ldb, int incb,
			    double *c, int ldc, int incc,
			    int ni, int nk, int nj, double alpha){
  double *tc1, *tc2;
  
  tc2 = c + ni * ldc;
  ldc -= nj * incc;
  ni = nj * incb;
  nj = nj * incc;

  while( c != tc2 ){
    tc1 = c + nj;
    while( c != tc1 ){
      *c += alpha * yl_blas_inline_ddot (nk, a, inca, b, ldb );
      c += incc;
      b += incb;
    }
    c += ldc;
    a += lda;
    b -= ni;
  }
}
/*C += A * B (kij)  ( . - _ ) few misses, more write back*/
void yl_cblas_dgemm_np_kij(double *a, int lda, int inca,
			    double *b, int ldb, int incb,
			    double *c, int ldc, int incc,
			    int ni, int nk, int nj){
  double *ta1, *ta2;
  
  ta2 = a + ni * lda;
  lda -= nk * inca;
  ni = nk * ldb;
  nk = nk * inca;

  while( a != ta2 ){
    ta1 = a + nk;
    while( a != ta1 ){
      yl_blas_inline_daxpy_ (nj, *a, b, incb, c, incc );
      a += inca;
      b += ldb;
    }
    c += ldc;
    a += lda;
    b -= ni;
  }
}

// special
void yl_cblas_dgemm_np_ijk_s(double *a, int lda,
			      double *b, int ldb,
			      double *c, int ldc,
			      int ni, int nk, int nj){
  double *tc1, *tc2;
  
  tc2 = c + ni * ldc;
  ldc -= nj;
  while( c != tc2 ){
    tc1 = c + nj;
    while( c != tc1 ){
      *c += yl_blas_inline_ddot (nk, a, 1, b, ldb );
      c++;
      b++;
    }
    c += ldc;
    a += lda;
    b -= nj;
  }
}


#define  BIGBZ 1024
#define  SMLBZ 10
void yl_cblas_dgemm_opt(double *a, int lda, int inca,
			 double *b, int ldb, int incb,
			 double *c, int ldc, int incc,
			 int ni, int nk, int nj,
			 double alpha, double beta, char opt){
  int bszj, bszk;
  YL_CBLAS_DMSCAL dmsf = (beta==0.) ? yl_blas_dmset : yl_blas_dmscal;
  // YL_CBLAS_DGEMM_NP dgemmnp;

  if(alpha==0.){
    dmsf(0, c, ldc, incc, ni, nj*incc);
    return;
  }
  
  if(opt & 1)        {bszj = SMLBZ; bszk = yfmin(BIGBZ, ldb); tswap(ldb, incb);}
  else               {bszj = yfmin(BIGBZ, ldb); bszk = SMLBZ;}
  if(opt & 2)        {bszk = SMLBZ; tswap(lda, inca);}

  //cout<<(int)opt<<' '<<bszj<<' '<<bszk<<endl;
  int jj, kk;
  int bszc;
  for(jj=0; jj<nj; jj+=bszj){
    bszc = yfmin ( bszj , nj-jj );
    dmsf(beta, c + jj * incc, ldc, incc, ni, bszc*incc);
    for(kk=0; kk<nk; kk+=bszk){
      yl_cblas_dgemm_np_ijk(a + kk * inca, lda, inca,
			    b + jj * incb + kk * ldb,  ldb,  incb,
			    c + jj * incc, ldc, incc,
			    ni, yfmin ( bszk , nk-kk ), bszc);
    }
  }
}

void yl_cblas_dgemm(double *a, int lda, int inca,
		    double *b, int ldb, int incb,
		    double *c, int ldc, int incc,
		    int ni, int nk, int nj,
  		    double alpha, double beta){
  int bszj, bszk;
  YL_CBLAS_DMSCAL dmsf = (beta==0.) ? yl_blas_dmset : yl_blas_dmscal;
  if(alpha==0.){
    if( ldc > incc )
      dmsf(0, c, ldc, incc, ni, nj*incc);
    else
      dmsf(0, c, incc, ldc, nj, ni*ldc);      
    return;
  }

  if( ldb < incb ) {
    bszj = SMLBZ; bszk = yfmin(BIGBZ, incb);
  } else {
    bszj = yfmin(BIGBZ, ldb); bszk = SMLBZ;
  }
  if( lda < inca ) bszk = SMLBZ;


  int jj, kk;
  int bszc;

  for(jj=0; jj<nj; jj+=bszj){
    bszc = yfmin ( bszj , nj-jj );
    dmsf(beta, c + jj * incc, ldc, incc, ni, bszc*incc);
    for(kk=0; kk<nk; kk+=bszk){
      if(alpha==1.0)      yl_cblas_dgemm_np_ijk(a + kk * inca, lda, inca,
						b + jj * incb + kk * ldb,  ldb,  incb,
						c + jj * incc, ldc, incc,
						ni, yfmin ( bszk , nk-kk ), bszc);
      else                yl_cblas_dgemm_npa_ijk(a + kk * inca, lda, inca,
						b + jj * incb + kk * ldb,  ldb,  incb,
						c + jj * incc, ldc, incc,
						ni, yfmin ( bszk , nk-kk ), bszc, alpha);

    }
  }
}


void yl_cblas_dgemm_n_ijk(double *a, int lda, int inca,
			   double *b, int ldb, int incb,
			   double *c, int ldc, int incc,
			   int ni, int nk, int nj, char opt){
  if(opt == 0){
    for(int i=0; i<ni; i++){
      for(int j=0; j<nj; j++){
	double sum = 0;
	for(int k=0; k<nk; k++){
	  sum += a[i*lda+k*inca] * b[k*ldb+j*incb];
	}
	c[i*ldc+j*incc] += sum;
      }
    }
  } else if(opt==1){
    for(int i=0; i<ni; i++){
      for(int j=0; j<nj; j++){
	double sum = 0;
	for(int k=0; k<nk; k++){
	  sum += a[i*lda+k*inca] * b[j*ldb+k*incb];
	}
	c[i*ldc+j*incc] += sum;
      }
    }
  }
  else if(opt==3){
    for(int i=0; i<ni; i++){
      for(int j=0; j<nj; j++){
	double sum = 0;
	for(int k=0; k<nk; k++){
	  sum += a[k*lda+i*inca] * b[j*ldb+k*incb];
	}
	c[i*ldc+j*incc] += sum;
      }
    }
  } else {
    for(int i=0; i<ni; i++){
      for(int j=0; j<nj; j++){
	double sum = 0;
	for(int k=0; k<nk; k++){
	  sum += a[k*lda+i*inca] * b[k*ldb+j*incb];
	}
	c[i*ldc+j*incc] += sum;
      }
    }
  }
}




void yl_cblas_dgemm2(double *a, int lda, int inca,
		     double *b, int ldb, int incb,
		     double *c, int ldc, int incc,
		     int ni, int nk, int nj){
  int i,j,k,jj, kk;
  
  const int bszj=5, bszk=5;
  if(1 && inca==1 && incb==1 && incc==1){
    for(jj=0; jj<nj; jj+=bszj){
      for(kk=0; kk<nk; kk+=bszk){
	for(i=0; i<ni; i++){
	  for(j=jj; j<yfmin(jj+bszj, nj); j++){
	    double sum=0;
	    for(k=kk; k<yfmin(kk+bszk,nk); k++){
	      sum += a[i*lda+k] * b[k*ldb+j];
	    }
	    c[i*ldc+j] += sum;
	  }
	}
      }
    }
  } else {
    for(jj=0; jj<nj; jj+=bszj){
      /*
	for(i=0; i<ni; i++)
	for(j=jj; j<yfmin(jj+bszj, nj); j++)
	c[i*ldc+j*incc] = 0;
      */
      for(kk=0; kk<nk; kk+=bszk){
	
	for(i=0; i<ni; i++){
	  for(j=jj; j<yfmin(jj+bszj, nj); j++){
	    double sum=0;
	    for(k=kk; k<yfmin(kk+bszk,nk); k++){
	      sum += a[i*lda+k*inca] * b[k*ldb+j*incb];
	    }
	    c[i*ldc+j*incc] += sum;
	  }
	}
      }
    }
  }
}

