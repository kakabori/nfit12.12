#include "mydll.h"
#include <math.h>
#include <stdio.h>
#include <iostream>

using namespace std;


void yl_blas_mprint(ostream& out, int m, int n, double *a, int lda, int inca, int trans){
  if(trans) {
    tswap(lda, inca);
  }
  out<<'{';
  for(int i=0; i<m; i++){
    out<<'{';
    for(int j=0; j<n; j++){
      out<<a[j*lda + i*inca]<<' ';
    }
    if(i==m-1) out<<"}}"<<endl;
    else out<<'}'<<endl;
  }
}

void yl_blas_vprint(ostream& out, int n, double *a, int inca){
  out<<"{";
  for(int i=0; i<n; i++){
    out<<a[i*inca]<<' ';
  }
  out<<"}"<<endl;
}
