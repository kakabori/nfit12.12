#include "mydll.h"

void NrDct::initialize(int tn){
  //  printf("sss0\n"); fflush(stdout);
  if(tn == n) return;
  free(sinarr); free(cosarr);
  sinarr = (double *) malloc(sizeof(double)*tn);
  cosarr = (double *) malloc(sizeof(double)*tn);
  //printf("sss1\n"); fflush(stdout);
  for(int i=0; i<tn; i++){
    sinarr[i] = sin(i*YFPI/tn);
    cosarr[i] = cos(i*YFPI/tn);
  }
  n = tn;
  p = rfftw_create_plan(tn, FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE);
  //printf("sss2\n"); fflush(stdout);
  free(work);
  //printf("sss3\n"); fflush(stdout);
  work = (double *) malloc(sizeof(double)*tn);
  //printf("sss4\n"); fflush(stdout);
}
/*
void NrDct::dct1(double*y){
  double wi=0, wr=1, wpi, wpr, wtemp, y1, y2;
  double sum = 0.5 * (y[0] - y[n]);

  wtemp = sin(0.5*YFPI/n);
  wpr = -2*wtemp*wtemp;
  wpi = sin(YFPI/n);

  y[0] = 0.5 * (y[0] + y[n]);
  for(int j=1; j<n/2; j++){
    wr = (wtemp=wr)*wpr-wi*wpi+wr;
    wi=wi*wpr+wtemp*wpi+wi;

    y1 = 0.5 * (y[j] + y[n-j]);
    y2 = y[j] - y[n-j];
    y[j] = y1 - wi * y2;
    y[n-j] = y1 + wi * y2;
    sum += wr * y2;
  }
  rfftw_one(p, y, work);
  y[1] = sum;
  y[0] = work[0];
  y[n] = work[n/2];
  for(int j=1; j<n/2; j++){
    y[j*2] = work[j];
    y[j*2+1] = y[j*2-1] - work[n-j];
  }
}
*/
void NrDct::dct1(double*y){
  int tn=n;
  int mid = (tn+1)/2;
  double y1, y2;
  double sum = 0.5 * (y[0] - y[tn]);

  y[0] = 0.5 * (y[0] + y[tn]);
  for(int j=1; j<mid; j++){
    y1 = 0.5 * (y[j] + y[tn-j]);
    y2 = y[j] - y[tn-j];
    y[j] = y1 - sinarr[j] * y2;
    y[tn-j] = y1 + sinarr[j] * y2;
    sum += cosarr[j] * y2;
  }
  rfftw_one(p, y, work);
  y[1] = sum;
  y[0] = work[0];
  y[tn] = work[tn/2];
  for(int j=1; j<mid; j++){
    y[j*2] = work[j];
    y[j*2+1] = y[j*2-1] - work[tn-j];
  }
}
