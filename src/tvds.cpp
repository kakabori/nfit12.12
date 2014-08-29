/*
This does transformation  qx<-->px
See tvds.h for qz<--pz
*/
// use DSscale for similar triangle
// use DSfit2 to fit 2 sample given S


#include <cmath>
#include <iostream>
#include <vector>

#include "tvds.h"

using namespace std;

/*  Intermediate steps */
double TvDs2::getl(double theta, double dspacing){//take the order num; return in pixel unit
  /*given \theta and D, return l (l is defined in the appendix) */
  theta=asin(theta*wavelength/2/dspacing); // theta_lipid
  theta=acos(cos(theta)*nindex); //calc theta_air
  return ((R*cos(theta)+tan(2*theta)*(s+R*sin(theta)))/(cos(alpha)-tan(2*theta)*sin(alpha))+bc2b)/pz;
}

double TvDs2::getSinTheta(double pixel){//length in mm (from beaker center to the pixel position
  /* given pixels, return sin(\theta) */
  double sinTheta=0, lstar, length=pixel*pz-bc2b;
  for(int i=0;i<20;i++){ // n=20 will be enough to get precise value of x0
    lstar=length/sqrt(1-sinTheta*sinTheta)*(cos(alpha)-tan(2*asin(sinTheta))*sin(alpha));
    sinTheta=(sqrt(s*s+2*(lstar-R)*lstar)-s)/2/lstar;
  }
  sinTheta=sin(acos( cos(asin(sinTheta))/nindex )); // calc bragg theta from apparent theta
  return sinTheta;
}

/* The two transformations */
double TvDs2::getqr(double rpix){
  /* given pixels, return qr */
  return 2*YFPI/wavelength*sin(atan((rpix-qrzero)*pz/s));
}

double TvDs2::getPixQr(double qr){
  /* given qr, return pixels */
  return s/pz*tan(asin(qr*wavelength/2/YFPI));
}
