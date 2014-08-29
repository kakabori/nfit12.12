#include "mydll.h"
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <iostream>

using namespace std;

void polyinterp(const double *xa, const double *ya,int n,double x,double *y,double *dy);


void Qromb::default_setup(){
  eps=1.0e-3;jmax=20; initNumOfIter=6; numForExterp=5;
  hForExterp = (double*) malloc(sizeof(double)*(jmax+1));
  sForExterp = (double*) malloc(sizeof(double)*(jmax+1));
  inith(1);
}

void Qromb::inith(double x){
  hForExterp[0]=x;
  for(int i=0; i<jmax; i++) hForExterp[i+1]=0.25*hForExterp[i];
}


double Qromb::trapzd(double a,double b,int n)
{
  register double x,s,del;
  double tnm;
  int j;
  
  if (n==0) {
    iteration=1;
    return (sum=0.5*(b-a)*(func(instance, a)+func(instance, b)));
  } else {
    tnm=iteration;
    del=(b-a)/tnm;
    x=a+0.5*del;
    for (s=0.0,j=1;j<=iteration;j++,x+=del) s += func(instance, x);
    iteration *= 2;
    return sum=0.5*(sum+(b-a)*s/tnm);
  }
}

double Qromb::calc(double a,double b)
{
  double ss, dss;
  for (int j=0;j<jmax;j++) {
    sForExterp[j]=trapzd(a,b,j);
    if (j > initNumOfIter) {
      polyinterp(&hForExterp[j-numForExterp],&sForExterp[j-numForExterp],numForExterp,0.0,&ss,&dss);
      if (fabs(dss) < eps*fabs(ss)) {
		return ss;
      }
    }
  }
  printf("Too many steps for integration");
  return ss;
}

double Qromb::trapzd(Integrand f, double a,double b,int n)
{
  register double x,s,del;
  double tnm;
  int j;
  
  if (n==0) {
    iteration=1;
    return (sum=0.5*(b-a)*(f(instance, a)+f(instance, b)));
  } else {
    tnm=iteration;
    del=(b-a)/tnm;
    x=a+0.5*del;
    for (s=0.0,j=1;j<=iteration;j++,x+=del) s += f(instance, x);
    iteration *= 2;
    return sum=0.5*(sum+(b-a)*s/tnm);
  }
}

double Qromb::calc(Integrand f, double a,double b)
{
  double ss, dss;
  for (int j=0;j<jmax;j++) {
    sForExterp[j]=trapzd(f,a,b,j);
    if (j > initNumOfIter) {
      polyinterp(&hForExterp[j-numForExterp],&sForExterp[j-numForExterp],numForExterp,0.0,&ss,&dss);
      if (fabs(dss) < eps*fabs(ss)) {
		return ss;
      }
    }
  }
  printf("Too many steps for integration");
  return ss;
}

void Qromb::custom_setup(double a){
	eps=a; initNumOfIter=8; numForExterp=5;
}


/*
  double myfunc(void *n, double a){
  return sin(a);
  }

  int main(){
  Qromb qromb;
  qromb.setIntegrand(myfunc);
  cout<<qromb.calc(0,3.1416/3)<<endl;
  double x=pow(4,12);
  for(int i=0;i<16; i++){
  x*=0.25;
  //   printf("%.16f,\n",qromb.hForExterp[i]);
  }
  
  return 0;
  }
*/
