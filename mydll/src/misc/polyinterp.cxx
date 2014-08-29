#include "mydll.h"
#include <cmath>
#include <cstdlib>
#include <cstdio>

/* assume x seperated by 1, z is x-x0*/
double polint3(double z,double *ya){
  return z*((z-1)/2*(ya[2]-ya[0])-(z-2)*(ya[1]-ya[0]))+ya[0];
}

double polint4t(double xmx0,double *ya){
  return (ya[0]*(3-xmx0)+ya[3]*xmx0)*(xmx0-1)*(xmx0-2)/6+
    (ya[2]*(1-xmx0)+ya[1]*(xmx0-2))*(xmx0-3)*xmx0/2;
}
double polint3_log(double xmxa,float *ya){
  return xmxa*50*((100*xmxa-2)*(ya[0]+ya[2]-2*ya[1])+ya[2]-ya[0])+ya[0];
}
double polint30_log(double x,float *ya){
  return (0.9772372212*ya[0]*(x-1)-42.93136751*ya[1]*x)*(x-1.023292992)+41.95413029*ya[2]*(x-1)*x;
}


//return the offset
int mydllNearNeighbour(double *xa, int incx, int n, double q){
	double* np = xa;
	double* ori = xa;
	double* end = xa+n*incx;
	double dif = fabs( q - *xa );
	double dift;

	while( (xa+=incx) < end){
		if ( (dift = fabs(q - *xa)) < dif ){
			np = xa;
			dif = dift;
		}
	}
	return (np-ori);
}


//
//Neville’s algorithm
//
void polyinterp(const double *xa, const double *ya,int n,double x,double *y,double *dy)
{
  int i,m,ns=0;
  double den,dif,dift,ho,hp;
  double *c=(double *)malloc(n*sizeof(double)),
	     *d=(double *)malloc(n*sizeof(double));
  dif=fabs(x-xa[0]);
  for (i=0;i<n;i++) {
    if ( (dift=fabs(x-xa[i])) < dif) {
      ns=i;
      dif=dift;
    }
    d[i]=c[i]=ya[i];
  }
  *y=ya[ns];
  for (m=1;m<n;m++) {
    for (i=0;i<n-m;i++) {
      ho=xa[i]-x;
      hp=xa[i+m]-x;
      //w=c[i+1]-d[i];
      //if ( (den=ho-hp) == 0.0){	printf("Error in routine POLINT");}
	  // den = w / den;
      den = (c[i+1]-d[i]) / (ho-hp);
      d[i]=hp*den;
      c[i]=ho*den;
    }
    *y += (*dy=(2*ns<(n-m) ? c[ns] : d[--ns]));
  }
  free(d);
  free(c);
}

