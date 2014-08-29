#include "mydll.h"
#include <cmath>
#include <cstdlib>
#include <cstdio>

double d_sign(double a, double b){
	return( b >= 0 ? fabs(a) : -fabs(a));
}


double bessi0(double x){
  double ax,ans;
  double y;
  if ((ax=fabs(x)) < 3.75) {
    y=x/3.75;
    y*=y;
    ans=1.0+y*(3.5156229+y*(3.0899424+y*(1.2067492
					 +y*(0.2659732+y*(0.360768e-1+y*0.45813e-2)))));
  } else {
    y=3.75/ax;
    ans=(exp(ax)/sqrt(ax))*(0.39894228+y*(0.1328592e-1
					  +y*(0.225319e-2+y*(-0.157565e-2+y*(0.916281e-2
									     +y*(-0.2057706e-1+y*(0.2635537e-1+y*(-0.1647633e-1
														  +y*0.392377e-2))))))));
  }
  return ans;
}


double YFfindmin(double *x,double *step, int *xpin,int varnum, CHIFunc chifunc, int iter, void *client){
  int j,i,k,adjust,limit,limit2;
  double abr,err;
  for(i=0; i<varnum; i++) printf("%d %g\n", xpin[i], x[i]);
  printf("\n%g\n", x[3]);
  err=(*chifunc)(client, x);
  for(j=0,limit=0;j<iter;j++){
    do{
      adjust=0;
      for(i=0;i<varnum;i++){
	if(!xpin[i]) continue;
	x[i]+=step[i];
	if((abr=(*chifunc)(client, x))<err) {err=abr;adjust++;}
	else{x[i]-=step[i];step[i]=-step[i];}
	limit2=0;
	do{
	  x[i]+=step[i];
	  if((abr=(*chifunc)(client, x))<err) {err=abr;adjust++;}
	  else {
	    x[i]-=step[i];
	    break;
	  }
	  if(limit2++>6000) {return(err);}
	}while(1);
      }
      if(limit++>100) {return(err);}
    }while(adjust>3);
    for(k=0;k<varnum;k++) step[k]/=3.;
  }
  return(err);
}

void free_vector(double *v,int nl,int nh){
  free((char*) (v+nl));
}


double nr_erf(double x)
{
  double t,z,ans;
  z=fabs(x);
  t=1.0/(1.0+0.5*z);
  ans=t*exp(-z*z-1.26551223+t*(1.00002368+t*(0.37409196+t*(0.09678418+
							   t*(-0.18628806+t*(0.27886807+t*(-1.13520398+t*(1.48851587+t*(-0.82215223+t*0.17087277)))))))));
  return  x >= 0.0 ? 1-ans : -1.0+ans;
}
#define MAXIT 100
#define EULER 0.5772156649
#define FPMIN 1.0e-30
#define EPS 1.0e-7
double expint(double x){
  int i;
  double a,b,c,d,del,fact,h,ans=0;
  if (x > 1.0) {
    b=x+1;
    c=1.0/FPMIN;
    d=1.0/b;
    h=d;
    for (i=1;i<=MAXIT;i++) {
      a = -i*(i);
      b += 2.0;
      d=1.0/(a*d+b);
      c=b+a/c;
      del=c*d;
      h *= del;
      if (fabs(del-1.0) < EPS) {
	    ans=h*exp(-x);
	    return ans;
      }
    }
    printf("continued fraction failed in expint");
  } else {
    ans = (-log(x)-EULER);
    fact=1.0;
    for (i=1;i<=MAXIT;i++) {
      fact *= -x/i;
      del = -fact/(i);
      ans += del;
      if (fabs(del) < fabs(ans)*EPS) return ans;
    }
    printf("series failed in expint");
  }
  return ans;
}
#undef MAXIT
#undef EULER
#undef FPMIN
#undef EPS
double bessj0(double x)
{
  double ax,z;
  double xx,y,ans,ans1,ans2;
  if ((ax=fabs(x)) < 8.0) {
    y=x*x;
    ans1=57568490574.0+y*(-13362590354.0+y*(651619640.7
					    +y*(-11214424.18+y*(77392.33017+y*(-184.9052456)))));
    ans2=57568490411.0+y*(1029532985.0+y*(9494680.718
					  +y*(59272.64853+y*(267.8532712+y*1.0))));
    ans=ans1/ans2;
  } else {
    z=8.0/ax;
    y=z*z;
    xx=ax-0.785398164;
    ans1=1.0+y*(-0.1098628627e-2+y*(0.2734510407e-4
				    +y*(-0.2073370639e-5+y*0.2093887211e-6)));
    ans2 = -0.1562499995e-1+y*(0.1430488765e-3
			       +y*(-0.6911147651e-5+y*(0.7621095161e-6
						       -y*0.934935152e-7)));
    ans=sqrt(0.636619772/ax)*(cos(xx)*ans1-z*sin(xx)*ans2);
  }
  return ans;
}




#define M1 259200
#define IA1 7141
#define IC1 54773
#define RM1 (1.0/M1)
#define M2 134456
#define IA2 8121
#define IC2 28411
#define RM2 (1.0/M2)
#define M3 243000
#define IA3 4561
#define IC3 51349

double ran1(int *idum)
{
  static long ix1,ix2,ix3;
  static double r[98];
  double temp;
  static int iff=0;
  int j;

  if (*idum < 0 || iff == 0) {
    iff=1;
    ix1=(IC1-(*idum)) % M1;
    ix1=(IA1*ix1+IC1) % M1;
    ix2=ix1 % M2;
    ix1=(IA1*ix1+IC1) % M1;
    ix3=ix1 % M3;
    for (j=1;j<=97;j++) {
      ix1=(IA1*ix1+IC1) % M1;
      ix2=(IA2*ix2+IC2) % M2;
      r[j]=(ix1+ix2*RM2)*RM1;
    }
    *idum=1;
  }
  ix1=(IA1*ix1+IC1) % M1;
  ix2=(IA2*ix2+IC2) % M2;
  ix3=(IA3*ix3+IC3) % M3;
  j=1 + ((97*ix3)/M3);
  if (j > 97 || j < 1) printf("RAN1: This cannot happen.");
  temp=r[j];
  r[j]=(ix1+ix2*RM2)*RM1;
  return temp;
}

#undef M1
#undef IA1
#undef IC1
#undef RM1
#undef M2
#undef IA2
#undef IC2
#undef RM2
#undef M3
#undef IA3
#undef IC3

double gammln(double xx)
{
  double x,tmp,ser;
  static double cof[6]={76.18009173,-86.50532033,24.01409822,
			-1.231739516,0.120858003e-2,-0.536382e-5};
  int j;

  x=xx-1.0;
  tmp=x+5.5;
  tmp -= (x+0.5)*log(tmp);
  ser=1.0;
  for (j=0;j<=5;j++) {
    x += 1.0;
    ser += cof[j]/x;
  }
  return -tmp+log(2.50662827465*ser);
}



#define PI 3.1415926535897932384626433832795

double poidev(double xm, int *idum)
{
  static double sq,alxm,g,oldm=(-1.0);
  double em,t,y;
  if (xm < 12.0) {
    if (xm != oldm) {
      oldm=xm;
      g=exp(-xm);
    }
    em = -1;
    t=1.0;
    do {
      em += 1.0;
      t *= ran1(idum);
    } while (t > g);
  } else {
    if (xm != oldm) {
      oldm=xm;
      sq=sqrt(2.0*xm);
      alxm=log(xm);
      g=xm*alxm-gammln(xm+1.0);
    }
    do {
      do {
	y=tan(PI*ran1(idum));
	em=sq*y+xm;
      } while (em < 0.0);
      em=floor(em);
      t=0.9*(1.0+y*y)*exp(em*alxm-gammln(em+1.0)-g);
    } while (ran1(idum) > t);
  }
  return em;
}

#undef PI
