#include <math.h>
#include <stdlib.h>
#include <mydll.h>

#define ROTATE(a,i,j,k,l) g=(*(a+i*lda+j*inca));h=(*(a+k*lda+l*inca));*(a+i*lda+j*inca)=g-s*(h+g*tau);*(a+k*lda+l*inca)=h+s*(g-h*tau);

void yl_cblas_nr_jacobi(int n, double *a, int lda, int inca, double *v, int ldv, int incv, double *d, int *nrot)
{
  int j,iq,ip,i;
  double tresh,theta,tau,t,sm,s,h,g,c,*b,*z;
  
  b = (double*) malloc( sizeof(double) * n);
  z = (double*) malloc( sizeof(double) * n);

  yl_blas_dmset(0, v, ldv, incv, n, n*incv);
  yl_blas_dset(n, 1, v, ldv+incv);

  for (ip=0;ip<n;ip++) {
    b[ip] = d[ip] = *(a + ip*(lda+inca));
    z[ip] = 0.0;
  }
  *nrot=0;
  for (i=1;i<=50;i++) {
    sm=0.0;
    for (ip=1;ip<n;ip++) {
      sm += yl_blas_dasum(n-ip, a+ip*(lda+inca)-lda, inca);
    }
    if (sm == 0.0) {
      free(z);
      free(b);
      return;
    }
    if (i < 4)
      tresh=0.2*sm/(n*n);
    else
      tresh=0.0;
    for (ip=0;ip<n-1;ip++) {
      for (iq=ip+1;iq<n;iq++) {
	double *ta = a+ip*lda+iq*inca;
	g=100.0*fabs(*ta);
	if (i > 4 && (double)(fabs(d[ip])+g) == (double)fabs(d[ip])
	    && (double)(fabs(d[iq])+g) == (double)fabs(d[iq]))
	  *ta=0.0;
	else if (fabs(*ta) > tresh) {
	  h=d[iq]-d[ip];
	  if ((double)(fabs(h)+g) == (double)fabs(h))
	    t=(*ta)/h;
	  else {
	    theta=0.5*h/(*ta);
	    t=1.0/(fabs(theta)+sqrt(1.0+theta*theta));
	    if (theta < 0.0) t = -t;
	  }
	  c=1.0/sqrt(1+t*t);
	  s=t*c;
	  tau=s/(1.0+c);
	  h=t*(*ta);
	  z[ip] -= h;
	  z[iq] += h;
	  d[ip] -= h;
	  d[iq] += h;
	  *ta=0.0;
	  for (j=0;j<ip;j++) {
	    ROTATE(a,j,ip,j,iq);
	  }
	  for (j=ip+1;j<iq;j++) {
	    ROTATE(a,ip,j,j,iq);
	  }
	  for (j=iq+1;j<n;j++) {
	    ROTATE(a,ip,j,iq,j);
	  }
	  for (j=0;j<n;j++) {
	    ROTATE(v,j,ip,j,iq);
	  }
	  ++(*nrot);
	}
      }
    }
    for (ip=0;ip<n;ip++) {
      b[ip] += z[ip];
      d[ip]=b[ip];
      z[ip]=0.0;
    }
  }
  free(z);
  free(b);
}
#undef ROTATE

