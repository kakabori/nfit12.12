#include <math.h>
#include <mydll.h>


void lnsrch(int n,
	    double* xold, double fold, double* g, double* p,
	    double *x, double *f,
	    double stpmax, int *check,
	    YLMinimizationFunc func, void* client, DfpminCtrl *ctrl){
  int i;
  double a,alam,alam2,alamin,b,disc,f2,rhs1,rhs2,slope,sum,temp,test,tmplam;

  *check=0;
  sum = yl_blas_dnrm2(n, p, 1);
  if (sum > stpmax) yl_blas_dscal(n, stpmax/sum, p, 1);
  slope = yl_blas_ddot(n, g, 1, p, 1);
  if (slope >= 0.0) printf("Roundoff problem in lnsrch.");
  test=0.0;
  for (i=0;i<n;i++) {
    temp = fabs(p[i]) / yl_max( fabs(xold[i]), 1.0 );
    if (temp > test) test=temp;
  }
  alamin=ctrl->tolx/test;
  alam2=alam=1.0;
  f2 = *f;
  for (;;) {
    for (i=0;i<n;i++) x[i]=xold[i]+alam*p[i];
    *f=(*func)(x,client,NULL);
    if (alam < alamin) {
      yl_blas_dcopy(n, xold, 1, x, 1); 
      *check=1;
      return;
    } else if (*f <= fold+ctrl->alf*alam*slope) return;
    else {
      if (alam == 1.0)
	tmplam = -slope/(2.0*(*f-fold-slope));
      else {
	rhs1 = *f-fold-alam*slope;
	rhs2=f2-fold-alam2*slope;
	a=(rhs1/(alam*alam)-rhs2/(alam2*alam2))/(alam-alam2);
	b=(-alam2*rhs1/(alam*alam)+alam*rhs2/(alam2*alam2))/(alam-alam2);
	if (a == 0.0) tmplam = -slope/(2.0*b);
	else {
	  disc=b*b-3.0*a*slope;
	  if (disc < 0.0) tmplam=0.5*alam;
	  else if (b <= 0.0) tmplam=(-b+sqrt(disc))/(3.0*a);
	  else tmplam=-slope/(b+sqrt(disc));
	}
	if (tmplam > 0.5*alam)
	  tmplam=0.5*alam;
      }
    }
    alam2=alam;
    f2 = *f;
    alam=yl_max(tmplam,0.1*alam);
  }
}


void dfpmin(double* p, int n, int *iter, double *fret,
	    YLMinimizationFunc func,
	    YLMinimizationDFunc dfunc, void* client, DfpminCtrl *ctrl){
  int check,i,its,j;
  int itmax = *iter;
  double den,fac,fad,fae,fp,stpmax,sum=0.0,sumdg,sumxi,temp,test;
  double *dg,*g,*hdg,**hessin,*pnew,*xi;
  
  tvec(&g,n*5);
  tmat(&hessin, n, n);
  dg = g + n;
  hdg = dg + n;
  pnew = hdg + n;
  xi = pnew + n;
  
  fp=(*func)(p, client, NULL);
  (*dfunc)(p,g,client,NULL);
  
  for (i=0;i<n;i++) {
    hessin[i][i]=1.0;
    xi[i] = -g[i];
    sum += p[i]*p[i];
  }
  stpmax = ctrl->stpmx * yl_max(sqrt(sum),(double)n);
  for (its=1;its<=itmax;its++) {
    *iter=its;
    lnsrch(n,p,fp,g,xi,pnew,fret,stpmax,&check,func, client, ctrl);
    fp = *fret;
    for (i=0;i<n;i++) {
      xi[i]=pnew[i]-p[i];
      p[i]=pnew[i];
    }
    test=0.0;
    for (i=0;i<n;i++) {
      temp=fabs(xi[i]) / yl_max(fabs(p[i]),1.0);
      if (temp > test) test=temp;
    }
    if (test < ctrl->tolx)  goto freeall;
    yl_blas_dcopy(n, g, 1, dg, 1);
    (*dfunc)(p,g,client,NULL);
    test=0.0;
    den=yl_max(*fret,1.0);
    for (i=0;i<n;i++) {
      temp=fabs(g[i])*yl_max(fabs(p[i]),1.0)/den;
      if (temp > test) test=temp;
    }
    if (test < ctrl->gtol) goto freeall;
    for (i=0;i<n;i++) dg[i]=g[i]-dg[i];
    
    
    for (i=0;i<n;i++) {
      hdg[i] = yl_blas_ddot(n, hessin[i], 1, dg, 1);
    }
    fac=fae=sumdg=sumxi=0.0;
    for (i=0;i<n;i++) {
      fac += dg[i]*xi[i];
      fae += dg[i]*hdg[i];
      sumdg += sqr(dg[i]);
      sumxi += sqr(xi[i]);
    }
    if (fac > sqrt(ctrl->eps*sumdg*sumxi)) {
      fac=1.0/fac;
      fad=1.0/fae;
      for (i=0;i<n;i++) dg[i]=fac*xi[i]-fad*hdg[i];
      for (i=0;i<n;i++) {
	for (j=i;j<n;j++) {
	  hessin[i][j] += fac*xi[i]*xi[j]
	    -fad*hdg[i]*hdg[j]+fae*dg[i]*dg[j];
	  hessin[j][i]=hessin[i][j];
	}
      }
    }
    for (i=0;i<n;i++) {
      xi[i] = -yl_blas_ddot(n, hessin[i], 1, g, 1);
    }
  }
  printf("too many iterations in dfpmin");
 freeall:
  free(g);
  freetmat(hessin);
}
