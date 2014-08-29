#include "mydll.h"
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include "cminpack.h"   //included for cminpack covar function
#define real __cminpack_real__  //included for cminpack covar function

double dpmpar(int i){
  if(i==1) return(2.22044604926e-16);
  else if(i==2) return(2.22507385852e-308);
  else return(1.79769313485e+308);
}

/*
the original enorm routine is more complex to deal with overflow and underflow
*/
double enorm(int n,double *x){
  double s2=0;
  for(int i=0;i<n;i++) s2+=x[i]*x[i];
  return(sqrt(s2));
}

int fdjac2(LmdifFunc fcn,int m,int n,double *x,double *fvec,double **fjac,
	   int ldfjac,double *epsfcn,double *wa,void *client){
  int i,j;
  double h,temp;
  for(j=0;j<n;j++){
    if( (h=epsfcn[j]*fabs(temp=x[j])) == 0.0 ) h=epsfcn[j];
    x[j]+=h;
    if((*fcn)(m,n,x,wa,client,(void *)2)) return 1;
    x[j]=temp;
    for(i=0;i<m;i++) fjac[j][i] = (wa[i]-fvec[i]) / h;
  }
  return 0;
}


void qrfac(int m,int n,double **a,int lda,int *ipvt,int lipvt,double *rdiag,
	   double *acnorm,double *wa){
  int i,j,k,kmax,minmn;
  double  ajnorm,epsmch,sum,temp;
  epsmch=dpmpar(2);
  for(j=0;j<n;j++){
    acnorm[j]=enorm(m,a[j]);
    rdiag[j]=acnorm[j];
    wa[j]=rdiag[j];
    ipvt[j]=j;
  }
  minmn=m>n?n:m;
  for(j=0;j<minmn;j++){
    kmax=j;
    for(k=j;k<n;k++) if(rdiag[k]>rdiag[kmax]) kmax=k;
    if(kmax!=j){
      for(i=0;i<m;i++){
	temp=a[j][i];a[j][i]=a[kmax][i];a[kmax][i]=temp;
      }
      rdiag[kmax]=rdiag[j];wa[kmax]=wa[j];
      k=ipvt[j];
      ipvt[j]=ipvt[kmax];
      ipvt[kmax]=k;
    }

    ajnorm=enorm(m-j,&a[j][j]);
    if(ajnorm==0) {rdiag[j]=-ajnorm;continue;}
    if(a[j][j]<0) ajnorm=-ajnorm;
    for(i=j;i<m;i++) a[j][i]/=ajnorm;
    a[j][j]+=1.0;
    for(k=j+1;k<n;k++){
      for(i=j,sum=0;i<m;i++) sum+=a[j][i]*a[k][i];
      temp=sum/a[j][j];
      for(i=j;i<m;i++) a[k][i]-=temp*a[j][i];
      if(rdiag[k]==0) continue;
      temp=a[k][j]/rdiag[k];
      rdiag[k]=rdiag[k]*sqrt(yfmax(0.0,1.0-temp*temp));
      if(0.05*(rdiag[k]/wa[k])*(rdiag[k]/wa[k])>epsmch) continue;
      rdiag[k]=enorm(m-j-1,&a[k][j+1]);wa[k]=rdiag[k];
    }
    rdiag[j]=-ajnorm;
  }
  return;
}


void qrsolv(int n,double **r,int ldr,int *ipvt,double *diag,double *qtb,
	    double *x,double *sdiag,double *wa){
  int i,j,k,nsing;
  double  cos,cotan,qtbpj,sin,sum,tan,temp;
  for(j=0;j<n;j++){
    for(i=j;i<n;i++) r[j][i]=r[i][j];
    x[j]=r[j][j];wa[j]=qtb[j];
  }
  for(j=0;j<n;j++){
    if(diag[ipvt[j]]!=0){
      for(k=j;k<n;k++) sdiag[k]=0;
      sdiag[j]=diag[ipvt[j]];qtbpj=0;
      for(k=j;k<n;k++){
	if(sdiag[k]==0) continue;
	if(fabs(r[k][k])<fabs(sdiag[k])){
	  cotan=r[k][k]/sdiag[k];
	  sin=0.5/sqrt(0.25+0.25*cotan*cotan);cos=sin*cotan;
	}
	else{
	  tan=sdiag[k]/r[k][k];
	  cos=0.5/sqrt(0.25+0.25*tan*tan);sin=cos*tan;
	}
	r[k][k]=cos*r[k][k]+sin*sdiag[k];
	temp=cos*wa[k]+sin*qtbpj;
	qtbpj=-sin*wa[k]+cos*qtbpj;wa[k]=temp;
	for(i=k+1;i<n;i++){
	  temp=cos*r[k][i]+sin*sdiag[i];
	  sdiag[i]=-sin*r[k][i]+cos*sdiag[i];r[k][i]=temp;
	}
      }
    }
    sdiag[j]=r[j][j];r[j][j]=x[j];
  }
  nsing=n;
  for(j=0;j<n;j++){
    if(sdiag[j]==0&&nsing==n) nsing=j;
    if(nsing<n) wa[j]=0;
  }
  for(k=1;k<=nsing;k++){
    for(j=nsing-k,sum=0,i=j+1;i<nsing;i++) sum+=r[j][i]*wa[i];
    wa[j]=(wa[j]-sum)/sdiag[j];
  }
  for(j=0;j<n;j++) x[ipvt[j]]=wa[j];
  return;
}

void lmpar(int n,double **r,int ldr,int *ipvt,double *diag,
	   double *qtb,double delta,double par,double *x,
	   double *sdiag,double *wa1,double *wa2){
  int i,iter=0,j,k,nsing;
  double  dxnorm,dwarf,fp,gnorm,parc,parl,paru,sum,temp;
  dwarf=dpmpar(2);
  nsing=n;
  for(j=0;j<n;j++){
    if(r[j][j]==0&&nsing==n) nsing=j;
    if(nsing<n) {wa1[j]=0;continue;}
    wa1[j]=qtb[j];
  }
  for(k=1;k<=nsing;k++){
    j=nsing-k;wa1[j]/=r[j][j];
    for(i=0,temp=wa1[j];i<j;i++) wa1[i]-=r[j][i]*temp;
  }
  for(j=0;j<n;j++) x[ipvt[j]]=wa1[j];
  for(j=0;j<n;j++) wa2[j]=diag[j]*x[j];
  dxnorm=enorm(n,wa2);
  fp=dxnorm-delta;
  if(fp<=0.1*delta){par=0;return;}
  parl=0;
  if(nsing>=n){
    for(j=0;j<n;j++) wa1[j]=diag[ipvt[j]]*wa2[ipvt[j]]/dxnorm;
    for(j=0;j<n;j++){
      for(i=0,sum=0;i<j;i++) sum+=r[j][i]*wa1[i];
      wa1[j]=(wa1[j]-sum)/r[j][j];
    }
    temp=enorm(n,wa1);parl=fp/delta/temp/temp;
  }
  for(j=0;j<n;j++){
    for(i=0,sum=0;i<=j;i++) sum+=r[j][i]*qtb[i];
    wa1[j]=sum/diag[ipvt[j]];
  }
  gnorm=enorm(n,wa1);
  paru=gnorm/delta;
  if(paru==0) paru=dwarf/yfmin(delta,0.1);
  par=yfmax(par,parl);par=yfmin(par,paru);
  if(par==0) par=gnorm/dxnorm;
  do{
    iter++;
    if(par==0) par=yfmax(dwarf,0.001*paru);
    temp=sqrt(par);
    for(j=0;j<n;j++) wa1[j]=temp*diag[j];
    qrsolv(n,r,ldr,ipvt,wa1,qtb,x,sdiag,wa2);
    for(j=0;j<n;j++) wa2[j]=diag[j]*x[j];
    dxnorm=enorm(n,wa2);
    temp=fp;
    fp=dxnorm-delta;
    if(fabs(fp) <= 0.1*delta || (parl == 0 && fp <= temp && temp < 0) || iter==10) return;
    for(j=0;j<n;j++){wa1[j]=diag[ipvt[j]]*(wa2[ipvt[j]]/dxnorm);}
    for(j=0;j<n;j++){
      wa1[j]/=sdiag[j];temp=wa1[j];
      for(i=j+1;i<n;i++) wa1[i]-=r[j][i]*temp;
    }
    temp=enorm(n,wa1);parc=((fp/delta)/temp)/temp;
    if(fp>0) parl=yfmax(parl,par);
    if(fp<0) paru=yfmin(paru,par);
    par=yfmax(parl,par+parc);
  }while(1);
}

Lmdif::Lmdif() 
{ 
  x = fvec = qtf = diag = NULL;
  wa1 = wa2 = wa3 = wa4 = NULL; 
  fjac = NULL;
  fjac_cminpack = NULL; 
  ipvt = NULL; n = 0; 
}

void Lmdif::cleanup(){
  int i;
  for(i=0;i<n;i++){
    free(fjac[i]);
    fjac[i]=NULL;
  }
  n=0;
  free(fjac);fjac=NULL;
  free(wa1);wa1=NULL;
  free(wa2);wa2=NULL;
  free(wa3);wa3=NULL;
  free(wa4);wa4=NULL;
  free(x);x=NULL;
  free(fvec);fvec=NULL;
  free(qtf);qtf=NULL;
  free(ipvt);ipvt=NULL;
  free(diag);diag=NULL;
  free(fjac_cminpack); fjac_cminpack = NULL;
}

/*
pfcn  : the function to maximize or minimize
m     :
n     :
maxfev: maximum function evaluation
ftol  : tolerance of
gtol  : tolerance of
xtol  : tolerance of
factor: a parameter in Lavenberg algrithom to control the initial jump step
epsfcn: for dirivative calculation.
client:
*/
void Lmdif::setup(LmdifFunc f, int _m, int _n, int _maxfev, double _ftol,
		  double _gtol, double _xtol, double _factor, double *_epsfcn, void *_client)
{
  int i;
  cleanup();

  pfcn = f;
  m = _m;
  n = _n;
  maxfev = _maxfev;
  ftol = _ftol;
  gtol = _gtol;
  xtol = _xtol;
  factor = _factor;
  epsfcn = _epsfcn;
  client = _client;

  ldfjac = m;


  fjac=(double **)malloc(sizeof(double *)*n);
  for(i=0;i<n;i++) fjac[i]=(double *)malloc(sizeof(double)*m);
  x=(double *)malloc(sizeof(double)*n);
  fvec=(double *)malloc(sizeof(double)*m);
  diag=(double *)malloc(sizeof(double)*n);
  qtf=(double *)malloc(sizeof(double)*n);
  ipvt=(int *)malloc(sizeof(int)*n);
  wa1=(double *)malloc(sizeof(double)*n);
  wa2=(double *)malloc(sizeof(double)*n);
  wa3=(double *)malloc(sizeof(double)*n);
  wa4=(double *)malloc(sizeof(double)*m);
  fjac_cminpack = (double *)malloc(sizeof(double) * (m*n)); //for cminpack covar function, added by KA, 9/20/2012
}

void Lmdif::init(double *b, int n){
  for(int i = 0; i < n; i++) x[i] = b[i];
}

void Lmdif::init(double *b[], int n){
  for(int i=0; i<n; i++) x[i] = *(b[i]);
}

#define NFEV_PLUS_ONE(a,b) {if((++a)==b) return 5;}
/*
this is the function we call, but fully understanding of the above two routines is sufficient for our work.
*/

int Lmdif::fit(){
  int i,iter=0,j,info=0,nfev=0;
  double actred,delta=0,dirder,epsmch,fnorm,fnorm1,gnorm,par=0,pnorm,prered,ratio,sum,temp,temp1,xnorm=0;
  epsmch=dpmpar(2);
  if(n<=0||m<n||ldfjac<m||ftol<0||xtol<0||gtol<0||maxfev<=0||factor<=0) {
    printf("wrong argument\n");
    return(0);
  }
  printf("maxfev=%d\n",maxfev);
  if((*pfcn)(m,n,x,fvec,client,(void *)0)) return 1; NFEV_PLUS_ONE(nfev,maxfev);
  fnorm=enorm(m,fvec);iter=1;
  do{
    if(fdjac2(pfcn,m,n,x,fvec,fjac,ldfjac,epsfcn,wa4,client)) return 1;
    qrfac(m,n,fjac,ldfjac,ipvt,n,wa1,wa2,wa3);
    if(iter==1){
      for(j=0;j<n;j++){diag[j]=wa2[j]; if(wa2[j]==0) diag[j]=1.0;}
      for(j=0;j<n;j++) wa3[j]=diag[j]*x[j];
      xnorm=enorm(n,wa3);delta=factor*xnorm;
      if(delta==0) delta=factor;
    }
    for(i=0;i<m;i++) wa4[i]=fvec[i];
    for(j=0;j<n;j++){
      if(fjac[j][j]!=0){
	for(i=j,sum=0;i<m;i++) sum+=fjac[j][i]*wa4[i];
	for(i=j;i<m;i++) wa4[i]-=fjac[j][i]*sum/fjac[j][j];
      }
      fjac[j][j]=wa1[j];qtf[j]=wa4[j];
    }
    gnorm=0;
    if(fnorm!=0)
      for(j=0;j<n;j++){
	if(wa2[ipvt[j]]==0) continue;
	for(i=0,sum=0;i<=j;i++) sum+=fjac[j][i]*(qtf[i]/fnorm);
	gnorm=yfmax(gnorm,fabs(sum/wa2[ipvt[j]]));
      }
    if(gnorm<=gtol) return(4);
    for(j=0;j<n;j++) diag[j]=yfmax(diag[j],wa2[j]);
    /** beginning of the inner loop. **/
    do{
    	lmpar(n,fjac,ldfjac,ipvt,diag,qtf,delta,par,wa1,wa2,wa3,wa4);
    	for(j=0;j<n;j++){
    		wa1[j]=-wa1[j];
    		wa2[j]=x[j]+wa1[j];
    		wa3[j]=diag[j]*wa1[j];
    	}
    	pnorm=enorm(n,wa3);
      	if(iter==1) delta=yfmin(delta,pnorm);
      	if((*pfcn)(m,n,wa2,wa4,client,(void *)1)) return 1;
      	NFEV_PLUS_ONE(nfev,maxfev);
      	fnorm1=enorm(m,wa4);
      	actred=-1.0;//
      	if(0.1*fnorm1<fnorm) actred=1.0-(fnorm1/fnorm)*(fnorm1/fnorm);
      	for(j=0;j<n;j++){
			wa3[j]=0;temp=wa1[ipvt[j]];
			for(i=0;i<=j;i++)
				wa3[i]+=fjac[j][i]*temp;
      	}
      	temp=enorm(n,wa3)/fnorm;
      	temp*=temp;
      	temp1=(sqrt(par)*pnorm)/fnorm;temp1*=temp1;
      	prered=temp+temp1/0.5;
      	dirder=-(temp+temp1);
      	ratio=0;
      	if(prered!=0) ratio=actred/prered;
      	if(ratio<=0.25){
			temp=actred>=0?0.5:0.5*dirder/(dirder+0.5*actred);
			if(0.1*fnorm1>=fnorm||temp<0.1) temp=0.1;
			delta=temp*yfmin(delta,pnorm/0.1);
			par/=temp;
      	}
      	else if(par!=0||ratio<0.75) ;
      	else {
      		delta=pnorm/0.5;
      		par*=0.5;
      	}
      	if(ratio>=0.0001){
			for(j=0;j<n;j++){
				x[j]=wa2[j];
				wa2[j]=diag[j]*x[j];
				printf("x[%d] = %g\n", j, x[j]);
			}
			for(i=0;i<m;i++) fvec[i]=wa4[i];
			xnorm=enorm(n,wa2);
			fnorm=fnorm1;iter++;
      	}
      	if(fabs(actred)<=ftol&&prered<=ftol&&0.5*ratio<=1.0) info=1;
      	if(delta<=xtol*xnorm) {info=2;}
      	if(fabs(actred)<=ftol&&prered<=ftol&&0.5*ratio<=1.0&&info==2) info=3;
      	if(info!=0) return(info);
      	if(fabs(actred)<=epsmch&&prered<=epsmch&&0.5*ratio<=1.0) {info=6;/*cout<<"\nactred="<<actred<<" prered="<<prered<<" ratio="<<ratio<<'\n';*/}
      if(delta<=epsmch*xnorm) {info=7;/*cout<<"\ndelta="<<delta<<'\n';*/}
      if(gnorm<=epsmch) {info=8;}
      if(info!=0) return(info);
    }while(ratio<0.0001);
  }while(1);
}

/** calls cminpack covar function and display covariant matrix on screen**/
void Lmdif::covar_cminpack()
{
    double fnorm = enorm(m, fvec);
    int i, j;
    printf("m = %d\n", m);
    printf("n = %d\n", n);
    printf("ldfjac = %d\n", ldfjac);
    for (i = 0; i < n; ++i)
    	printf("x[%d] = %g\n", i, x[i]);
    printf("fnorm = %15.7g\n", fnorm);
    double covfac = fnorm*fnorm / (m-n);
    int counter = -1;
    for (j = 0; j < n; ++j) {
        for (i = 0; i < m; ++i) {
            fjac_cminpack[i + j*ldfjac] = fjac[j][i];
            counter++;
            //printf("%d %g %g\n", counter, fjac[j][i], fjac_cminpack[i + j*ldfjac]);
        }
    }
    /*
    printf("Finished copying fjac to fjac_cminpack\n");
    printf("fjac_cminpack[0] = %g\n", fjac_cminpack[0]);
    printf("fjac_cminpack[1] = %g\n", fjac_cminpack[1]);
    printf("fjac_cminpack[%d] = %g\n", ldfjac, fjac_cminpack[ldfjac]);
    printf("fjac_cminpack[%d + 1] = %g\n", ldfjac, fjac_cminpack[ldfjac+1]);
    */
    /** increment ipvt[i] by 1 to match the convention used in cminpack; otherwise, covar causes segfault **/
    for (i = 0; i < n; ++i)
    	ipvt[i] = ipvt[i]+1;
    covar(n, fjac_cminpack, ldfjac, ipvt, ftol, wa1);
    printf("covariance (using covar)\n");
    for (i = 0; i < n; ++i) {
        for (j = 0; j < n; ++j)
            printf("%s%15.7g", j%n==0?"\n   ":"", (double)fjac_cminpack[i*ldfjac + j]*covfac);
    }
    printf("\n");
}
