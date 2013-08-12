/***********************************************************
code borrowed from numerical recipe
look for the chapter that explains SVD 
************************************************************/


#include <math.h>
#include <stdlib.h>
#define dfprec float

static dfprec at,bt,ct;
#define PYTHAG(a,b) ((at=fabs(a)) > (bt=fabs(b)) ? \
(ct=bt/at,at*sqrt(1.0+ct*ct)) : (bt ? (ct=at/bt,bt*sqrt(1.0+ct*ct)): 0.0))

static dfprec maxarg1,maxarg2;
#define MAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ?\
	(maxarg1) : (maxarg2))
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
void svdcmp(dfprec *a,int m,int n,dfprec *w,dfprec *v)
{
  int flag,i,its,j,jj,k,l,nm;
  dfprec c,f,h,s,x,y,z;
  dfprec anorm=0.0,g=0.0,scale=0.0;
  dfprec *rv1;
  rv1=(dfprec *)malloc(sizeof(dfprec)*n);
  l = 1;
  for (i=0;i<n;i++) {
    dfprec *ai;
    ai=a+i*m;
    l=i+1;
    rv1[i]=scale*g;
    g=s=scale=0.0;
    if (i < m) {
      for (k=i;k<m;k++) scale += fabs(ai[k]);
      if (scale) {
	for (k=i;k<m;k++) {
	  ai[k] /= scale;
	  s += ai[k]*ai[k];
	}
	f=ai[i];
	g = -SIGN(sqrt(s),f);
	h=f*g-s;
	ai[i]=f-g;
	if (i != n-1) {
	  for (j=l;j<n;j++) {
	    dfprec *aj;
	    aj=a+j*m;
	    for (s=0.0,k=i;k<m;k++) s += ai[k]*aj[k];
	    f=s/h;
	    for (k=i;k<m;k++) aj[k] += f*ai[k];
	  }
	}
	for (k=i;k<m;k++) ai[k] *= scale;
      }
    }
    w[i]=scale*g;
    g=s=scale=0.0;
    if (i < m && i != n-1) {
      dfprec *atmp;
      atmp=a+l*m;
      for (k=l;k<n;k++,atmp+=m) scale += fabs(atmp[i]);
      if (scale) {
	atmp=a+l*m;
	for (k=l;k<n;k++,atmp+=m) {
	  atmp[i] /= scale;
	  s += atmp[i]*atmp[i];
	}
	atmp=a+l*m;
	f=atmp[i];
	g = -SIGN(sqrt(s),f);
	h=f*g-s;
	atmp[i]=f-g;
	for (k=l;k<n;k++,atmp+=m) rv1[k]=atmp[i]/h;
	if (i != m-1) {
	  for (j=l;j<m;j++) {
	    atmp=a+l*m;
	    for (s=0.0,k=l;k<n;k++,atmp+=m) s += atmp[j]*atmp[i];
	    atmp=a+l*m;
	    for (k=l;k<n;k++,atmp+=m) atmp[j] += s*rv1[k];
	  }
	}
	atmp=a+l*m;
	for (k=l;k<n;k++,atmp+=m) atmp[i] *= scale;
      }
    }
    anorm=MAX(anorm,(fabs(w[i])+fabs(rv1[i])));
  }
  for (i=n-1;i>=0;i--) {
    dfprec *vi;
    vi=v+i*n;
    if (i < n-1) {
      if (g) {
	for (j=l;j<n;j++)
	  vi[j]=((a+m*j)[i]/(a+m*l)[i])/g;
	for (j=l;j<n;j++) {
	  for (s=0.0,k=l;k<n;k++) s += (a+m*k)[i]*(v+n*j)[k];
	  for (k=l;k<n;k++) (v+n*j)[k] += s*vi[k];
	}
      }
      for (j=l;j<n;j++) (v+n*j)[i]=vi[j]=0.0;
    }
    vi[i]=1.0;
    g=rv1[i];
    l=i;
  }
  for (i=n-1;i>=0;i--) {
    dfprec *ai;
    ai=a+i*m;
    l=i+1;
    g=w[i];
    if (i < n-1)
      for (j=l;j<n;j++) (a+m*j)[i]=0.0;
    if (g) {
      g=1.0/g;
      if (i != n-1) {
	for (j=l;j<n;j++) {
	  for (s=0.0,k=l;k<m;k++) s += ai[k]*(a+m*j)[k];
	  f=(s/ai[i])*g;
	  for (k=i;k<m;k++) (a+m*j)[k] += f*ai[k];
	}
      }
      for (j=i;j<m;j++) ai[j] *= g;
    } else {
      for (j=i;j<m;j++) ai[j]=0.0;
    }
    ++ai[i];
  }
  for (k=n-1;k>=0;k--) {
    dfprec *vk;
    vk=v+k*n;
    for (its=1;its<=30;its++) {
      flag=1;
      for (l=k;l>=0;l--) {
	nm=l-1;
	if (fabs(rv1[l])+anorm == anorm) {
	  flag=0;
	  break;
	}
	if (fabs(w[nm])+anorm == anorm) break;
      }
      if (flag) {
	c=0.0;
	s=1.0;
	for (i=l;i<=k;i++) {
	  f=s*rv1[i];
	  if (fabs(f)+anorm != anorm) {
	    dfprec *anm, *ai;
	    anm=a+nm*m;
	    ai=a+m*i;
	    g=w[i];
	    h=PYTHAG(f,g);
	    w[i]=h;
	    h=1.0/h;
	    c=g*h;
	    s=(-f*h);
	    for (j=0;j<m;j++) {
	      y=anm[j];
	      z=ai[j];
	      anm[j]=y*c+z*s;
	      ai[j]=z*c-y*s;
	    }
	  }
	}
      }
      z=w[k];
      if (l == k) {
	if (z < 0.0) {
	  w[k] = -z;
	  for (j=0;j<n;j++) vk[j]=(-vk[j]);
	}
	break;
      }
      x=w[l];
      nm=k-1;
      y=w[nm];
      g=rv1[nm];
      h=rv1[k];
      f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
      g=PYTHAG(f,1.0);
      f=((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;
      c=s=1.0;
      for (j=l;j<=nm;j++) {
	dfprec *vj, *vi;
	vj=v+j*n;
	vi=vj+n;
	i=j+1;
	g=rv1[i];
	y=w[i];
	h=s*g;
	g=c*g;
	z=PYTHAG(f,h);
	rv1[j]=z;
	c=f/z;
	s=h/z;
	f=x*c+g*s;
	g=g*c-x*s;
	h=y*s;
	y=y*c;
	for (jj=0;jj<n;jj++) {
	  x=vj[jj];
	  z=vi[jj];
	  vj[jj]=x*c+z*s;
	  vi[jj]=z*c-x*s;
	}
	z=PYTHAG(f,h);
	w[j]=z;
	if (z) {
	  z=1.0/z;
	  c=f*z;
	  s=h*z;
	}
	f=(c*g)+(s*y);
	x=(c*y)-(s*g);
	vj=a+j*m;//use vj for aj
	vi=vj+m;//use vi for ai
	for (jj=0;jj<m;jj++) {
	  y=vj[jj];
	  z=vi[jj];
	  vj[jj]=y*c+z*s;
	  vi[jj]=z*c-y*s;
	}
      }
      rv1[l]=0.0;
      rv1[k]=f;
      w[k]=x;
    }
  }
  free(rv1);
}

#undef SIGN
#undef MAX
#undef PYTHAG

void svbksb(dfprec *u,dfprec *w,dfprec *v,int m,int n,dfprec *b,dfprec *x)
{
  int jj,j,i;
  dfprec s,*tmp;
  tmp=(dfprec *)malloc(sizeof(dfprec)*n);
  for (j=0;j<n;j++) {
    dfprec *uj;
    uj=u+j*m;
    s=0.0;
    if (w[j]) {
      for (i=0;i<m;i++) s += uj[i]*b[i];
      s /= w[j];
    }
    tmp[j]=s;
  }
  for (j=0;j<n;j++) {
    dfprec *vj;
    vj=v+j;
    s=0.0;
    for (jj=0;jj<n;jj++,vj+=n) s += (*vj)*tmp[jj];
    x[j]=s;
  }
  free(tmp);
}

#define TOL 1.0e-5

void linfit_pre(dfprec *x, dfprec *y, dfprec *sig, int ndata, int ma, int xn, dfprec *u, dfprec *b,void (*funcs)(dfprec *,dfprec *,int)){
  int i,j;
  dfprec *afunc;
  dfprec tmp;
  afunc=(dfprec *)malloc(sizeof(dfprec)*ma);
  for (i=0;i<ndata;i++) {
    (*funcs)(x+i*xn,afunc,ma);
    tmp=1.0/sig[i];
    for (j=0;j<ma;j++) (u+j*ndata)[i]=afunc[j]*tmp;
    b[i]=y[i]*tmp;
  }
  free(afunc);
}

dfprec linfit_post(dfprec *x, dfprec *y, dfprec *sig, int ndata, int ma, int xn, dfprec *a, void (*funcs)(dfprec *,dfprec *,int)){
  dfprec chisq=0.0, sum, *afunc, tmp;
  int i,j;
  afunc=(dfprec *)malloc(sizeof(dfprec)*ma);
  for (i=0;i<ndata;i++) {
    (*funcs)(x+i*xn,afunc,ma);
    for (sum=0.0,j=0;j<ma;j++) sum += a[j]*afunc[j];
    chisq += (tmp=(y[i]-sum)/sig[i],tmp*tmp);
  }
  free(afunc);
  return chisq;
}
void svdfit(int ndata,dfprec *a,int ma,dfprec *u,dfprec *v,dfprec *w, dfprec *b){
  int j;
  dfprec wmax,thresh;
  svdcmp(u,ndata,ma,w,v);
  wmax=0.0;
  for (j=0;j<ma;j++)
    if (w[j] > wmax) wmax=w[j];
  thresh=TOL*wmax;
  for (j=0;j<ma;j++)
    if (w[j] < thresh) w[j]=0.0;
  svbksb(u,w,v,ndata,ma,b,a);
}
#undef TOL



