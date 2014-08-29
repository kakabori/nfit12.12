#include <math.h>
#include <stdlib.h>
#include "mydll.h"

double pythag(double a, double b){
  double at,bt,ct;
  return ((at=fabs(a)) > (bt=fabs(b)) ? (ct=bt/at,at*sqrt(1.0+ct*ct)) : (bt ? (ct=at/bt,bt*sqrt(1.0+ct*ct)): 0.0));
}


#define a_ptr(a_ld, a_inc) (a + lda * (a_ld) + inca * (a_inc))
#define aic_ptr(aic_ld) ( aic + (aic_ld) * lda )

#define w_ptr(w_inc) ( w + (w_inc) * incw )

#define v_ptr(v_ld, v_inc) (v + ldv * (v_ld) + incv * (v_inc))
#define vir_ptr(vir_inc) ( vir + (vir_inc) * incv )


/************************************************

    A           =         U         W           VT

    n                            |            ||||||
  ||||||               ||||||     |           ||||||
  ||||||               ||||||      |          ||||||
m ||||||        =      ||||||       |         ||||||
  ||||||               ||||||        |        ||||||
                                      |       ||||||

*************************************************/

void yl_cblas_nr_svdcmp( int m, int n,
			 double *a, int lda, int inca,
			 double *w, int incw,
			 double *v, int ldv, int incv)
{
  int flag,i,its,j,k,l,nm;
  double c,f,h,s,x,y,z;
  double anorm=0.0,g=0.0,scale=0.0;
  double *rv1;
  rv1=(double *)malloc(sizeof(double)*n);
  for (i=0;i<n;i++) {
    double *aii;
    aii = a_ptr(i,i);
    l = i+1;
    rv1[i] = scale*g;
    g=s=scale=0.0;
    if (i < m) {
      scale += yl_blas_dasum(m-i, aii, inca);
      if (scale) {
	double *aend, *acur;
	acur = aii;
	aend = aii+(m-i)*inca;
	while (acur != aend) {
	  *acur /= scale;
	  s += *acur * *acur;
	  acur += inca;
	}
	f = *aii;
	g = -d_sign(sqrt(s),f);
	h = f*g-s;
	*aii = f-g;
	if (i != n-1) {
	  for (j=l;j<n;j++) {
	    double *aji;
	    aji = a_ptr(j,i);
	    s = yl_blas_ddot( m-i, aii, inca, aji, inca);
	    f=s/h;
	    yl_blas_daxpy( m-i, f, aii, inca, aji, inca);
	  }
	}
	yl_blas_dscal(m-i, scale, aii, inca);
      }
    }
    *w_ptr(i) = scale*g;
    g=s=scale=0.0;
    if (i < m && i != n-1) {
      double *ali;
      ali = a_ptr(l, i);
      scale += yl_blas_dasum( n-l, ali, lda );
      if (scale) {
	double *aend, *acur;
	acur = ali;
	aend = ali + (n - l)*lda;
	while (acur != aend) {
	  *acur /= scale;
	  s += *acur * *acur;
	  acur += lda;
	}
	f = *ali;
	g = -d_sign(sqrt(s),f);
	h = f*g-s;
	*ali = f-g;
	yl_blas_daxty( n-l, 1/h, ali, lda, rv1+l, 1);
	if (i != m-1) {
	  for (j=l;j<m;j++) {
	    double *alj;
	    alj = a_ptr(l,j);
	    s = yl_blas_ddot(n-l, ali, lda, alj, lda);
	    yl_blas_daxpy( n-l, s, rv1+l, 1, alj, lda);
	  }
	}
	yl_blas_dscal( n-l, scale, ali, lda );
      }
    }
    anorm = yfmax(anorm,(fabs(*w_ptr(i))+fabs(rv1[i])));
  }
  for (i=n-1;i>=0;i--) {
    double *vir;
    double *aic;

    vir = v + i*ldv;
    aic = a + i*inca;
    if (i < n-1) {
      if (g) {
	//can be better?
	for (j=l;j<n;j++){
	  *vir_ptr(j) = *aic_ptr(j) / *aic_ptr(l) / g;
	}
	for (j=l;j<n;j++) {
	  s = yl_blas_ddot(n-l, a_ptr(l,i), lda, v_ptr(j,l), incv);
	  yl_blas_daxpy( n-l, s, vir_ptr(l), incv, v_ptr(j, l), incv);
	}
      }
      yl_blas_dset(n-l, 0, vir_ptr(l), incv);
      yl_blas_dset(n-l, 0, v_ptr(l,i), ldv);
    }
    *vir_ptr(i) = 1.0;
    g=rv1[i];
    l=i;
  }
  for (i=yfmin(n,m)-1;i>=0;i--) {
    double *aii;
    aii = a_ptr(i,i);
    l=i+1;
    g = *w_ptr(i);
    yl_blas_dset(n-l, 0, a_ptr(l,i), lda);
    if (g) {
      g=1.0/g;
      if (i != n-1) {
	for (j=l;j<n;j++) {
	  s = yl_blas_ddot( m-l, a_ptr(i,l), inca, a_ptr(j,l), inca);
	  f = s / *aii * g;
	  yl_blas_daxpy( m-i, f, aii, inca, a_ptr(j,i), inca);
	}
      }
      yl_blas_dscal(m-i, g, aii, inca);
    } else {
      yl_blas_dset(m-i, 0, aii, inca);
    }
    ++(*aii);
  }
  for (k=n-1;k>=0;k--) {
    for (its=1;its<=30;its++) {
      flag=1;
      for (l=k;l>=0;l--) {
	nm=l-1;
	if (fabs(rv1[l])+anorm == anorm) {
	  flag=0;
	  break;
	}
	if (fabs( *w_ptr(nm) ) + anorm == anorm) break;
      }
      if (flag) {
	c=0.0;
	s=1.0;
	for (i=l;i<=k;i++) {
	  f = s*rv1[i];
	  if (fabs(f)+anorm != anorm) {
	    double *anm, *ai, *aend;
	    g = *w_ptr(i);
	    h = pythag(f,g);
	    *w_ptr(i) = h;
	    h = 1.0/h;
	    c = g*h;
	    s = (-f*h);

	    anm = a + nm * lda;
	    ai  = a + i  * lda;
	    aend = anm + m *inca;
	    while ( anm != aend ) {
	      y = *anm;
	      z = *ai;
	      *anm = y*c+z*s;
	      *ai  = z*c-y*s;
	      anm += inca;
	      ai  += inca;
	    }
	  }
	}
      }
      z = *w_ptr(k);
      if (l == k) {
	if (z < 0.0) {
	  *w_ptr(k) = -z;
	  yl_blas_dscal( n, -1, v+k*ldv, incv );
	}
	break;
      }
      x = *w_ptr(l);
      nm= k-1;
      y = *w_ptr(nm);
      g = rv1[nm];
      h = rv1[k];
      f = ((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
      g = pythag(f,1.0);
      f = ((x-z)*(x+z)+h*((y/(f+d_sign(g,f)))-h))/x;
      c = s = 1.0;
      for (j=l;j<=nm;j++) {
	double *vj, *vi, *vend;
	i  = j+1;
	g  = rv1[i];
	y  = *w_ptr(i);
	h  = s*g;
	g  = c*g;
	z  = pythag(f,h);
	rv1[j]=z;
	c  = f/z;
	s  = h/z;
	f  = x*c+g*s;
	g  = g*c-x*s;
	h  = y*s;
	y  = y*c;

	vj = v+j*ldv;
	vi = vj+ldv;
	vend = vj + n * incv;
	while ( vj != vend) {
	  x = *vj;
	  z = *vi;
	  *vj = x*c+z*s;
	  *vi = z*c-x*s;
	  vj += incv;
	  vi += incv;
	}
	z=pythag(f,h);
	*w_ptr(j) = z;
	if (z) {
	  z=1.0/z;
	  c=f*z;
	  s=h*z;
	}
	f=(c*g)+(s*y);
	x=(c*y)-(s*g);
	vj = a+j*lda;//use vj for aj
	vi = vj+lda;//use vi for ai
	vend = vj + m * inca; // use vend for aend
	while (vj != vend) {
	  y = *vj;
	  z = *vi;
	  *vj = y*c+z*s;
	  *vi = z*c-y*s;
	  vj += inca;
	  vi += inca;
	}
      }
      rv1[l]=0.0;
      rv1[k]=f;
      *w_ptr(k)=x;
    }
  }
  free(rv1);
}


/************************************
    A    x    X  =  B

    n         |
  ||||||      |     |
  ||||||      |     |
m |||||| x    |  =  |
  ||||||      |     |
              |

only x is written on output
**************************************/


void yl_cblas_nr_svbksb(int m,int n,
			double *u, int ldu, int incu,
			double *w, int incw,
			double *v, int ldv, int incv,
			double *b, int incb,
			double *x, int incx)
{
  int j;
  double s,*tmp;
  tmp=(double *)malloc(sizeof(double)*n);
  for (j=0;j<n;j++) {
    double *uj;
    uj=u+j*ldu;
    s=0.0;
    if ( *w_ptr(j) ) {
      s += yl_blas_ddot(m, uj, incu, b, incb);
      s /= *w_ptr(j);
    }
    tmp[j]=s;
  }
  for (j=0;j<n;j++) {
    double *vj;
    vj=v+j*incv;
    s=0.0;
    s += yl_blas_ddot(n, tmp, 1, vj, ldv);
    *(x+j*incx) = s;
  }
  free(tmp);
}






void yl_cblas_nr_svdcmp(double *a,int m,int n,double *w,double *v)
{
  int flag,i,its,j,jj,k,l,nm;
  double c,f,h,s,x,y,z;
  double anorm=0.0,g=0.0,scale=0.0;
  double *rv1;
  rv1=(double *)malloc(sizeof(double)*n);
  for (i=0;i<n;i++) {
    double *ai;
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
	g = -d_sign(sqrt(s),f);
	h=f*g-s;
	ai[i]=f-g;
	if (i != n-1) {
	  for (j=l;j<n;j++) {
	    double *aj;
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
      double *atmp;
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
	g = -d_sign(sqrt(s),f);
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
    anorm = yfmax(anorm,(fabs(w[i])+fabs(rv1[i])));
  }
  for (i=n-1;i>=0;i--) {
    double *vi;
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
  for (i=yfmin(n,m)-1;i>=0;i--) {
    double *ai;
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
    double *vk;
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
	    double *anm, *ai;
	    anm=a+nm*m;
	    ai=a+m*i;
	    g=w[i];
	    h=pythag(f,g);
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
      g=pythag(f,1.0);
      f=((x-z)*(x+z)+h*((y/(f+d_sign(g,f)))-h))/x;
      c=s=1.0;
      for (j=l;j<=nm;j++) {
	double *vj, *vi;
	vj=v+j*n;
	vi=vj+n;
	i=j+1;
	g=rv1[i];
	y=w[i];
	h=s*g;
	g=c*g;
	z=pythag(f,h);
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
	z=pythag(f,h);
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

void yl_cblas_nr_svbksb(double *u,double *w,double *v,int m,int n,double *b,double *x)
{
  int jj,j,i;
  double s,*tmp;
  tmp=(double *)malloc(sizeof(double)*n);
  for (j=0;j<n;j++) {
    double *uj;
    uj=u+j*m;
    s=0.0;
    if (w[j]) {
      for (i=0;i<m;i++) s += uj[i]*b[i];
      s /= w[j];
    }
    tmp[j]=s;
  }
  for (j=0;j<n;j++) {
    double *vj;
    vj=v+j;
    s=0.0;
    for (jj=0;jj<n;jj++,vj+=n) s += (*vj)*tmp[jj];
    x[j]=s;
  }
  free(tmp);
}
