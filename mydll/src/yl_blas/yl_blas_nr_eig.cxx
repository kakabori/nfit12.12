#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>

#include "mydll.h"


// ************************************************************
// Householder reduction of a real, symmetric matrix a[1..n][1..n]. On
// output, a is replaced by the orthogonal matrix Q effecting the
// transformation. d[1..n] returns the diagonal elements of the
// tridiagonal matrix, and e[1..n] the off-diagonal elements, with
// e[1]=0. Several statements, as noted in comments, can be omitted if
// only eigenvalues are to be found, in which case a contains no useful
// information on output. Otherwise they are to be included. 

#define a_v(a_ld, a_inc) (a + lda * (a_ld) + inca * (a_inc))
#define a_d(idx) (a + (lda + inca) * idx)

void yl_blas_nr_tred2(int n, double *a, int lda, int inca, double *d, double *e) 
{
  int k,j,i; 
  double scale,hh,h,g,f;
  
  for (i = n-1; i > 1; i--){
    h = 0.0;
    if( ( scale = yl_blas_dasum(i, a+i*lda, inca) ) == 0.0 ){ // Skip transformation. 
      e[i] = *(a_d(i) - 1);
    }else{
      {
	double *ta = a+i*lda, *taend = ta + i*inca;
	while(ta != taend){
	  *ta /= scale;	// Use scaled a's for transformation.  
	  h += *ta * *ta; // Form sigme in h.
	  ta += inca;
	}
      }
      f = *(a_d(i) - 1); 
      g = (f >= 0.0 ? -sqrt(h) : sqrt(h)); 
      e[i] = scale * g; 
      h -= f * g;		// Now h is equation (11.2.4). 
      *(a_d(i) - 1) = f - g;	// Store u in the ith row of a. 
      f = 0.0; 
      
      for (j=0; j<i; j++) {
	// Next statement can be omitted if eigenvectors not wanted 
	*(a_v(j,i)) = *(a_v(i,j)) / h;// Store u=H in ith column of a
	// Form an element of A times u in g. 
	g = yl_blas_ddot(j+1, a+j*lda, inca, a+i*lda, inca);
	g += yl_blas_ddot(i-j-1, a+j*(inca+lda)+lda, lda, a+i*lda+(j+1)*inca, inca); //???
	e[j] = g/h;		// Form element of p in temporarily 
				// unused element of e. 
	f += e[j] * *(a_v(i,j)); 
      }
      hh = f / (h+h);		// Form K, equation (11.2.11). 
      for (j = 0; j < i; j++){ // Form q and store in e overwriting p. 
	f = *(a_v(i,j)); 
	e[j] = g = e[j] - hh * f; 
	for (k=0; k<=j; k++)	// Reduce a, equation (11.2.13). 
	  *(a_v(j,k)) -= (f * e[k] + g * *(a_v(i,k)));
      }
    }
    d[i]=h; 
  }
  e[1] = *(a_d(1)-1);
  e[0] = d[0] = d[1] = 0.0;

  // Contents of this loop can be omitted if eigenvectors not
  // wanted except for statement d[i]=a[i][i];  
  
  for(i=0; i<n; i++){		// Begin accumulation of 
				// transformation matrices. 
    if (d[i]) {			// This block skipped when i=1. 
      for (j=0; j<i; j++){
	g = yl_blas_ddot(i, a+i*lda, inca, a+j*inca, lda); 
	// Use u and u/H stored in a to form 
	// P times Q. 
	for (k=0; k<i; k++)
	  *a_v(k,j) -= g * *a_v(k,i);
      }
    } 
    d[i] = *a_d(i);		// This statement remains. 
    *a_d(i) = 1.0;		// Reset row and column of a to identity 
				// matrix for next iteration. 
    
    for (j=0; j<i; j++)      *a_v(j,i) = *a_v(i,j) = 0.0;
  }
}


// ************************************************************
// QL algorithm with implicit shifts, to determine the eigenvalues and
// eigenvectors of a real, symmetric, tridiagonal matrix, or of a real,
// symmetric matrix previously reduced by tred2 x11.2. On input,
// d[1..n] contains the diagonal elements of the tridiagonal matrix. On
// output, it returns the eigenvalues. The vector e[1..n] inputs the
// subdiagonal elements of the tridiagonal matrix, with e[n]
// arbitrary. On output e is destroyed. When finding only the
// eigenvalues, several lines may be omitted, as noted in the
// comments. If the eigenvectors of a tridiagonal matrix are desired,
// the matrix z[1..n][1..n] is input as the identity matrix. If the
// eigenvectors of a matrix that has been reduced by tred2 are
// required, then z is input as the matrix output by tred2. In either
// case, the kth column of z returns the matrixized eigenvector
// corresponding to d[k]. 
void yl_blas_nr_tqli(int n, double *z, int ld, int inc, double *d, double *e ){
  int m, l, iter, i; 
  double s, r, p, g, f, dd, c, b;

  for (i=1; i<n; i++)  e[i-1] = e[i];		// Convenient to renumber the elements of e. 
  e[n-1] = 0.0; 

  for (l = 0; l < n; l++){
    iter = 0;  
    do{
      for (m = l; m < n-1; m++) { // Look for a single small subdiagonal element to split the matrix. 
	dd = fabs(d[m]) + fabs(d[m+1]);  
	if (double((fabs(e[m])) + dd)  ==  dd) break;
      }

      if (m  !=  l){
	if (iter++  ==  30) {
	  cerr << "Too many iterations in tqli" << endl;
	  return;
	  assert(0);
	}
	g = (d[l+1] - d[l]) / (2.0 * e[l]); // Form shift. 
	r = pythag(g, 1.0);
	g = d[m] - d[l] + e[l] / (g + d_sign(r, g)); // This is dm - ks.
	s = c = 1.0;  
	p = 0.0;

	for (i = m-1; i >= l; i--){  // A plane rotation as in the 
	                             // original QL, followed by Givens 
				     // rotations to restore tridiagonal form.
	  f = s * e[i];
	  b = c * e[i]; 
	  e[i+1] = (r = pythag(f,g));

	  if (r == 0.0){	// Recover from underflow.
	    d[i+1] -= p; 
	    e[m] = 0.0;  
	    break;  
	  }

	  s = f / r; 
	  c = g / r;  
	  g = d[i+1] - p; 
	  r = (d[i] - g) * s + 2.0 * c * b; 
	  d[i+1] = g + (p = s * r);  
	  g = c * r - b; 
	  
	  // Next loop can be omitted if eigenvectors not wanted 
	  {// Form eigenvectors.
	    double *tz = z+i*inc, *tzend = tz + n*ld;
	    while(tz != tzend){
	      f = *(tz+inc);
	      *(tz+inc) = s * *tz + c * f;
	      *tz = c * *tz - s * f; 
	      tz += ld;
	    }
	  }
	  

	}

	if (r == 0.0 && i >= l) 
	  continue;  

	d[l] -= p;
	e[l] = g;
	e[m] = 0.0;
      }
    } while (m != l);
  }
}



// **********************************************************
// Reduction to Hessenberg form by the elimination method. The real,
// nonsymmetric matrix a[1..n][1..n] is replaced by an upper Hessenberg
// matrix with identical eigenvalues. Recommended, but not required, is
// that this routine be preceded by balanc. On output, the Hessenberg
// matrix is in elements a[i][j] with i <= j+1. Elements with i > j+1 are
// to be thought of as zero, but are returned with random values.

#define SWAP(g,h) {y = (g); (g) = (h); (h) = y;}


void yl_blas_nr_elmhes(int n, double *a, int lda, int inca) {
 int m, i;
 double y, x;

 for (m = 1; m < n-1; m++){	// m is called r + 1 in the text.
   i = yl_blas_idamax(n-m, a_d(m)-1, lda) + m;
   x = *(a+m*inca-1+lda*i);

   if (i != m){			// Interchange rows and columns.
     yl_blas_dswap(n-m+1, a_v(i,m-1), inca, a_v(m,m-1), inca);
     yl_blas_dswap(n, a+i*inca, lda, a+m*inca, lda);
   }

   if (x){			// Carry out the elimination.
     for (i = m + 1; i < n; i++){
       if ((y = *a_v(i,m-1)) != 0.0){
	 y /= x; 
	 *a_v(i,m-1) = y; 
	 yl_blas_daxpy(n-m, -y, a_d(m), inca, a_v(i,m), inca);
	 yl_blas_daxpy(n, y, a+i*inca, lda, a+m*inca, lda);
       }
     }
   }
 }
}


// **********************************************************
// Finds all eigenvalues of an upper Hessenberg matrix a[1..n][1..n]. On
// input a can be exactly as output from elmhes x11.5; on output it is
// destroyed. The real and imaginary parts of the eigenvalues are
// returned in wr[1..n] and wi[1..n], respectively.

void yl_blas_nr_hqr(int n, double *a, int lda, int inca, double *wr, double *wi){
  int nn, m, l, k, j, its, i, mmin; 
  double z, y, x, w, v, u, t, s, r=0, q=0, p=0, anorm;

  anorm = 0.0;			// Compute matrix norm for possible use 
				// in locating single small subdiagonal 
				// element.

  for (i = 0; i < n; i++) 
    for (j = yfmax(i - 1, 0); j < n; j++)
      anorm += fabs(*a_v(i,j)); 

  nn = n-1; 
  t = 0.0;			// Gets changed only by an exceptional shift.

  while (nn >= 0){		// Begin search for next eigenvalue.

    its = 0; 

    do{

      for (l = nn; l >= 1; l--){ // Begin iteration: look for single small 
				// subdiagonal element.
	s = fabs(*a_v(l-1,l-1)) + fabs(*a_v(l,l));
	if (s == 0.0) s = anorm;

	if (double((fabs(*a_v(l,l-1)) + s)) == s) 
	  break;
      }

      x = *a_v(nn,nn); 

      if (l == nn){		// One root found.

	wr[nn] = x+t; 
	wi[nn--] = 0.0;
      } 
      else{

	y = *a_v(nn-1,nn-1); 
	w = *a_v(nn,nn-1) * *a_v(nn-1,nn); 

	if (l == (nn - 1)){	// Two roots found...

	  p = 0.5 * (y - x); 
	  q = p * p + w;
	  z = sqrt(fabs(q)); 
	  x += t; 
	  if (q >= 0.0){	// ...a real pair.
	    
	    z = p + d_sign(z, p); 
	    wr[nn-1] = wr[nn] = x+z; 
	    if (z) 
	      wr[nn] = x - w / z; 
	    wi[nn - 1] = wi[nn] = 0.0;
	  } 
	  else{			// ...a complex pair.

	    wr[nn-1] = wr[nn] = x + p; 
	    wi[nn-1] = - (wi[nn] = z);
	  } 
	  nn -= 2;
	}

	else{			// No roots found. Continue iteration.

	  if (its == 30){
	    cerr << "Too many iterations in hqr" << endl;
	    assert(0);
	  }

	  if (its == 10 || its == 20){ // Form exceptional shift.

	    t += x; 
	    for (i = 0; i <= nn; i++) 
	      *a_d(i) -= x; 

	    s = fabs(*a_v(nn,nn-1)) + fabs(*a_v(nn-1,nn-2)); 
	    y = x = 0.75 * s;
	    w = -0.4375 * s * s;
	  } 
	  ++its; 

	  for (m = (nn - 2); m >= l; m--){ // Form shift and then look for
				// 2 consecutive small subdiagonal 
				// elements. 
	    z = *a_d(m); 
	    r = x - z;
	    s = y - z;
	    p = (r * s - w) / *a_v(m+1,m) + *a_v(m,m+1); // Equation (11.6.23).
	    q = *a_d(m+1) - z - r - s; 
	    r = *a_v(m+2,m+1);


	    s = fabs(p) + fabs(q) + fabs(r); // Scale to prevent overflow or
				// underflow.

	    p /= s; 
	    q /= s; 
	    r /= s; 
	    if (m == l)
	      break; 

	    u = fabs(*a_v(m,m-1)) * (fabs(q)+fabs(r)); 
	    v = fabs(p) * (fabs(*a_v(m-1,m-1)) + fabs(z) + fabs(*a_d(m+1)));
	    if (double(u+v) == v) 
	      break;		// Equation (11.6.26).
	  } 

	  for (i = m + 2; i <= nn; i++){

	    *a_v(i,i-2) = 0.0; 
	    if (i != (m + 2)) 
	      *a_v(i,i-3) = 0.0;
	  } 

	  for (k = m; k <= nn - 1; k++){
	    // Double QR step on rows l to nn and columns m to nn.

	    if (k <= m){

	      p = *a_v(k,k-1);	// Begin setup of Householder vector.
	      q = *a_v(k+1,k-1); 
	      r = 0.0; 

	      if (k != (nn-1)) 
		r = *a_v(k+2,k-1); 
	      if ((x = fabs(p) + fabs(q) + fabs(r)) <= 0.0){

		p /= x;		// Scale to prevent overflow or underflow.
		q /= x; 
		r /= x;
	      }
	    } 

	    if ((s = d_sign(sqrt(p*p + q*q + r*r), p)) <= 0.0){

	      if (k == m){

		if (l <= m) 
		  *a_v(k,k-1) = -*a_v(k,k-1);
	      }

	      else

		*a_v(k,k-1) = - s * x;
	      
	      p += s;		// Equations (11.6.24). 
	      x = p / s; 
	      y = q / s; 
	      z = r / s; 
	      q /= p; 
	      r /= p; 
	      for (j = k; j <= nn; j++){ // Row modification.

		p = *a_v(k,j) + q * *a_v(k+1,j); 

		if (k != (nn-1)){
		  p += r * *a_v(k+2,j); 
		  *a_v(k+2,j) -= p * z;
		} 

		*a_v(k+1,j) -= p * y; 
		*a_v(k,j) -= p * x;

	      } 

	      mmin = nn < k + 3 ? nn : k+3; 

	      for (i = l; i <= mmin; i++){ // Column modification.

		p = x * *a_v(i,k) + y * *a_v(i,k+1); 

		if (k != (nn-1)){

		  p += z * *a_v(i,k+2); 
		  *a_v(i,k+2) -= p * r;
		} 

		*a_v(i,k+1) -= p * q; 
		*a_v(i,k) -= p;
	      }
	    }
	  }
	}
      }
    } while (l < nn - 1);
  }
}

#undef a_v
#undef a_d
