#include "mydll.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "yl_subplex.h"

isubc isubc_1;

//#define isubc_1 isubc_

/*  Originally Coded by Tom Rowan */

int yl_subplex_calcc(int ns,       /*subspace dimension*/
                     double *s,    /*>=ns*(ns+3) used to store simplex*/
                     int ih,       /*index to vertex with highest function value*/
                     int inew,     /*index to new point*/
                     bool updatc,  /*true: update; false: calculate centroid from scratch*/
                     double *c)    /*centroid of the simplex without vertex with */
                                   /*highest function value. OUTPUT: new centroid*/
{
  int i;

  if (updatc) {
    if (ih == inew) return 0;
    for (i = 0; i < ns; ++i)
      c[i] += (s[i + inew * ns] - s[i + ih * ns]) / ns;
  } else {
    memset(c, 0, sizeof(double)*ns); /* memset is better than yl_blas_copy_ */
    for (i = 0; i <= ns; ++i)
      if (i != ih) yl_blas_daxpy(ns, 1, s+i*ns, 1, c, 1);
    yl_blas_dscal(ns, 1.0/ns, c, 1);
  }
  return 0;
}


/* start creates the initial simplex for simplx minimization. */

void yl_subplex_start(int n,     /*problem dimension */
		       double *x,			  /*current best point */
		       double *step,		  /*stepsizes for corresponding components of x */
		       int ns,				  /*subspace dimension */
		       int *ips,			  /*permutation vector */
		       double *s,			  /*OUTPUT: first ns+1 columns contain initial simplex */
		       bool *small)			  /*OUTPUT:logical flag */
  /*            = .true.  : coincident points */
  /*            = .false. : otherwise */

{
  int j;

  for (j = 0; j < ns; ++j)	s[j] = x[ips[j]];

  for (j = 0; j < ns; ++j) {
    double *tp;
    tp = s + (j+1)*ns;
    yl_blas_dcopy(ns, s, 1, tp, 1);
    tp[j] = s[j] + step[ips[j]];
  }

  /* check for coincident points */

  for (j = 0; j < ns; ++j) {
    if (s[j + (j+1)*ns] == s[j]) {
      *small = true;
      return;
    }
  }
  *small = false;
}

/*checked*/

/* evalf evaluates the function f at a point defined by x */
/* with ns of its components replaced by those in xs. */

int yl_subplex_evalf(
		      YLMinimizationFunc f,
		      int ns,					/* subspace dimension */
		      int *ips,				/*permutation vector*/
		      double *xs,				/*double precision ns-vector to be mapped to x*/
		      int n,					/*problem dimension*/
		      double *x,				/*double precision n-vector*/
		      double *sfx,				/*OUTPUT: signed value of f evaluated at x */
		      int *nfe,				/*number of function evaluations OUTPUT: incremented*/
		      YLSUBPLEXctrl* ctrl)
{
  int i;
  double fx;
  bool newbst;
  /* mapping subspace ns-vector to x */
  for (i = 0; i < ns; ++i)    x[ips[i]] = xs[i];

  ctrl->newx = ( isubc_1.neww || ctrl->irepl != 2 );
  /* evaluate the function at the new postion x */
  fx = (*f)(x,ctrl->client,NULL);
  //  if(!ctrl->minf) fx = -fx; //add by YL

  if (ctrl->irepl == 0) {
    *sfx = fx;
  } else if (isubc_1.neww) {
    *sfx = fx;
    newbst = ( fx < ctrl->ftest );
    if (ctrl->initx || newbst) {
      if (ctrl->irepl == 1) {
	yl_subplex_fstats(fx, 1, true, ctrl);
      }
      ctrl->ftest = fx;
      isubc_1.sfbest = *sfx;
    }
  } else {
    if (ctrl->irepl == 1) {
      yl_subplex_fstats(fx, 1, false, ctrl);
      fx = ctrl->fxstat[ctrl->ifxsw - 1];
    }
    ctrl->ftest = fx + isubc_1.fbonus * ctrl->fxstat[3];
    *sfx = ctrl->ftest;
    isubc_1.sfbest = fx;
  }
  ++(*nfe);
  return 0;
} /* evalf_ */

/* fstats modifies the common /usubc/ variables nfxe,fxstat. */
/*   fx     - most recent evaluation of f at best x */
/*   ifxwt  - int weight for fx */
/*   reset  - bool switch */
/*            = .true.  : initialize nfxe,fxstat */
/*            = .false. : update nfxe,fxstat */
int yl_subplex_fstats(double fx,
		       int ifxwt,
		       bool reset,
		       YLSUBPLEXctrl* ctrl)
{
  double d1, d2, d3;
  double fscale;
  int nsv;
  double f1sv;

  if (reset) {
    ctrl->nfxe = ifxwt;
    ctrl->fxstat[0] = fx;
    ctrl->fxstat[1] = fx;
    ctrl->fxstat[2] = fx;
    ctrl->fxstat[3] = 0.;
  } else {
    nsv = ctrl->nfxe;
    f1sv = ctrl->fxstat[0];
    ctrl->nfxe += ifxwt;
    ctrl->fxstat[0] += ifxwt * (fx - ctrl->fxstat[0]) / ctrl->nfxe;
    ctrl->fxstat[1] = yfmax(ctrl->fxstat[1],fx);
    ctrl->fxstat[2] = yfmin(ctrl->fxstat[2],fx);
    fscale = yfmax( yfmax(fabs(ctrl->fxstat[1]), fabs(ctrl->fxstat[2])), 1.);
    d1 = ctrl->fxstat[3] / fscale;
    d2 = (ctrl->fxstat[0] - f1sv) / fscale;
    d3 = (fx - ctrl->fxstat[0]) / fscale;
    ctrl->fxstat[3] = fscale * sqrt(((nsv - 1) * (d1 * d1) + nsv * (
								      d2 * d2) + ifxwt * (d3 * d3)) / (ctrl->nfxe - 1));
  }
  return 0;
} /* fstats_ */


/* newpt performs reflections, expansions, contractions, and */
/* shrinkages (massive contractions) by computing: */
/* xbase + coef * (xbase - xold) */
/* The result is stored in xnew if neww .eq. .true., */
/* in xold otherwise. */
/* use :  coef .gt. 0 to reflect */
/*        coef .lt. 0 to expand, contract, or shrink */
/* input */
/*   ns     - number of components (subspace dimension) */
/*   coef   - one of four simplex method coefficients */
/*   xbase  - double precision ns-vector representing base point */
/*   xold   - double precision ns-vector representing old point */
/*   neww    - bool switch */
/*            = .true.  : store result in xnew */
/*            = .false. : store result in xold, xnew is not */
/*                        referenced */

/* output */
/*   xold   - unchanged if neww .eq. .true., contains neww */
/*            point otherwise */
/*   xnew   - double precision ns-vector representing neww */
/*            point if  neww .eq. .true., not referenced */
/*            otherwise */
/*   small  - bool flag */
/*            = .true.  : coincident points */
/*            = .false. : otherwise */

int yl_subplex_newpt(int ns,
		      double coef,
		      double *xbase,
		      double *xold,
		      bool neww,
		      double *xnew,
		      bool *small)
{
  int i;
  double xoldi;
  bool eqbase = true;
  bool eqold = true;

  if (neww) {
    for (i = 0; i < ns; ++i) {
      xnew[i] = xbase[i] + coef * (xbase[i] - xold[i]);
      eqbase = ( eqbase && xnew[i] == xbase[i] );
      eqold = ( eqold && xnew[i] == xold[i] );
    }
  } else {
    for (i = 0; i < ns; ++i) {
      xoldi = xold[i];
      xold[i] = xbase[i] + coef * (xbase[i] - xold[i]);
      eqbase = ( eqbase && xold[i] == xbase[i] );
      eqold = ( eqold && xold[i] == xoldi );
    }
  }
  *small = ( eqbase || eqold );
  return 0;
}

/* order determines the indices of the vertices with the */
/* lowest, second highest, and highest function values. */
/*   npts   - number of points in simplex */
/*   fs     - double precision vector of function values of simplex */
/*   il     - index to vertex with lowest function value */
/* output */
/*   il     - new index to vertex with lowest function value */
/*   is     - new index to vertex with second highest function value */
/*   ih     - new index to vertex with highest function value */

int yl_subplex_order(int npts,
		      double *fs,
		      int *il,
		      int *is,
		      int *ih)
{
  int iend;
  int i, j, il0;

  il0 = *il;
  j = (il0+1) % npts;
  if (fs[j] >= fs[*il]) {
    *ih = j;
    *is = il0;
  } else {
    *ih = il0;
    *il = *is = j;
  }
  iend = il0 + npts;
  for (i = il0 + 2; i < iend; ++i) {
    j = i % npts;
    if (fs[j] >= fs[*ih]) {
      *is = *ih;
      *ih = j;
    } else if (fs[j] > fs[*is]) *is = j;
    else if (fs[j] < fs[*il])	*il = j;
  }
  return 0;
}

/* setstp sets the stepsizes for the corresponding components of the solution vector. */
int yl_subplex_setstp(int *nsubs,/*   nsubs  - number of subspaces */
		       int n,/*   n      - number of components (problem dimension) */
		       double *deltax,/*   deltax - vector of change in solution vector */
		       double *step,/*   step   - new stepsizes */
		       YLSUBPLEXctrl* ctrl)
{
  int i;
  double stpfac;

  if (*nsubs > 1) {
    stpfac = yfmin(
		 yfmax(yl_blas_dasum(n, deltax, 1) / yl_blas_dasum(n, step, 1),
		     ctrl->omega),
		 1. / ctrl->omega);
  } else {
    stpfac = ctrl->psi;
  }
  yl_blas_dscal(n, stpfac, step, 1);

  /*     reorient simplex */
  for (i = 0; i < n; ++i) {
    if (deltax[i] != (float)0.) {
      step[i] = d_sign(step[i], deltax[i]);
    } else {
      step[i] = -step[i];
    }
  }
  return 0;
} 


/* sortd uses the ?shakersort? method to sort an array of keys */
/* in decreasing order. The sort is performed implicitly by */
/* modifying a vector of indices. */

/* For nearly sorted arrays, sortd requires O(n) comparisons. */
/* for completely unsorted arrays, sortd requires O(n**2) */
/* comparisons and will be inefficient unless n is small. */
void yl_subplex_sortd(int n,/*   n      - number of components */
		       double *xkey,/*   xkey   - double precision vector of keys */
		       int *ix)/*   ix     - int vector of indices */
  /*   ix     - indices satisfy xkey(ix(i)) .ge. xkey(ix(i+1)) */
  /*            for i = 1,...,n-1 */
{
  int ixip1, i, ilast, iswap, ifirst, ixi;
  ifirst = 0;
  iswap = 0;
  ilast = n - 2;
  while (ifirst <= ilast) {
    for (i = ifirst; i <= ilast; ++i) {
      ixi = ix[i];
      ixip1 = ix[i + 1];
      if (xkey[ixi] < xkey[ixip1]) {
	ix[i] = ixip1;
	ix[i + 1] = ixi;
	iswap = i;
      }
    }
    ilast = iswap - 1;
    for (i = ilast; i >= ifirst; --i) {
      ixi = ix[i];
      ixip1 = ix[i + 1];
      if (xkey[ixi] < xkey[ixip1]) {
	ix[i] = ixip1;
	ix[i + 1] = ixi;
	iswap = i;
      }
    }
    ifirst = iswap + 1;
  }
}
