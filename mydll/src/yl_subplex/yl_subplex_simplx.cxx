#include "mydll.h"
#include "yl_subplex.h"
#include <stdio.h>

//extern struct {
//  double fbonus, sfstop, sfbest;
//  bool neww;
//} isubc_1;
extern isubc isubc_1;

//#define isubc_1 isubc_

/*                                         Coded by Tom Rowan */
/* simplx uses the Nelder-Mead simplex method to minimize the */
/* function f on a subspace. */
/* input */
/*   f      - function to be minimized, declared external in */
/*            calling routine */
/*   n      - problem dimension */
/*   step   - stepsizes for corresponding components of x */
/*   ns     - subspace dimension */
/*   ips    - permutation vector */
/*   maxnfe - maximum number of function evaluations */
/*   cmode  - bool switch */
/*            = .true.  : continuation of previous call */
/*            = .false. : first call */
/*   x      - starting guess for minimum */
/*   fx     - value of f at x */
/*   nfe    - number of function evaluations */
/*   s      - double precision work array of dimension .ge. */
/*            ns*(ns+3) used to store simplex */
/*   fs     - double precision work array of dimension .ge. */
/*            ns+1 used to store function values of simplex */
/*            vertices */

/* output */
/*   x      - computed minimum */
/*   fx     - value of f at x */
/*   nfe    - incremented number of function evaluations */
/*   iflag  - error flag */
/*            = -1 : maxnfe exceeded */
/*            =  0 : simplex reduced by factor of psi */
/*            =  1 : limit of machine precision */
/*            =  2 : reached fstop */

int yl_subplex_simplx(YLMinimizationFunc  f,
		       int n,
		       double *step,
		       int ns,
		       int *ips,
		       int maxnfe,
		       bool *cmode,
		       double *x,
		       double *fx,
		       int *nfe,
		       double *s,
		       double *fs,
		       int *iflag,
		       YLSUBPLEXctrl* ctrl)
{
  static int inew;
  static int npts;
  int  j;
  bool small;
  static int itemp;
  static double fc, fe;
  static int ih, il;
  static double fr;
  static int is;
  static bool updatc;
  static double dum, tol;



  /* Function Body */
  if (*cmode) {
    goto L50;
  }
  npts = ns + 1;
  itemp = (ns + 2)*ns;
  updatc = false;
  yl_subplex_start(n, x, step, ns, ips, s, &small);
  if (small) {
    *iflag = 1;
    //printf("return pos1\n");
    return 0;
  }
  /*debug
    printf("pass check point 1\n");
  */
  if (ctrl->irepl > 0) {
    isubc_1.neww = false;
    yl_subplex_evalf(f, ns, ips, s, n, x, fs, nfe, ctrl);
  } else {
    *fs = *fx;
  }
  isubc_1.neww = true;
  for (j = 1; j < npts; ++j) {
    yl_subplex_evalf(f, ns, ips, &s[j*ns], n, x, fs+j, nfe, ctrl);
  }
  il = 0;//???
  yl_subplex_order(npts, fs, &il, &is, &ih);
  tol = ctrl->psi * yl_blas_dist(ns, s+ih*ns,1, s+il*ns,1);

  /* main loop */

 L20:
  yl_subplex_calcc(ns, s, ih, inew, updatc, s+npts*ns);
  updatc = true;
  inew = ih;

  /* reflect */

  yl_subplex_newpt(ns, ctrl->alpha, s+npts*ns, s+ih*ns, true, s+itemp, &small);
  if (small) {
    goto L40;
  }
  yl_subplex_evalf(f, ns, ips, s+itemp, n, x, &fr, nfe, ctrl);
  if (fr < fs[il]) {

    /* expand */
    yl_subplex_newpt(ns, -ctrl->gamma, s+npts*ns, s+itemp, true, s+ih*ns, &small);
    if (small) {
      goto L40;
    }
    yl_subplex_evalf(f, ns, ips, s+ih*ns, n, x, &fe, nfe, ctrl);
    if (fe < fr) {
      fs[ih] = fe;
    } else {
      yl_blas_dcopy(ns, s+itemp, 1, s+ih*ns, 1);
      fs[ih] = fr;
    }
  } else if (fr < fs[is]) {

    /* accept reflected point */

    yl_blas_dcopy(ns, s+itemp, 1, s+ih*ns, 1);
    fs[ih] = fr;
  } else {

    /* contract */

    if (fr > fs[ih]) {
      yl_subplex_newpt(ns, -ctrl->beta, s+npts*ns, s+ih*ns, true, s+itemp, &small);
    } else {
      yl_subplex_newpt(ns, -ctrl->beta, s+npts*ns, s+itemp, false, &dum, &small);
    }
    if (small) {
      goto L40;
    }
    yl_subplex_evalf(f, ns, ips, s+itemp, n, x, &fc, nfe, ctrl);
    if (fc < yfmin(fr,fs[ih])) {
      yl_blas_dcopy(ns, s+itemp, 1, s+ih*ns, 1);
      fs[ih] = fc;
    } else {

      /* shrink simplex */
      for (j = 0; j < npts; ++j) {
	if (j != il) {
	  yl_subplex_newpt(ns, -ctrl->delta, s+il*ns, s+j*ns, false, &dum, &small);
	  if (small) {
	    goto L40;
	  }
	  yl_subplex_evalf(f, ns, ips, s+j*ns, n, x, fs+j, nfe, ctrl);
	}
      }
    }
    updatc = false;
  }
  yl_subplex_order(npts, fs, &il, &is, &ih);

  /*       check termination */

 L40:
  if (ctrl->irepl == 0)    *fx = fs[il];
  else    *fx = isubc_1.sfbest;
  
 L50:
  if (ctrl->nfstop > 0 && *fx <= isubc_1.sfstop && ctrl->nfxe >= ctrl->nfstop) {
    *iflag = 2;
  } else if (*nfe >= maxnfe) {
    *iflag = -1;
  } else if (yl_blas_dist(ns, s+ih*ns,1, s+il*ns,1) <= tol || small) {
    *iflag = 0;
  } else {
    goto L20;
  }
  
  /*     end main loop, return best point */
  for (j = 0; j < ns; ++j) {
    x[ips[j]] = s[j + il*ns];
  }
  return 0;
} /* simplx_ */

