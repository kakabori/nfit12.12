#include "mydll.h"
#include "yl_subplex.h"
#include <math.h>
#include <stdio.h>

//extern struct {
//    double fbonus, sfstop, sfbest;
//    bool neww;
//} isubc_1;
extern isubc isubc_1;

//#define isubc_1 isubc_



/*                                         Coded by Tom Rowan */
/* subplx uses the subplex method to solve unconstrained */
/* optimization problems.  The method is well suited for */
/* optimizing objective functions that are noisy or are */
/* discontinuous at the solution. */

/* subplx sets default optimization options by calling the */
/* subroutine subopt.  The user can override these defaults */
/* by calling subopt prior to calling subplx, changing the */
/* appropriate common variables, and setting the value of */
/* mode as indicated below. */

/* By default, subplx performs minimization. */

/* input */

/*   f      - user supplied function f(n,x) to be optimized, */
/*            declared external in calling routine */

/*   n      - problem dimension */

/*   tol    - relative error tolerance for x (tol .ge. 0.) */

/*   maxnfe - maximum number of function evaluations */

/*   mode   - int mode switch with binary expansion */
/*            (bit 1) (bit 0) : */
/*            bit 0 = 0 : first call to subplx */
/*                  = 1 : continuation of previous call */
/*            bit 1 = 0 : use default options */
/*                  = 1 : user set options */

/*   scale  - scale and initial stepsizes for corresponding */
/*            components of x */
/*            (If scale(1) .lt. 0., */
/*            abs(scale(1)) is used for all components of x, */
/*            and scale(2),...,scale(n) are not referenced.) */

/*   x      - starting guess for optimum */

/*   work   - double precision work array of dimension .ge. */
/*            2*n + nsmax*(nsmax+4) + 1 */
/*            (nsmax is set in subroutine subopt. */
/*            default: nsmax = min(5,n)) */

/*   iwork  - int work array of dimension .ge. */
/*            n + int(n/nsmin) */
/*            (nsmin is set in subroutine subopt. */
/*            default: nsmin = min(2,n)) */

/* output */

/*   x      - computed optimum */
/*   fx     - value of f at x */
/*   nfe    - number of function evaluations */
/*   iflag  - error flag */
/*            = -2 : invalid input */
/*            = -1 : maxnfe exceeded */
/*            =  0 : tol satisfied */
/*            =  1 : limit of machine precision */
/*            =  2 : fstop reached (fstop usage is determined */
/*                   by values of options minf, nfstop, and */
/*                   irepl. default: f(x) not tested against */
/*                   fstop) */
/*            iflag should not be reset between calls to */
/*            subplx. */


int yl_subplex_subplx(YLMinimizationFunc f,
		       int n,
		       double tol,
		       int maxnfe,
		       int mode,
		       double *scale,
		       double *x,
		       double *fx,
		       int *nfe,
		       double *work,
		       int *iwork,
		       int *iflag,
		       void* _ctrl)
{
  /* Initialized data */
  YLSUBPLEXctrl* ctrl = (YLSUBPLEXctrl*) _ctrl;
  static const double bnsfac[6]	/* was [3][2] */ = { -1.,-2.,0.,1.,0.,2. };

  /* Local variables */
  int i;
  static bool cmode;
  static int nsubs, ipptr;
  static int ns, insfnl, ifsptr;
  static double scl;
  static int ins;
  static double sfx;

  
  if (mode % 2 == 0) {
    /* first call, check input */
    if (n < 1 || tol < 0. || maxnfe < 1) return *iflag = iflagINPUT;
    if (scale[0] > 0.) {
      for (i = 0; i < n; ++i)
	if (scale[i] == 0) return *iflag = iflagINPUT;
    } else {
      scl = fabs(scale[0]);
      if (scl == 0) return *iflag = iflagINPUT;
    }
    if (ctrl->alpha <= 0. ||
	ctrl->beta <= 0.  ||
	ctrl->beta >= 1.  ||
	ctrl->gamma <= 1. ||
	ctrl->delta <= 0. ||
	ctrl->delta >= 1. ||
	ctrl->psi <= 0.   ||
	ctrl->psi >= 1.   ||
	ctrl->omega <= 0. ||
	ctrl->omega >= 1. ||
	ctrl->nsmin < 1   ||
	ctrl->nsmax < ctrl->nsmin ||
	n < ctrl->nsmax   ||
	n < ((n - 1) / ctrl->nsmax + 1) * ctrl->nsmin ||
	ctrl->irepl < 0   ||
	ctrl->irepl > 2   ||
	ctrl->ifxsw < 1   ||
	ctrl->ifxsw > 3   ||
	ctrl->bonus < 0.  ||
	ctrl->nfstop < 0) {
      return *iflag = iflagINPUT;
    }

    /*       initialization */

    ifsptr = 2*n + ctrl->nsmax * (ctrl->nsmax + 3);
    if (*scale > 0.) {
      yl_blas_dcopy(n, scale, 1, work, 1);
      yl_blas_dcopy(n, scale, 1, work+n, 1);
      //printf("subplx: branch 1\n");
    } else {
      yl_blas_dcopy(n, &scl, 0, work, 1);
      yl_blas_dcopy(n, &scl, 0, work+n, 1);
      //printf("subplx: branch 2 n=%d\n",n);
      //printf("scl = %g work=%g\n",scl, work[n]);
    }
    for (i = 0; i < n; ++i) iwork[i] = i;
    *nfe = 0;
    ctrl->nfxe = 1;
    if (ctrl->irepl == 0) {
      isubc_1.fbonus = 0.;
    } else {
      isubc_1.fbonus = bnsfac[ctrl->ifxsw - 1] * ctrl->bonus;
    }
    if (ctrl->nfstop == 0) {
      isubc_1.sfstop = 0.;
    } else{
      isubc_1.sfstop = ctrl->fstop;
    }
    ctrl->ftest = 0.;
    cmode = false;
    isubc_1.neww = true;
    ctrl->initx = true;
    yl_subplex_evalf(f, 0, iwork, NULL, n, x, &sfx, nfe, ctrl);
    ctrl->initx = false;

  } else {
    /*       continuation of previous call */
    switch(*iflag){

    case iflagFSTOP:
      isubc_1.sfstop = ctrl->fstop;
      cmode = true;
      goto L70;

    case iflagMAXNFE:
      cmode = true;
      goto L70;
		
    case iflagTOL:
      cmode = false;
      goto CHKterm;
		
    default:
      return *iflag;
    }
  }

  /*     subplex loop */
  // do {
 L40:
  for (i = 0; i < n; ++i) work[i] = fabs(work[i]);
  yl_subplex_sortd(n, work, iwork);
  yl_subplex_partx(n, iwork, work, &nsubs, iwork+n, ctrl);

  //printf("nsubs=%d\n", nsubs);
  //for(i=0;i<nsubs;i++) printf("dim=%d\n", *(iwork+n+i));
    
  yl_blas_dcopy(n, x, 1, work, 1);
  ins = n;
  insfnl = n + nsubs - 1;
  ipptr = 0;
    
  /* simplex loop */

 L60:
  ns = iwork[ins];
 L70:
  yl_subplex_simplx(f, n, work+n, ns, &iwork[ipptr], maxnfe, &cmode, x,
		     &sfx, nfe, work+2*n, &work[ifsptr], iflag, ctrl);

  cmode = false;
  if (*iflag != 0) {
    //      *fx = ctrl->minf ? sfx : -sfx;
    *fx = sfx;
    return *iflag;
  }
  if (ins < insfnl) {
    ++ins;
    ipptr += ns;
    goto L60;
  }
    
  /*  end simplex loop */
  for (i = 0; i < n; ++i) work[i] = x[i] - work[i];

  /* check termination */
 CHKterm:
  for (i = 0; i < n; ++i) {
    if ( yfmax(fabs(work[i]), fabs(work[i+n]) * ctrl->psi) / yfmax( fabs(x[i]), 1.) > tol) {
      yl_subplex_setstp(&nsubs, n, work, work+n, ctrl);
      goto L40;
    }
  }
  /* end subplex loop */

  *iflag = iflagTOL;
  //    *fx = ctrl->minf ? sfx : -sfx;
  *fx = sfx;
  return *iflag;
}

// work[0..n] scale || x
// iwork[0..n] permutation
// iwork[n..n+nsubs] nsvals
