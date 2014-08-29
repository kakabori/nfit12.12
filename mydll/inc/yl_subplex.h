#ifndef GUARD_YL_SUBPLEX_H
#define GUARD_YL_SUBPLEX_H

#include "mydll.h"

/* *********************************************************** */
/* simplex method strategy parameters */
/* alpha, beta, gamma, delta          */
/* *********************************************************** */
/* subplex method strategy parameters */
/* psi, omega, nsmin, msmax*/
/* *********************************************************** */
/* nsmin and nsmax specify a range of subspace dimensions. */
/* In addition to satisfying  1 <= nsmin <= nsmax <= n, */
/* nsmin and nsmax must be chosen so that n can be expressed */
/* as a sum of positive ints where each of these ints */
/* ns(i) satisfies   nsmin <= ns(i) .ge. nsmax. */
/* Specifically, */
/*     nsmin*ceil(n/nsmax) <= n   must be true. */
/* nsmin  - subspace dimension minimum */
/* nsmax  - subspace dimension maximum */
/************************************************************** */
/* *********************************************************** */
/* subplex method special cases */
/* *********************************************************** */
/* nelder-mead simplex method with periodic restarts */
/*   nsmin = nsmax = n */
/* *********************************************************** */
/* nelder-mead simplex method */
/*   nsmin = nsmax = n, psi = small positive */
/* *********************************************************** */

/* irepl, ifxsw, and bonus deal with measurement replication. */
/* Objective functions subject to large amounts of noise can */
/* cause an optimization method to halt at a false optimum. */
/* An expensive solution to this problem is to evaluate f */
/* several times at each point and return the average (or max */
/* or min) of these trials as the function value.  subplx */
/* performs measurement replication only at the current best */
/* point. The longer a point is retained as best, the more */
/* accurate its function value becomes. */

/* The common variable nfxe contains the number of function */
/* evaluations at the current best point. fxstat contains the */
/* mean, max, min, and standard deviation of these trials. */

/* irepl  - measurement replication switch */
/* irepl  = 0, 1, or 2 */
/*        = 0 : no measurement replication */
/*        = 1 : subplx performs measurement replication */
/*        = 2 : user performs measurement replication */
/*              (This is useful when optimizing on the mean, */
/*              max, or min of trials is insufficient. Common */
/*              variable initx is true for first function */
/*              evaluation. newx is true for first trial at */
/*              this point. The user uses subroutine fstats */
/*              within his objective function to maintain */
/*              fxstat. By monitoring newx, the user can tell */
/*              whether to return the function evaluation */
/*              (newx = .true.) or to use the new function */
/*              evaluation to refine the function evaluation */
/*              of the current best point (newx = .false.). */
/*              The common variable ftest gives the function */
/*              value that a new point must beat to be */
/*              considered the new best point.) */

/* ifxsw  - measurement replication optimization switch */
/* ifxsw  = 1, 2, or 3 */
/*        = 1 : retain mean of trials as best function value */
/*        = 2 : retain max */
/*        = 3 : retain min */

/* Since the current best point will also be the most */
/* accurately evaluated point whenever irepl > 0, a bonus */
/* should be added to the function value of the best point */
/* so that the best point is not replaced by a new point */
/* that only appears better because of noise. */
/* subplx uses bonus to determine how many multiples of */
/* fxstat(4) should be added as a bonus to the function */
/* evaluation. (The bonus is adjusted automatically by */
/* subplx when ifxsw or minf is changed.) */

/* bonus  - measurement replication bonus coefficient */
/*          bonus .ge. 0 (normally, bonus = 0 or 1) */
/*        = 0 : bonus not used */
/*        = 1 : bonus used */

/* nfstop = 0 : f(x) is not tested against fstop */
/*        = 1 : if f(x) has reached fstop, subplx returns */
/*              iflag = 2 */
/*        = 2 : (only valid when irepl > 0) */
/*              if f(x) has reached fstop and */
/*              nfxe > nfstop, subplx returns iflag = 2 */

/* fstop  - f target value */
/*          Its usage is determined by the value of nfstop. */
/*********************************************************** */

void yl_subplex_start(int n,
		   double *x,
		   double *step,
		   int ns,
		   int *ips,
		   double *s,
		   bool *small);
int yl_subplex_evalf(YLMinimizationFunc f,
		      int ns,
		      int *ips,
		      double *xs,
		      int n,
		      double *x,
		      double *sfx,
		      int *nfe,
		      YLSUBPLEXctrl* ctrl);
int yl_subplex_fstats(double fx,
		       int ifxwt,
		       bool reset,
		       YLSUBPLEXctrl* ctrl);
int yl_subplex_order(int npts,
		   double *fs,
		   int *il,
		   int *is,
		   int *ih);
int yl_subplex_calcc(int ns,
		      double *s,
		      int ih,
		      int inew,
		      bool updatc,
		      double *c);
int yl_subplex_newpt(int ns,
		   double coef,
		   double *xbase,
		   double *xold,
		   bool neww,
		   double *xneww,
		   bool *small);
int yl_subplex_simplx(YLMinimizationFunc f,
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
		       YLSUBPLEXctrl* ctrl);
void yl_subplex_sortd(int n,
		   double *xkey,
		   int *ix);
void yl_subplex_partx(int n,
		       int *ip,
		       double *absdx,
		       int *nsubs,
		       int *nsvals,
		       YLSUBPLEXctrl* ctrl);
int yl_subplex_setstp(int *nsubs,
		       int n,
		       double *deltax,
		       double *step,
		       YLSUBPLEXctrl* ctrl);

enum iflagConditions
{
	iflagINPUT  = -2,
	iflagMAXNFE = -1,
	iflagTOL    = 0,
	iflagEPS	= 1,
	iflagFSTOP  = 2
};

struct isubc {
  double fbonus, sfstop, sfbest;
  bool neww;
};

#endif
