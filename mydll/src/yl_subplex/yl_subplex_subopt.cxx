#include "mydll.h"
#include <math.h>
#include <stdio.h>

extern struct {
    double alpha, beta, gamma, delta, psi, omega;
    int nsmin, nsmax, irepl, ifxsw;
    double bonus, fstop;
    int nfstop, nfxe;
    double fxstat[4], ftest;
    bool minf, initx, newx;
} usubc_;

#define usubc_1 usubc_

int yl_subplex_subopt_(int n)
{
printf("*************subopt used*************\n");

/*                                         Coded by Tom Rowan */
/* subopt sets options for subplx. */
/* input */
/*   n      - problem dimension */
/* *********************************************************** */
/* simplex method strategy parameters */
/* *********************************************************** */
/* alpha  - reflection coefficient */
/*          alpha > 0 */

    usubc_1.alpha = 1.;

/* beta   - contraction coefficient */
/*          0 < beta < 1 */

    usubc_1.beta = .5;

/* gamma  - expansion coefficient */
/*          gamma > 1 */

    usubc_1.gamma = 2.;

/* delta  - shrinkage (massive contraction) coefficient */
/*          0 < delta < 1 */

    usubc_1.delta = .5;

/* *********************************************************** */
/* subplex method strategy parameters */
/* *********************************************************** */

/* psi    - simplex reduction coefficient */
/*          0 < psi < 1 */

    usubc_1.psi = .25;

/* omega  - step reduction coefficient */
/*          0 < omega < 1 */

    usubc_1.omega = .1;

/* nsmin and nsmax specify a range of subspace dimensions. */
/* In addition to satisfying  1 <= nsmin <= nsmax <= n, */
/* nsmin and nsmax must be chosen so that n can be expressed */
/* as a sum of positive ints where each of these ints */
/* ns(i) satisfies   nsmin <= ns(i) .ge. nsmax. */
/* Specifically, */
/*     nsmin*ceil(n/nsmax) <= n   must be true. */

/* nsmin  - subspace dimension minimum */

    usubc_1.nsmin = yfmin(2,n);

/* nsmax  - subspace dimension maximum */

    usubc_1.nsmax = yfmin(5,n);

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

    usubc_1.irepl = 0;

/* ifxsw  - measurement replication optimization switch */
/* ifxsw  = 1, 2, or 3 */
/*        = 1 : retain mean of trials as best function value */
/*        = 2 : retain max */
/*        = 3 : retain min */

    usubc_1.ifxsw = 1;

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

    usubc_1.bonus = 1.;

/* nfstop = 0 : f(x) is not tested against fstop */
/*        = 1 : if f(x) has reached fstop, subplx returns */
/*              iflag = 2 */
/*        = 2 : (only valid when irepl > 0) */
/*              if f(x) has reached fstop and */
/*              nfxe > nfstop, subplx returns iflag = 2 */

    usubc_1.nfstop = 0;

/* fstop  - f target value */
/*          Its usage is determined by the value of nfstop. */

/* minf   - bool switch */
/*        = .true.  : subplx performs minimization */
/*        = .false. : subplx performs maximization */

    usubc_1.minf = true;
    usubc_1.minf = false;
    return 0;
} /* subopt_ */

