#include "mydll.h"
#include <math.h>
#include "yl_subplex.h"
/* partx partitions the vector x by grouping components of */
/* similar magnitude of change. */
/*   n      - number of components (problem dimension) */
/*   ip     - permutation vector */
/*   absdx  - vector of magnitude of change in x */
/*   nsvals - int array dimensioned .ge. int(n/nsmin) */

/* output */
/*   nsubs  - number of subspaces */
/*   nsvals - int array of subspace dimensions */

void yl_subplex_partx(int n,
		       int *ip,
		       double *absdx,
		       int *nsubs,
		       int *nsvals,
		       YLSUBPLEXctrl* ctrl)
{
  int i1;
  int i, nleft, nused;
  double as1max=0, gapmax, asleft, as1, as2;
  int ns1, ns2;
  double gap;

  *nsubs = 0;
  nused = 0;
  nleft = n;

  asleft = *absdx;
  for (i = 1; i < n; ++i) 	asleft += absdx[i];

  while(nused < n) {

    i1 = ctrl->nsmin-1;
    for (as1 = 0., i = 0; i < i1; ++i)
      as1 += absdx[ip[nused + i]];

    gapmax = -1.;
    for(ns1 = i1,
	  i1 = yfmin(ctrl->nsmax, nleft); ns1 < i1; ++ns1)
      {
	as1 += absdx[ip[nused + ns1]];
	ns2 = nleft - ns1 - 1;

	if (ns2 > 0) {
	  if (ns2 >= ((ns2 - 1) / ctrl->nsmax + 1) * ctrl->nsmin) {
	    as2 = asleft - as1;
	    gap = as1 / (ns1+1) - as2 / ns2;
	    if (gap > gapmax) {
	      gapmax = gap;
	      nsvals[*nsubs] = (ns1+1);
	      as1max = as1;
	    }
	  }
	} else {
	  if (as1 / (ns1+1) > gapmax) {
	    nsvals[*nsubs] = (ns1+1);
	    ++(*nsubs);
	    return;
	  }
	}
      }
    nused += nsvals[*nsubs];
    nleft = n - nused;
    asleft -= as1max;
    ++(*nsubs);
  }
} 

