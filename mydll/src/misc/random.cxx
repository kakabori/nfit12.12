#include "mydll.h"
#include <cmath>
#include <cstdlib>
#include <cstdio>

void locate(float *xx, int n, float x, int *j)
{
  int ju,jm,jl;
  bool ascnd;

  jl=0;
  ju=n;
  ascnd = (xx[n-1] >= xx[0]);
  while (ju-jl > 1) {
    jm=(ju+jl) >> 1;
    if (x >= xx[jm] == ascnd)
      jl=jm;
    else
      ju=jm;
  }
  *j=jl;
}

