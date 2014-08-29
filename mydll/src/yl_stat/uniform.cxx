#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <values.h>
#include <stack>
#include <signal.h>
#include <unistd.h>

#include <typeinfo>
#include <math.h>
#include <fstream>
#include <iostream>
#include <string>
#include "mydll.h"

// ***********************************************************************
// 
// uniform - class for uniform distribution
// 
// ***********************************************************************

// **************************************************
// *** constructor

uniform::uniform(int __dimension) :  dimension(__dimension){
  assert(dimension > 0);
  limits = NULL;
  allocate();
}

// **************************************************
// *** memory allocation

void uniform::allocate(){
  assert(dimension > 0);
  realloc(limits, dimension * sizeof(dduple) );
  for (int d = 0; d < dimension; d++){
    limits[d].x = 0.0;
    limits[d].y = 1.0;
  }
}

void uniform::deallocate(){
  free(limits);
}


// **************************************************
// *** destructor
uniform::~uniform(){
  deallocate();
}

// **************************************************
// *** direct element access

double & uniform::lower(int i){
  assert(i >= 0 && i < dimension);
  return limits[i].x;
}

double & uniform::upper(int i){
  assert(i >= 0 && i < dimension);
  return limits[i].y;
}

// **************************************************
// *** generate uniformly distributed samples

void uniform::generateSample(cmatrix &data, int numSamples ){
  assert(data.rows() >= numSamples);
  assert(numSamples > 0);
  assert(data.cols() == dimension);
  // --------------------------------------------------
  // generate samples

  for (int d = 0; d < dimension; d++){ // loop through dimensions
    double llim = limits[d].x;
    double range = (limits[d].y - llim);
    for (int i = 0; i < numSamples; i++){	// loop through samples
      data.value(i,d) = llim + (urand0() * range);
    }
  }
}




float urand0(){
  return rand()/(float)RAND_MAX;
}

float gasdev(){
  static int iset=0;
  static float gset;
  float fac,rsq,v1,v2;

  if  (iset == 0) {
    do {
      v1=2.0*urand0()-1.0;
      v2=2.0*urand0()-1.0;
      rsq=v1*v1+v2*v2;
    } while (rsq >= 1.0 || rsq == 0.0);
    fac=sqrt(-2.0*log(rsq)/rsq);
    gset=v1*fac;
    iset=1;
    return v2*fac;
  } else {
    iset=0;
    return gset;
  }
}


int bisearch(float *xx, int n, float x){
  int ju,jm,jl;
  bool ascnd;

  jl=0;
  ju=n;
  ascnd = (xx[n-1] >= xx[0]);
  while (ju-jl > 1) {
    jm=(ju+jl) >> 1;
    if ( (x >= xx[jm]) == ascnd)
      jl=jm;
    else
      ju=jm;
  }
  return jl;
}
