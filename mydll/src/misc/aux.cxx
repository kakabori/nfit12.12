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


using namespace std;

// check if number in interval
bool boundCheck(double value, double min, double max){
  assert(min <= max);
  return (value >= min && value <= max);
}


// Euclieadn distance
double L2Dist(double x1, double y1, double x2, double y2){
  return sqrt(sqr(x1 - x2) + sqr(y1 - y2));
}


// Truncation and Cyclic Truncation
double cyclicTruncate(double x, double min, double max){
  double delta = max - min;
  assert (delta > 0.0);
  x = fmod(x, delta);
  for (; x <  min; ) x += delta;
  for (; x >= max; ) x -= delta;
  return x;
}


double regularTruncate(double x, double min, double max){
  assert (max >= min);
  if (x < min) x = min;
  else if (x > max) x = max;
  return x;
}
