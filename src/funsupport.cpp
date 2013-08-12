#include <math.h>
#include <stdio.h>
#include <queue>
#include <vector>
#include <algorithm>
#include <iostream>

#include "funsupport.h"

using std::vector;

FunSupport::FunSupport(){
  /*
    initialize pointers to indicate that memory has
    not been allocated
  */
  acc = NULL;
  spline = NULL;
  setSplineTypeToCubic();
}

bool FunSupport::on_same_line( double dx, double dy, double y2){
  /* test accuracy */
  if( dx > maxdx) return false;
  return ( dx < mindx || dy < abserr || dy < fabs(y2*relerr) );
}

void FunSupport::_findPoints(double initLeft, double initRight, double initVleft, double initVright)
{
  /*
    A helper function that iteratively finds support points
  */

  using namespace std;

  double left, right, vleft, vright;
  double middle, vmiddle;
  queue<double> recursionQueue;
  vector<double>::iterator iter;

/*This function had to be re-written to support multithreaded calculations,
  because the vast number of stack frames it generated would overflow the
  stack when this function was run in multiple threads concurrently. This
  approach transfers the memory burden to the heap, which on a modern machine
  with a large amount of ram can be quite large.
*/


/* Start by setting up a queue with the initial right and left end points and
   their corresponding function values.
*/
  recursionQueue.push(initLeft);
  recursionQueue.push(initRight);
  recursionQueue.push(initVleft);
  recursionQueue.push(initVright);

  while(! recursionQueue.empty()) {
    /*If the function is not linear between two end points within
      specified tolerance, split the range in half, and add two sets of
      points and function values to the back of the queue.
    */
    left = recursionQueue.front();
    recursionQueue.pop();
    right = recursionQueue.front();
    recursionQueue.pop();
    vleft = recursionQueue.front();
    recursionQueue.pop();
    vright = recursionQueue.front();
    recursionQueue.pop();
    middle = (left + right) / 2;
    vmiddle = func( middle, para );

    //cout << left << " " << middle << " " << right << " "  << vleft << " "  << vmiddle
    //<< " "  << vright << endl;


    if( !on_same_line(middle-left, fabs(vleft+vright-vmiddle-vmiddle) / 2 , vmiddle)) {
      recursionQueue.push(left);
      recursionQueue.push(middle);
      recursionQueue.push(vleft);
      recursionQueue.push(vmiddle);

      recursionQueue.push(middle);
      recursionQueue.push(right);
      recursionQueue.push(vmiddle);
      recursionQueue.push(vright);
       /*If the function is linear within tolerance over the given range, add the
         middle and right points to the list of interpolation points.*/
    } else {
      x.push_back(middle);
      x.push_back(right);
    }
  }

    /*The points are not guaranteed to be in sorted order, so x values must be
    sorted. Finally corresponding function values are added to the vector of y
    values.*/
  sort(x.begin(), x.end());
  for(iter = x.begin(); iter != x.end(); iter ++) {
    y.push_back(func(*iter, para));
  }
}

/* find the support points given the range (left, right) */
void FunSupport::findPoints(double left, double right)
{
  x.resize(0);
  y.resize(0);
  double vleft = func(left, para);
  double vright = func(right, para);
  x.push_back(left);

  _findPoints(left, right, vleft, vright);

  if(acc!=NULL) gsl_interp_accel_free(acc);
  acc = gsl_interp_accel_alloc();
  if(spline!=NULL) gsl_spline_free (spline);
  spline = (type == enum_cubic) ? gsl_spline_alloc(gsl_interp_cspline, x.size()):
                                  gsl_spline_alloc(gsl_interp_linear, x.size());
  gsl_spline_init(spline, &x[0], &y[0], x.size());
}

void FunSupport::print(FILE *fp){
  /* for debugging */
  for(unsigned int i=0; i<x.size(); i++){
    fprintf(fp, "%g %g\n", x[i], y[i]);
  }
  fprintf(fp, "\n");
  fflush(fp);
}

double FunSupport::val(double xv){
  /* return the interpolated (splined) value */
  return gsl_spline_eval(spline, xv, acc);
}

void FunSupport::getPoints(vector<double>& _x, vector<double>& _y)
{
	vector<double>::iterator iter;
	for (iter = x.begin(); iter != x.end(); iter++) {
  	_x.push_back(*iter);
  }
  for (iter = y.begin(); iter != y.end(); iter++) {
  	_y.push_back(*iter);
  }
}
