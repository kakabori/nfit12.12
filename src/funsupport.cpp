#include <math.h>
#include <stdio.h>
#include <queue>
#include <vector>
#include <algorithm>
#include <stdexcept>
#include <iostream>
#include "funsupport.h"

using std::vector; using std::cout; using std::endl;
using std::domain_error;


/* This function is a predicate used to sort xyvec */
bool compare(const Two_Doubles& x, const Two_Doubles& y)
{
  return x.x < y.x;
}


// initialize pointers to indicate that memory has not been allocated
FunSupport::FunSupport(){
  acc = NULL;
  spline = NULL;
  setSplineTypeToCubic();
}


/* test accuracy */
bool FunSupport::on_same_line( double dx, double dy, double y2){
  if( dx > maxdx) return false;
  return ( dx < mindx || dy < abserr || dy < fabs(y2*relerr) );
}


// A helper function that iteratively finds support points
/*This function had to be re-written to support multithreaded calculations,
  because the vast number of stack frames it generated would overflow the
  stack when this function was run in multiple threads concurrently. This
  approach transfers the memory burden to the heap, which on a modern machine
  with a large amount of ram can be quite large.
*/
void FunSupport::_findPoints(double initLeft, double initRight, double initVleft, double initVright)
{
  using namespace std;

  double left, right, vleft, vright;
  double middle, vmiddle;
  queue<double> recursionQueue;
  Two_Doubles tmp;

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
    vmiddle = func(middle, para);
    
    //cout << middle << " " << vmiddle << endl;

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
    // If the function is linear within tolerance over the given range, add the
    // middle and right points to the list of interpolation points.
    } else {
      tmp.x = middle;
      tmp.y = vmiddle;
      xyvec.push_back(tmp);
      tmp.x = right;
      tmp.y = vright;
      xyvec.push_back(tmp);
    }
  }

  // The points are not guaranteed to be in sorted order, so x values must be
  // sorted. 
  sort(xyvec.begin(), xyvec.end(), compare);
}


/* find the support points given the range (left, right) */
void FunSupport::findPoints(double left, double right) {
  Two_Doubles tmp;
  double vleft = func(left, para);
  double vright = func(right, para);
  xyvec.resize(0);
  tmp.x = left;
  tmp.y = vleft;
  xyvec.push_back(tmp);

  _findPoints(left, right, vleft, vright);
  
  // Copy content of xyvec to x and y vectors so that they can be used in
  // gsl_spline object
  vector<double> x, y;
  vector<Two_Doubles>::iterator it;
  for (it = xyvec.begin(); it != xyvec.end(); it++) {
    x.push_back(it->x);
    y.push_back(it->y);
  }

  if(acc!=NULL) gsl_interp_accel_free(acc);
  acc = gsl_interp_accel_alloc();
  if(spline!=NULL) gsl_spline_free (spline);
  spline = (type == enum_cubic) ? gsl_spline_alloc(gsl_interp_cspline, x.size()):
                                  gsl_spline_alloc(gsl_interp_linear, x.size());
  gsl_spline_init(spline, &x[0], &y[0], x.size());
  xmin = x[0];
  xmax = x.back();
}


/* return the interpolated (splined) value */
double FunSupport::val(double xv) 
{ 
  if (xv < xmin || xv > xmax) {
    cout << "xv: " << xv << " xmin: " << xmin << " xmax: " << xmax << endl;
    throw domain_error("The requested x value is outside of the range.");
  } 
  
  return gsl_spline_eval(spline, xv, acc);
}


void FunSupport::getPoints(vector<double>& x, vector<double>& y)
{
	vector<Two_Doubles>::iterator it;
	for (it = xyvec.begin(); it != xyvec.end(); it++) {
  	x.push_back(it->x);
  	y.push_back(it->y);
  }
}
