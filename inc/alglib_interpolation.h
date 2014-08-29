#ifndef GUARD_ALGLIB_INTERPOLATION_H
#define GUARD_ALGLIB_INTERPOLATION_H

#include "~/programs/alglib/src/stdafx.h"
#include "~/programs/alglib/src/interpolation.h"
#include <vector>

class Alglib_CubicSpline2D {
public:
  void initialize(double (*f)(double, double, void*), void *p) {Function = f; parameter = p;}
  void set_DistanceError_Bounds(double a, double b, double c, double d, double e, double f) {
    min_dx = a, max_dx = b, min_dy = c, max_dy = d, abserr = e; relerr = f;}
  bool onStraightLine(double, double, double, double);
  int buildGridPoints(double, double, double, std::vector<double>&);
  void buildInterpolant(double, double, double, double);
  double evaluate(double, double);
private:
  double max_dx, max_dy; //maximum seperation between neighbouring points allowed by the user
  double min_dx, min_dy; //minimum seperation between neighbouring points allowed by the user
  double abserr; //absolute error tolerance
  double relerr; //relative error tolerance
  void *parameter;
  double (*Function)(double, double, void*); //the function to be interpolated takes two doubles for the independent variables and one void pointer that points to the parameters to evaluate the function
  alglib::real_1d_array xarray, yarray, farray;
  spline2dinterpolant spline;
};

class Alglib_LinearSpline3D {
public:

private:
};

#endif
