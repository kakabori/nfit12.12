#ifndef GUARD_FUNSUPPORT_H
#define GUARD_FUNSUPPORT_H

#include <vector>
#include <gsl/gsl_spline.h>


/*
  given a function f(x, other fixed parameters), we want to approximate
  the function with a set of discrete values (f_i, x_i).
  When f(x) is needed, an interpolated value will be returned.

  usage:
  1. init() : give the class the function pointer and extra parameters
  2. setol(): set the tolerance of the approximation
  3. findPoints(): find a set of discrete values to calculate f
  4. val()  : get the interpolated value given an x
*/

struct Two_Doubles {
  double x, y;
};
 

class FunSupport {
public:
	enum SplineType_Enum {
		enum_linear,
		enum_cubic
	};
  FunSupport();
  void setSplineTypeToLinear() {type = enum_linear;}
  void setSplineTypeToCubic() {type = enum_cubic;}
  bool on_same_line(double, double, double);
  void _findPoints( double, double, double, double);
  void findPoints(double a, double b);
  void print(FILE *fp);
  void init(double (*f)(double, void*), void *p){func=f;para=p;}
  void setol(double a, double r, double m1, double m2){abserr=a; relerr=r; mindx=m1; maxdx=m2;}
  double val(double);
	void getPoints(std::vector<double>&, std::vector<double>&);
	double get_lower_bound() {return xmin;}
	double get_upper_bound() {return xmax;}
private:
  double maxdx; // maximum seperation between neighbouring points
  double mindx; // minimum seperation between neighbouring points
  double abserr; // absolute error tolerance
  double relerr; // relative error tolerance
  std::vector<Two_Doubles> xyvec;
  void *para;
  gsl_interp_accel *acc;
  gsl_spline *spline;
  double (*func)(double, void*);
  SplineType_Enum type;
  double xmin, xmax; // Define the input range to the interpolating function
};

#endif
