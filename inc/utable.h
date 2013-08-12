#ifndef GUARD_UTABLE_H
#define GUARD_UTABLE_H

//#include <math.h>
//#include <stdlib.h>
//#include <tcl.h>
//#include <tk.h>

#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>
#include <vector>
#include <string>
#include <fstream>

double bessel_J0(double);
double expint(double);
double expIntegral(double);
double interp_Newton_2nd(double, 
                         double, double, double,
                         double, double, double);
double interp_Newton_3rd(double, 
                         double, double, double, double,
                         double, double, double, double);

/* data structure for u_n table */

class Utable {
public:
	Utable();
	~Utable(){cleanup();}
	void readUtableFile(const std::string&);
	void readUtableFile(const char *);
	void writeUtableFile(const std::string&, int, int);
	void writeUtableFile(const char *, int, int);
	double getHeightDiffFunc(int, double, double, double, double);
	static double s_integrandWrapper(double, void *);
	double calculateIntegration(int, double);
	double integrand(double);
	void save_t(std::ofstream&, const std::vector<double>&);
	void save_integration(std::ofstream&, const std::vector<double>&, int);
	void cleanup_splines_accs();
	void cleanup();
	double cailleApprox(int, double);
	double interp_utable(int, double);
	double gsl_interp_utable(int, double);
	double getIntegration(int, double);
private:
  int _n;
  double _t;
	double _tmax;
	double _nmax;
	std::vector<double> utab;
	int jump;
  // gsl integration and interpolation objects
  gsl_integration_workspace *workspace;
  std::vector<gsl_interp_accel *> accs;
  std::vector<gsl_spline *> splines;
  
  friend class ModelCalculator;
};

#endif

