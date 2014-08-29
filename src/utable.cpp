/*  This calculates the correlation functions in Eq. 3.31
and 3.32 in the thesis */
#include <iostream>
#include <sstream>
#include <iomanip>
#include <math.h>
#include <assert.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_expint.h>
#include "utable.h"

using std::string; using std::ifstream; using std::getline;
using std::istringstream; using std::vector; using std::ofstream;
using std::cout; using std::endl; using std::cerr;
using std::setw; using std::showpoint; using std::noshowpoint;

#define WORKSPACE_SIZE 100000
#define KEY 6
#define SMALLNUM 0.000001
#define PI 3.1415926535897932384626433832795

Utable::Utable()
{
  //Turn off the error handler in gsl_integration_qag
  gsl_set_error_handler_off();
	workspace = gsl_integration_workspace_alloc(WORKSPACE_SIZE);
}

/* Called by the destructor */
void Utable::cleanup()
{
	gsl_integration_workspace_free(workspace);
	cleanup_splines_accs();
}

/* Cleans up the two vectors, splines and accs */
void Utable::cleanup_splines_accs()
{
	typedef vector<gsl_spline *>::size_type sz;
	for (sz i = 0; i < splines.size(); i++) {
		if(splines[i] != NULL) gsl_spline_free (splines[i]);
	}
	splines.clear();
	
	typedef vector<gsl_interp_accel *>::size_type sz2;
	for (sz2 i = 0; i < accs.size(); i++) {
		if(accs[i] != NULL) gsl_interp_accel_free(accs[i]);
	}
	accs.clear();
}

/******************************************************************************
This function calculates the integration for the height difference function.

t represents r/xi, the dimensionless length, that appears in the theory
n is an integer, the sum index

The integration is done from 1e-10 to infinity, not from 0 to inf.
This is to avoid the numerically undefined t = 0 point, where
division by 0 occurs. The integrand itself is finite at t = 0,
so in the future, one can probably calculate the integrand
exactly by doing some calculus, which will allow the integration
to extend down to t = 0.
******************************************************************************/
double Utable::calculateIntegration(int n, double t)
{
	double result, abserr, epsabs = 0, epsrel = 1e-10;
	_n = n;
	_t = t;
	gsl_function F;
	F.function = &s_integrandWrapper;
	F.params = this;
  double lowerLimit = 1e-10;
  double upperLimit = 1e10;
	//gsl_integration_qagiu(&F, 1e-10, epsabs, epsrel, WORKSPACE_SIZE,
	//                      workspace, &result, &abserr);
	gsl_integration_qag(&F, lowerLimit, upperLimit, epsabs, epsrel,
                      WORKSPACE_SIZE, KEY, workspace, &result, &abserr);
  return result;
}

double Utable::s_integrandWrapper(double x, void *param)
{
	Utable *p = (Utable *)param;
  return p->integrand(x);
}

double Utable::integrand(double x)
{
	double ret;
	ret = bessel_J0(sqrt(2 * x) * _t);
	ret = ret * pow(sqrt(1 + x * x) - x, 2 * _n);
	ret = 1 - ret;
	ret = ret / x / sqrt(1 + x * x);
	return ret;
}

/******************************************************************************
writeUtableFile calculates the integration defined inside 
the height difference function and output as a text file.

It evaluates the integration for 0 <= t <= 1.01 at constant 
t step size of 0.01. For t >= 1, the integration is evaluated 
at constant log_10(t) step size defined by log10tstep. The
minimum log_10(t) is always equal to 0. The maximum log_10(t) 
is equal to log10t_numSteps * log10tstep.

There is an overlap range between t = 1 and 1.02. This is for
the conveience of interpolation. The table will be interpolated
with second order polynomials. For 0.99 < t < 1, the function 
values at t = 0.99, 1.00, and 1.01 will be used. For t >= 1, 
the function valued evaluated at log_10(t) will be used.

Note that t = r / xi, which is a dimensionless length.
******************************************************************************/
void Utable::writeUtableFile(const string& filename,
                             int log10t_numSteps = 600, int nmax = 1000)
{
	writeUtableFile(filename.c_str(), log10t_numSteps, nmax);             
}

void Utable::writeUtableFile(const char *filename,
                             int log10t_numSteps = 600, int nmax = 1000)
{ 
	const double log10tstep = 0.01;
	vector<double> x; // vector to store values of t
	double log10t; // temporary variable

	// build x vector from t = 0 to 1.01
	double t = 0; 
	for (int i = 0; i < 102; i++) {
		x.push_back(t);
		t += 0.01;
	}
	
	// Build x vector from log10(t) = 0 to log10tmax,
	// where log10tmax = log10t_numSteps * log10tstep
	// First, always evaluate at t = 1 
	t = 1;
	x.push_back(t);
	for (int i = 0; i < log10t_numSteps; i++) {
		log10t = log10(t) + log10tstep;
		t = pow(10, log10t);
		x.push_back(t);
	}
	
	ofstream myfile;
	myfile.open(filename);
	
	// save t values as the first line in the file
	save_t(myfile, x);
	// save calculated integration as a matrix in the file
	save_integration(myfile, x, nmax);
	
	myfile.close();
}

void Utable::save_t(ofstream& myfile, const vector<double>& x)
{
	typedef vector<double>::size_type sz;
	for (sz i = 0; i < x.size()-1; i++) myfile << setw(14) << x[i];
	// save the last element in x without a trailing white space
	myfile << setw(14) << x[x.size()-1] << endl;
}

void Utable::save_integration(ofstream& myfile, 
                              const vector<double>& x, int nmax)
{
	typedef vector<double>::size_type sz;
	for (int i = 0; i < nmax; i++) {
		for (sz j = 0; j < x.size()-1; j++) {
			myfile << showpoint << setw(14) << getIntegration(i, x[j]);
		}
		// do the last one separately to avoid a trailing white space
		myfile << setw(14) << getIntegration(i, x[x.size()-1]) << endl;
	}
	
	// do the last line separately to avoid an empty line at the end
	for (sz i = 0; i < x.size()-1; i++) {
		myfile << setw(14) << getIntegration(nmax, x[i]);
	}
	myfile << setw(14) << getIntegration(nmax, x[x.size()-1]) 
	       << noshowpoint;
}

/* Decide whether to do the integration or use the Caille approx. */
double Utable::getIntegration(int n, double t)
{
	if (n > 30) {
		if (t > 10) {
			return cailleApprox(n, t);
		} else {
			return calculateIntegration(n, t);
		}
	} else {
		if (t > 1000) {
			return cailleApprox(n, t);
		} else {
			return calculateIntegration(n, t);
		}
	}
}

/* approximation in Ning Lei's thesis, Eq. 2.4.87 */
double Utable::cailleApprox(int n, double t)
{
	//Euler's constant = 0.5772156649
	double t2 = t * t;
	double ret = 2 * 0.5772156649 + log(t2);
	return n == 0 ? ret : ret + expIntegral(t2 / 4 / n);
}

/******************************************************************************
readUtableFile reads an ASCII formatted file for the table of
the integration defined inside the height difference function. 
It expects as the first line,
the values of t = r / xi at which the integration was evaluated.
From the second line, the file should contain the integration values
for n = 0, n = 1, ..., n = nmax, line by line.

After finished reading the file, this function will build a vector
whose elements are gsl_spline interpolants for each n as well as
a vector whose elements are gsl_accel objects.

Note that the interpolation will be done as a function of log_10(t).
******************************************************************************/
void Utable::readUtableFile(const string& filename)
{
	readUtableFile(filename.c_str());
}

void Utable::readUtableFile(const char *filename) 
{
	// clean up the accel and spline objects stored in the vectors
	cleanup_splines_accs();
	utab.clear();
	
	// read in the utable file
	ifstream myfile;
	myfile.open(filename);
	if (!myfile) {
	  cerr << filename << ": file not found. The program exits." << endl;
		exit(1);
	}

	// read in the first line and build log_10(t) vector for spline
	vector<double> x;
	string line;
	double d;
	if (getline(myfile, line)) {
		istringstream iss(line);
		while (iss >> d) x.push_back(log10(d + SMALLNUM));
	}
	// the second to the last element in the first line is the 
	// largest t value for the table
	_tmax = pow(10, x[x.size()-2]) - SMALLNUM;
	cout << "_tmax: " << _tmax << endl;
	assert(x.size() == (102+601));
	jump = x.size();
	x.erase(x.begin()+100);
	x.erase(x.begin()+101);
	
	// build splines and accs vectors
	// read line by line, build an interpolant, and then store
	vector<double> y;
	gsl_interp_accel *acc;
  gsl_spline *spline;
	while (getline(myfile, line)) {
		istringstream iss(line);
		while (iss >> d) {
			y.push_back(d);
			utab.push_back(d);
		}
		y.erase(y.begin()+100);
		y.erase(y.begin()+101);
		assert(x.size() == y.size());
  	spline = gsl_spline_alloc(gsl_interp_cspline, x.size());
		gsl_spline_init(spline, &x[0], &y[0], x.size());
		splines.push_back(spline);
		
		acc = gsl_interp_accel_alloc();
  	accs.push_back(acc);
  	
		y.clear();
	}
	// the largest value of n is the total number of lines - 2
	_nmax = splines.size() - 1;
	cout << "_nmax: " << _nmax << endl;

	myfile.close();
}

/******************************************************************************
This function returns the height difference function evaluated at 
the input n and r, given the input lambda, eta, and D.

There are three situations:
1) the input n and r is within the range of the utab.dat, so the 
   height difference function gets interpolated. Generally, for 
   small r, the Caille approximation does not work, so the 
   precomputed table is available.

2) the input n or r is outside the range of the utab.dat and 
   the input r is very close to zero (t = r/xi < 0.01), where
   xi = sqrt(lambda * D). In this case, the Caille approximation
   fails, but the interpolation is unavailable, so the full
   integration gets calculated on the fly. This should not 
   reduce the performance of NFIT because this case occurs
   very infrequently, namely only near r = 0.
   
3) the input n or r is outside the range of the utab.dat but
   the input r is sufficietly large enough that the Caille 
   approximation works (t = r/xi > 0.01). 
******************************************************************************/
double Utable::getHeightDiffFunc(int n, double r, double lambda, 
                                 double eta, double D)
{ 
	double ret;
	double t = r / sqrt(lambda * D);
	
	// If t is less than tmax or nmax, then interpolation is available
	if (t < _tmax && n <= _nmax) {
		ret = interp_utable(n, t);
		//ret = gsl_interp_utable(n, t);
	// the input n or t is out of range, so resort to the Caille approx.
	} else {
		// If t is small, the Caille approximation is not accurate,
		// so do the integration on the fly
		if (t < 0.01) {
			//cout << "full blown calculation, r: " << r << " n: " << n << endl;
			ret = calculateIntegration(n, t);
		} else {
			//cout << "Caille, r: " << r << " n: " << n << endl;
			ret = cailleApprox(n, t);
		}
	}
	return ret * eta * D * D / 2 / PI / PI;
}

double Utable::interp_utable(int n, double t)
{
	int ex = 0;
	if (t >= 1) {
		t = log10(t);
		ex = 102;
	}
	int i = int(t/0.01);
	return interp_Newton_2nd(t, 
	                         0.01*i, 0.01*(i+1), 0.01*(i+2), 
	                         utab[i+ex+n*jump], 
	                         utab[i+ex+1+n*jump], 
	                         utab[i+ex+2+n*jump]);
}

double Utable::gsl_interp_utable(int n, double t)
{
	double log10t = log10(t + SMALLNUM);
	return gsl_spline_eval(splines[n], log10t, accs[n]);
}

/******************************************************************************
zeroth order Bessel function of the first kind, 
taken from Numerical Recipes
******************************************************************************/
double bessel_J0(double x)
{
	double ax, z;
	double xx, y, ans, ans1, ans2;

	if ((ax = fabs(x)) < 8.0) {
		y = x * x;
		ans1 = 57568490574.0 + y * (-13362590354.0 + y * (651619640.7 + y 
		* (-11214424.18 + y * (77392.33017 + y * (-184.9052456)))));
		ans2 = 57568490411.0 + y * (1029532985.0 + y * (9494680.718 + y 
		* (59272.64853 + y * (267.8532712 + y * 1.0))));
		ans = ans1 / ans2;
	} else {
		z = 8.0 / ax;
		y = z * z;
		xx = ax - 0.785398164;
		ans1 = 1.0 + y * (-0.1098628627e-2 + y * (0.2734510407e-4 + y 
		* (-0.2073370639e-5 + y * 0.2093887211e-6)));
		ans2 = -0.1562499995e-1 + y * (0.1430488765e-3 + y * (-0.6911147651e-5 +
		y * (0.7621095161e-6 - y * 0.934935152e-7)));
		ans = sqrt(0.636619772 / ax) * (cos(xx) * ans1 - z * sin(xx) * ans2);
	}
	return ans;
}

/******************************************************************************
Returns the exponential integral. Use GSL special function.
******************************************************************************/
double expIntegral(double x)
{
	// for x = 0, exponential integral is not defined
	// for x > 700, gsl_sf_expint_E1 results in underflow error
	return (x > 0 && x < 700) ? gsl_sf_expint_E1(x) : 0;
}

/* 2nd order Newton's divided difference interpolation formula */
double interp_Newton_2nd(double t, 
                         double t0, double t1, double t2,
                         double P0, double P1, double P2)
{
	double P01 = (P1 - P0) / (t1 - t0);
	double P12 = (P2 - P1) / (t2 - t1);
	double P012 = (P12 - P01) / (t2 - t0);
	return P0 + (t-t0)*P01 + (t-t0)*(t-t1)*P012;
}

/* 3rd order Newton's divided difference interpolation formula */
double interp_Newton_3rd(double t, 
                         double t0, double t1, double t2, double t3,
                         double P0, double P1, double P2, double P3)
{
	double P01 = (P1 - P0) / (t1 - t0);
	double P12 = (P2 - P1) / (t2 - t1);
	double P012 = (P12 - P01) / (t2 - t0);
	double P23 = (P3 - P2) / (t3 - t2);
	double P123 = (P23 - P12) / (t3 - t1);
	double P0123 = (P123 - P012) / (t3 - t0);
	return P0 + (t-t0)*P01 + (t-t0)*(t-t1)*P012 + (t-t0)*(t-t1)*(t-t2)*P0123;
}
