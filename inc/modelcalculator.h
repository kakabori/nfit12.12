#ifndef GUARD_MODELCALCULATOR_H
#define GUARD_MODELCALCULATOR_H

#include <math.h>
#include <gsl/gsl_integration.h>
//#include <tcl.h>
#include <vector>
#include <string>
#include "Para.h"
#include "tvds.h"
#include "nfit.h"
#include "funsupport.h"
#include "utable.h"
//#include "alglib/src/stdafx.h"
//#include "alglib/src/interpolation.h"

//using namespace alglib;
#define PI 3.1415926535897932384626433832795

double Pr(double, double);
double finiteSizeFactor(double);
double bessel_J0(double x);
double getLowerLimit(double, double, double);
double getUpperLimit(double, double, double);
enum Var stringToVarEnum(const std::string &);

/****************************************************************************************
ModelCalculator class: 
This is the class that calculates the theoretical structure factor.
FuncLmdif should bind to this class object for chi square minimization
****************************************************************************************/
class ModelCalculator {
public:
  ModelCalculator();
  //This constructor provides a means of cloning a ModelCalculator object
  ModelCalculator(const ModelCalculator& mc);
  ~ModelCalculator(){cleanup();}
  void cleanup();
  
  // setters and getters
  void setModelParameter(double, const std::string&);
  void setModelParameter(double, const char *);
  void setpara(double, int);
  void paraSet(Para *p);

  // Utable utilities
  void read_in_utable(const char *f){utable.readUtableFile(f);}
  
  // FunSupport utilities
  void resetTOL();

	// functions to build the interpolating function for the structure factor
  void setSliceParameter(double qz);
  void QxSlice(double qxlow, double qxhigh);
  double getCCDStrFct(double);

	// core functions that calculate the structure factor
	static double s_rotatedWrapper(double, void *);
	static double s_qyIntegrandWrapper(double, void *);
	double StrFct(double);
  double sumn(double);
  double Hr(double);
  double Hz(int);
  static double s_StrFctWrapper(double, void *);
  static double s_Fr(double, void *);
  static double s_sumnWrapper(double, void *);
  static double s_HrWrapper(double, void *);
  static double s_HrIntegrand(double, void *);
	double rotated(double);
	double beamConvolutedStrFct(double);
	static double s_convIntegrand(double, void *);
	double beamProfile(double);
  
  // a bunch of functions that are useful in debugging
  void get_spStrFctPoints(std::vector<double>&, std::vector<double>&);
  void get_spsumnPoints(std::vector<double>&, std::vector<double>&);
  void get_spHrPoints(std::vector<double>&, std::vector<double>&);
  void get_spRotatedPoints(std::vector<double>&, std::vector<double>&);
  void getStrFct(double, double, double, std::vector<double>&, std::vector<double>&);
  void getSumn(double, double, double, std::vector<double>&, std::vector<double>&);
  void getHr(double, double, std::vector<double>&, std::vector<double>&);
  void getRotatedSF(double, double, double, std::vector<double>&, std::vector<double>&);
private:
  double Kc, B, dspacing, T;
  double avgLr, avgMz;
	double edisp;
	double mosaic;
  double wavelength, pixelSize, bFWHM, sdistance, beamSigma;

  // cutoff for r integration and sum over n
  double cutoff_r, cutoff_n;
  // scaling constants necessary to get value from Utable
  double in, out;
	// qz value for the qx slice
	double qz;
  // qx lower and upper limits for CCD structure factor interpolant
  double qxHighLimit, qxLowLimit;
  // gsl workspace
  gsl_integration_workspace *workspace;
  // interpolation support object for sum over n
  FunSupport spsumn;
  // interpolation support object for the structure factor
  FunSupport spStrFct;
  // interpolation support object for effective finite size factor in r
  FunSupport spHr;
	// interpolation support object for rotated structure factor
  FunSupport spRotated;
  Utable utable;
  double curr_r, currQr, currQx, currQy;
	// should be called after parameters are changed. considered as updater
  void avgLrChanged();
  void avgMzChanged();
  void KcBDTChanged();
  void set_beamSigma();
};

#endif
