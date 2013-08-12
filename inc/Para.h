#ifndef GUARD_PARA_H
#define GUARD_PARA_H

#include <tcl.h>
#include "tvds.h"
#include "nfit.h"

// The total number of fitting parameters, defined in Para
#define TOT_NUM_FIT_PARAMS 15

enum Var stringToVarEnum(const std::string &);

/****************************************************************************************
parameter class
****************************************************************************************/
struct Para {
public:
  double Kc; // bending modulus
  double B; // bulk modulus
  double avgLr; // average domain size in-plane
  double avgMz; // average domain size out-of-plane
  double D; // D-spacing
	double mosaic; // mosaic spread
  double edisp; // energy dispersion and beam divergence
  double bFWHM; // beam full width half maximum in pixel
  TvDs2 setup;   // piece borrowed from DS
  double T; // temperature
  double beamSigma; // beamFWHM in q-space (inverse Angstrom)

  double lambda; // fluctuation parameter
  double eta; // fluctuation parameter

  static int s_numParams; // the total number of fitting parameters
  double epsfcn[TOT_NUM_FIT_PARAMS];
  int idx[TOT_NUM_FIT_PARAMS];
  double *xp[TOT_NUM_FIT_PARAMS];
  int nfit; // the number of free parameters for a current fit

  void setNfit(Tcl_Interp *interp, Tcl_Obj *const objv[], int n);
  void _xp(int ix, double **p, double *dp);
  bool _check(int ix);
  void setValue(double, int);
  void setValue(double, const char *);
  void setValue(double, const std::string&);
  void setLambdaEta();
  double getValue(int);
  void set_beamSigma();
};

#endif
