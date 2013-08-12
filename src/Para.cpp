#include <tcl.h>
#include <iostream>

#include "Para.h"
#include "nfit.h"

using namespace std;
#define PI 3.1415926535897932384626433832795
int Para::s_numParams = TOT_NUM_FIT_PARAMS;

/****************************************************************************************
This function sets a model parameter to an input value.
A model parameter is specified by the name string.
See ModelCalculator::stringToVarEnum for a list of strings that can be
passed as an input name.

When a new parameter needs to be added to the model, add a new case,
add a new if statement in ModelCalculator::stringToVarEnum, and
add a new enum value in enum Var

value: an input value
name: an input name specifying which parameter to be set to an input value
****************************************************************************************/
void Para::setValue(double value, const char *_name)
{
	string name(_name);
	setValue(value, name);
}

void Para::setValue(double value, const string& name)
{
  setValue(value, stringToVarEnum(name));
}

/****************************************************************************************
an accessor function
Given an input index, a corresponding variable gets set to an input value
When either Kc, B, D, or T gets modifiedy, lambda and eta get recalculated.
This function accepts int instead of enum itself so that a user of this object
does not have to define a separate enum.
****************************************************************************************/
void Para::setValue(double value, int index)
{
  switch (index) {
    case Var_Kc:                       Kc = value;  setLambdaEta(); break;
    case Var_B:                         B = value;  setLambdaEta(); break;
    case Var_avgLr:                 avgLr = value;                  break;
    case Var_avgMz:                 avgMz = value;                  break;
    case Var_D:                         D = value;  setLambdaEta(); break;
    case Var_mosaic:               mosaic = value;                  break;
    case Var_edisp:                 edisp = value;                  break;
    case Var_bFWHM:                 bFWHM = value; set_beamSigma(); break;
    case Var_sdistance:           setup.s = value; set_beamSigma(); break;
    case Var_bc2b:             setup.bc2b = value;                  break;
    case Var_wavelength: setup.wavelength = value; set_beamSigma(); break;
    case Var_pixelSize:          setup.pz = value; set_beamSigma(); break;
    case Var_qxzero:         setup.qrzero = value;                  break;
    case Var_nindex:         setup.nindex = value;                  break;
    case Var_T:                         T = value;  setLambdaEta(); break;
    case Var_beamSigma:         beamSigma = value;                  break;
    case Var_WrongInput: cerr << "Var_WrongInput at Para::setValue" << endl; break;
    default: cerr << "wrong index input to Para_setValue" << endl; break;
  }
}

double Para::getValue(int index)
{
  double value = 0;
  switch (index) {
    case Var_Kc:         value = Kc; break;
    case Var_B:          value = B; break;
    case Var_avgLr:      value = avgLr; break;
    case Var_avgMz:      value = avgMz; break;
    case Var_D:          value = D; break;
    case Var_mosaic:     value = mosaic; break;
    case Var_edisp:      value = edisp; break;
    case Var_bFWHM:      value = bFWHM; break;
    case Var_sdistance:  value = setup.s; break;
    case Var_bc2b:       value = setup.bc2b; break;
    case Var_wavelength: value = setup.wavelength; break;
    case Var_pixelSize:  value = setup.pz; break;
    case Var_qxzero:     value = setup.qrzero; break;
    case Var_nindex:     value = setup.nindex; break;
    case Var_T:          value = T; break;
    case Var_beamSigma:  value = beamSigma; break;
    case Var_WrongInput: cerr << "Var_WrongInput at Para::getValue" << endl; break;
    default: cerr << "wrong index input to Para_getValue" << endl; break;
  }
  return value;
}

/****************************************************************************************
This function updates values of lambda and eta.
Should be called whenever a value of Kc, B, D, or T gets modified
****************************************************************************************/
void Para::setLambdaEta()
{
  lambda = 1e16 * sqrt(Kc/B) / D;
  eta = 2.17 * (T+273.15) / sqrt(Kc*B) / D / D;
}

void Para::set_beamSigma()
{
	double deltaQPerPixel = 2 * PI * setup.pz / setup.wavelength / setup.s;
	beamSigma = deltaQPerPixel * bFWHM;
}

/****************************************************************************************
If Kc, B, lxctr, and mctr are free parameters, tcl procedure named setnfit executes
paraset nfit p 0 1 2  4
objv points to the fourth item, 0
n = 4, which is the number of free parameters
****************************************************************************************/
void Para::setNfit (Tcl_Interp *interp, Tcl_Obj *const objv[], int n)
{
  int ix;
  nfit = 0;
  //cout << "n =" << n << endl;
  while ( --n >= 0 ) {
    if ( Tcl_GetIntFromObj(interp, objv[n], &ix) != TCL_OK ) continue;
    //cout << "ix = " << ix << endl;
    //cout << "n = " << n << endl;
    if ( ! _check(ix) ) continue;
    //cout << "nfit = " << nfit << endl;
    idx[nfit] = ix;
    _xp(ix, xp + nfit, epsfcn + nfit);
    nfit += 1;
  }
}

bool Para::_check (int ix)
{
  for (int i = 0; i < nfit; i++) {
    if (ix == idx[i]) return false;
    else continue;
  }
  return ix < Para::s_numParams;
}

void Para::_xp (int ix, double **p, double *dp){
  double *ptr = (double*)this;
  *p = ptr + ix;
  *dp = 0.01;
  if(ix == Var_D)          *dp = 0.0001;
  else if(ix == Var_bc2b)  *dp = 0.001;
}

/******************************************************************************
This function converts an input string to a correspoinding enum value.
See the definition of enum Var in nfit.h

inString: an input string
******************************************************************************/
Var stringToVarEnum(const string& inString)
{
  if (inString == "Kc") return Var_Kc;
  else if (inString == "B") return Var_B;
  else if (inString == "D") return Var_D;
  else if (inString == "T") return Var_T;
  else if (inString == "Lr") return Var_avgLr;
  else if (inString == "Mz") return Var_avgMz;
	else if (inString == "edisp") return Var_edisp;
	else if (inString == "mosaic") return Var_mosaic;
  else if (inString == "wavelength") return Var_wavelength;
  else if (inString == "pixelSize") return Var_pixelSize;
	else if (inString == "bFWHM") return Var_bFWHM;
	else if (inString == "s") return Var_sdistance;
	else if (inString == "bc2b") return Var_bc2b;
	else if (inString == "qxzero") return Var_qxzero;
	else if (inString == "nindex") return Var_nindex;
	else if (inString == "beamSigma") return Var_beamSigma;
  return Var_WrongInput;
}
