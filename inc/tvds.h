#ifndef _TVDS_H_

#define _TVDS_H_



#include "mydll.h"



class TvDs2

{

 public:

  /*

    parameters related with optical setup

    They are defined in the thesis.

  */

  double s;

  double bc2b;// x0, beeker center to bottom of ccd

  double alpha;

  double R;

  double wavelength;

  double pz;

  double qrzero;

  double nindex;


  /* in-line pz-->qz */
  double getqz(double zpix){return 4*YFPI/wavelength*getSinTheta(zpix);}

  double getqr(double rpix);

  double getPixQr(double qr);

  double getl(double, double);

  double getSinTheta(double);

};





#endif

