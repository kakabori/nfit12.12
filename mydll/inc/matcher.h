#ifndef __MATCHER_HH
#define __MATCHER_HH

#include "mydll.h"

class Matcher {
public:
  cmatrix aut;
  cmatrix w;
  cmatrix vt;
  cmatrix b;

  cmatrix autV;
  int m;
  


  Matcher(int maxm);
  void  svdMatch(cmatrix &x);
  double  dfpMatch(cmatrix &x);
  double  match(cmatrix &x);
  void  match2(cmatrix &x, double thresh, void**);
  void  removeOutlier(cmatrix &x, double thresh,void**);
  void  match2(cmatrix &x, double thresh, int*);
  void  removeOutlier(cmatrix &x, double thresh,int*);
  void  push(double rx, double ry, double x, double y);
  void  setsize(int maxm);
  void  reset();
  double  error(double theta, double dx, double dy);
  void  derror(double theta, double dx, double dy , double *df);


};

#endif
