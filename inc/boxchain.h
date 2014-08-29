#ifndef _BOXCHAIN_HEADER_
#define _BOXCHAIN_HEADER_

#include "toadmisc.h"
#include <vector>

using namespace std;

enum Boption {RESIZEBOX, MOVEBOX, NEWBOX};

struct Box{
  /*
    We do not specify top/bottom, left/right of the box,
    because the axes of BLT graph window ($imgWN) in which
    the boxes are drawn can be descending or ascending.
   */
  int xc[2]; /* the two x values define the box */
  int yc[2]; /* the two y values define the box */
  char name[4]; /* name of the box */
  inline double nearx(double);
  inline double neary(double);
  inline int width(){return abs(xc[0]-xc[1]);}
  inline int height(){return abs(yc[0]-yc[1]);}
  inline double midxc(){return (xc[0]+xc[1])/2.0;} /* x of the center */
  inline double midyc(){return (yc[0]+yc[1])/2.0;} /* y of the center */
  inline int minxc(){return yfmin(xc[0],xc[1]);}   /* the smaller value in xc[] */
  inline int minyc(){return yfmin(yc[0],yc[1]);}   /* the smaller value in yc[] */
  inline void set(int x0, int y0, int x1, int y1){xc[0]=x0; xc[1]=x1; yc[0]=y0; yc[1]=y1;}
  inline void nearcorner(double x, double y, int& xi, int& yi);
};

inline double Box::nearx(double x){
  /* the distance from x to xc[] */
  double d0=fabs(x-xc[0]);
  double d1=fabs(x-xc[1]);
  return (d0<d1)?d0:d1;
}
inline double Box::neary(double y){
  /* the distance from y to yc[] */
  double d0=fabs(y-yc[0]);
  double d1=fabs(y-yc[1]);
  return (d0<d1)?d0:d1;
}
inline void Box::nearcorner(double x, double y, int& xi, int& yi){
  /*
    find the corner (xi, yi) close to point (x,y)
    xi, yi are the indices correspond to xc[] and yc[]
  */
  double d0=fabs(x-xc[0]);
  double d1=fabs(x-xc[1]);
  xi=(d0<d1)?0:1;
  d0=fabs(y-yc[0]);
  d1=fabs(y-yc[1]);
  yi=(d0<d1)?0:1;
}

class BoxVec{
  Tcl_Interp* itp; /* interpreter for callback */
  Tcl_HashTable bvhash; /* Hash table that stores Box pointers, name of boxes are the keys*/
  double xhalo, yhalo;
  Box* curb;
  int option;
  int xi, yi, pinx, piny;
  int nextname;
  vector<int> recyclename;
 public:
  void step1(double,double,char*);
  void step2(int,int);
  void step3();
  void sethalo(double x, double y){xhalo=x; yhalo=y;}
  void init(Tcl_Interp* interp);
  void newbox(int, int, int, int, char*);
  void del(char*);
  Box* boxPtr(char*);

  void updateBox(Box* bptr);
};

#endif
