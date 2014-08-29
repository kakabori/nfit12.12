#ifndef _TVMAIN_HEADER_
#define _TVMAIN_HEADER_

#include <vector>
#include <stdarg.h>
#include "boxchain.h"

using namespace std;

class YFImgBase;
enum ColorStyle { GRAYSCALE1, GRAYSCALE2 };

class ColorT{
  unsigned char pix[4];
public:
  void set(unsigned char p0, 
	   unsigned char p1, 
	   unsigned char p2, 
	   unsigned char p3){pix[0]=p0; pix[1]=p1; pix[2]=p2; pix[3]=p3;}
};

struct ColorSetup{
  ColorSetup(){ colorTable=NULL; min=max=saturation=bottom=filter=0;mode=1;}
  ~ColorSetup(){ free(colorTable); }
  ColorT *colorTable;
  ColorT negcolor;
  ColorT satcolor;
  ColorT filcolor;
  float min, max;
  // float scaling;
  float saturation;
  float bottom;
  float filter;
  int mode;
  void buildColorTable(ColorStyle opt);
};


struct GlobalData{
  Tcl_Interp *interp;
  Tk_PhotoHandle tv_photoHd;
  ColorSetup csetup;
  int plotopt;
  int xswath;
  int yswath;
  vector< YFImgBase* > imgPtrVec;
  Tcl_HashTable imgNT;
  //Graph* graphPtr;
  GlobalData(){Tcl_InitHashTable(&imgNT, TCL_STRING_KEYS);}
};

typedef struct PhotoMaster{
  Tk_ImageMaster tkMaster;
  Tcl_Interp *interp;
} PhotoMaster;

#endif
