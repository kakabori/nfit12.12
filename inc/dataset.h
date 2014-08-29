#ifndef _DataSet_h_
#define _DataSet_h_

#include <tcl.h>
#include <float.h>
#include <vector>
#include "modelcalculator.h"
#include "tvds.h"

/* This performs the 'dataset' command in the UserManual
   Caveat: A qzslice in this subroutine is what we called
   qrslice or qxslice in the thesis.
  The definition is that a slice has qz constant.
*/



struct DataPoint{
  /*
    a data point in a QzSlice.
    it only helps to construct the QzSlice structure
	where qz value is specified
    qx    : qx value
    inte  : the intensity
    sigma : estimated error on the intensity
    cal   : calculated (theoretical) value of the intensity
  */
  double qx, inte, sigma, cal;
};

struct QzSlice{
  /*
    qz: qz value of the slice
    scale: the scaling factor for the slice, will be
	  updated after each fitting iteration
    bias: the background adjustment. (cz parameter in the thesis)
    sdata: the array of data points
   */
  double qz;
  double scale;
  double bias;
  vector<DataPoint> sdata;

  // added by KA on 9/24/2012, used for calculating uncertainties in the scaling factor and cz parameter
  double sigma_scale;
  double sigma_bias;
};


struct Instruction {
  /*
    the instruction on how dataset should be loaded;
    this might be the source of multi-dataset bug,
    I will look into this later.
   */
  int qrstart, qrend;
  double basecnt;
  vector<int> qzstart;
  vector<int> qzend;
  int qrnum(){return qrend-qrstart+1;}
  int qztotnum();
  int qzrange(int &s, int &e);
};

class YFImgBase;
template <class T>
class TVImg;

typedef struct Data dataset;
struct Data {
  /*
    This represents the dataset objects

    qs : array of qzslices
  */
  size_t n;
  vector<QzSlice> qs;
  vector<double> pre;

  bool flag;
  Instruction *ins;

  TVImg<float> *fitimg;
  TVImg<float> *filterimg;

  void readin(YFImgBase *ptr, Instruction*);
  void setimg( TVImg<float> * );
  void filter( TVImg<float> * );
  void print(FILE *fp);
  void print(char *);
  double calvalue(int _qr, int _qz);
  double datvalue(int _qr, int _qz);
  void exportTiff(char *fname);
  void writeimg();
  void writefrm(char *, Para*);
  void info();
  Data(){ ins = NULL;}
};
#endif
