#ifndef GUARD_TVIMG_H
#define GUARD_TVIMG_H


#include <cstdlib>
#include <cstring>
#include <vector>
#include <cmath>
#include <iostream>
#include <typeinfo>
#include <tcl.h>
#include <tk.h>
#include <blt.h>
#include <tiff.h>
#include <tiffio.h>
#include "toad.h"
#include "dataset.h"
#include "boxchain.h"
#include "tvLinFit2d.h"
#include "fileTools.h"

using namespace std;
typedef struct Data dataset;

extern double xlow, xhigh;
extern int descending;
extern int verticalfit;
extern int csix, csiy;

enum DataType {yfUINT16, yfINT16, yfFLOAT, yfINT32};
enum YFImgOpt {yfIMGDEL, yfIMGCREAT};

void funcs_n(float x, float y, float *afunc, int n);
void funcs_1(float x, float y, float *afunc, int n);
void YFTcl_EvalEx(Tcl_Interp* interp, size_t sz, const char*fmt, ...);
extern GlobalData gData;

typedef void* GenericFunc(void*, void*);

union DataPtr{
  void* voidp;
  unsigned short* uint16p;
  short* int16p;
  int* int32p;
  float* floatp;
};

class YFImgBase {
public:
  int width;
  int height;
  int size;
  int type;
  DataPtr dp;
  YFRect<int> rect;
  char plotflag;
  Blt_Vector *xvPtr, *yvPtr;
  char xvname[32], yvname[32];
  char color[16];
  char* name;
  float plotpixels;
  Tcl_Interp* hostInterp;
  vector<GenericFunc*> callbackFuncPtr;
  void setDim(int w, int h){width=w; height=h; size=w*h;}
  void updateCoords();
  int imgsize(){return size;}
  int imgwidth(){return width;}
  int imgheight(){return height;}
  virtual ~YFImgBase(){}
  virtual int openImg(void* fp, int, int) = 0;
  virtual int openImgASCII(const char *) = 0;
  virtual YFImgBase& rotate(double,int,int,double,double,double,double) = 0;
  virtual YFImgBase& nc_rotate90(int,int,int,double,double) = 0;
  virtual YFImgBase& affine(double, double, double, double) = 0;
  virtual void report(FILE *fp) = 0;
  YFRect<int> coords(){return rect;}
  virtual void cleanup(Tcl_Interp*) = 0;
  virtual void ppmBlock(ColorT*, ColorSetup*) = 0;
  virtual void plotCrossSection(int, int, int) = 0;
  virtual void flipx() = 0;
  virtual void flipy() = 0;
  virtual void transpose() = 0;
  virtual void shift(YFDuple<int> s) = 0;
  void shiftposition(YFDuple<int> rv){rect.setll(rect.ll()+rv); updateCoords();}
  void setInterp(Tcl_Interp*);
  void UpdatePlot(int);
  void plotHide(int);
  virtual int exportTiff(char*) = 0;
  virtual void fit2d(Box*, Box*, int, int, int) = 0;
  void setname(char*);
  //void configurePlot(Tcl_Obj *const objv[], int);
  char* getname();
  void show(Tk_PhotoHandle photoHd, ColorSetup* cstp);
  virtual inline float pixData(int, int) = 0;
  virtual YFImgBase& operator+=(double f) = 0;
  virtual YFImgBase& operator*=(double f) = 0;
  virtual YFImgBase& operator+=(YFImgBase& a) = 0;
  virtual YFImgBase& operator*=(YFImgBase& a) = 0;
  virtual YFImgBase& operator-=(YFImgBase& a) = 0;
  virtual YFImgBase& operator/=(YFImgBase& a) = 0;
  virtual double getval(int a, int b) = 0;
  virtual void setval(int a, int b, float v) = 0;
  double getvalsafe(int a, int b);
  virtual void stat(Box*, double*) = 0;
  virtual void resize(int, int) = 0;
  virtual void bgadjust(Data *dp)=0;
  virtual void scadjust(double *sc)=0;
  int getSize();
  void getValue(int j, double *x, double *y);
};

// dummy template struct
template <typename T> struct id { typedef T type; };

template <class T> 
class TVImg : public YFImgBase {
private:
  T* data;
  int readASCII(const char *, id<float>);
  int openImgASCII(const char *, id<float>);
  // g++ complains unless the following dummy methods are defined
  int openImgASCII(const char *, id<int>) { return TCL_ERROR; }
  int openImgASCII(const char *, id<short>) { return TCL_ERROR; }
  int openImgASCII(const char *, id<short unsigned>) { return TCL_ERROR; }
protected:
  static void* geometryChanged(void*, void*);
public:
  TVImg() { initialize(); }
  void initialize();
  
  void plotCrossSection(int x, int y, int opt);
  void ppmBlock(ColorT* pixelPtr, ColorSetup* csetupPtr);
  
  void report(FILE *fp);
  void cleanup(Tcl_Interp*);
  double getval(int a, int b){ return descending==0 ? data[a+width*(height-1-b)] : data[a+width*b];}
  void setval(int a, int b, float v){ data[a+width*b] = (T)v ;}
  double getvalsafe(int a, int b);
  void fit2d(Box*, Box*, int, int, int);
  void stat(Box*, double*);
  void resize(int, int);
  
  TVImg<T>* crop(YFRect<int> rect);
  YFImgBase& operator+=(double f);
  YFImgBase& operator*=(double f);
  YFImgBase& operator+=(YFImgBase& a);
  YFImgBase& operator*=(YFImgBase& a);
  YFImgBase& operator-=(YFImgBase& a);
  YFImgBase& operator/=(YFImgBase& a);

  YFImgBase& affine(double, double, double, double);
  YFImgBase& nc_rotate90(int,int,int,double,double);
  YFImgBase& rotate(double,int,int,double,double,double,double);
  
  void flipx();
  void flipy();
  void shift(YFDuple<int> s);
  void transpose();
  float pixData(int x, int y);
  void shiftcontent(YFDuple<int> rv);
  int exportTiff(char*);
  void bgadjust(Data *dp);
  void scadjust(double *sc);
  int getSize();
  void getValue(int j, double *x, double *y);
  
  int openTiff(void *, int, int);
  int openImgASCII(const char *f) { return openImgASCII(f, id<T>()); }
  int openImg(void *, int, int);
};


/******************************************************************************
Implementation of TVImg<T> methods
******************************************************************************/
//template <class T> TVImg<T>::TVImg()
//{
//  data = NULL;
//  plotflag = 0;
//  xvPtr = yvPtr = NULL;
//  callbackFuncPtr.push_back(TVImg<T>::geometryChanged);
//}

template <class T> void TVImg<T>::initialize()
{
  data = NULL;
  plotflag = 0;
  xvPtr = yvPtr = NULL;
  callbackFuncPtr.push_back(TVImg<T>::geometryChanged);
}

template<class T>
void* TVImg<T>::geometryChanged(void* client, void* instance){
  /* provide callback mechanism */
  TVImg<T>* img=(TVImg<T>*)instance;
  img->updateCoords();
  return NULL;
}

// Only works with TVImg<float>
template <class T>
int TVImg<T>::readASCII(const char *path, id<float>)
{
  vector<float> tmpVec;
  int imgwidth = readMatrixFromFile(path, tmpVec);
  int imglength = tmpVec.size() / imgwidth;

  cout << "imgwidth: " << imgwidth 
       << " imglength: " << imglength << endl;
  setDim(imgwidth, imglength);

	dp.floatp = (float *)data;
	type = yfFLOAT;
	
	// reallocate memory
  free(data);
  data = (float *)malloc(sizeof(float) * tmpVec.size());

  // read the data into data array
  typedef vector<float>::size_type sz;
  for (sz i = 0; i < tmpVec.size(); i++) data[i] = tmpVec[i];

  typedef vector<GenericFunc *>::size_type sz2;
  for(sz2 i = 0; i < callbackFuncPtr.size(); i++)
    callbackFuncPtr[i](NULL, this);

  return TCL_OK;
}

/******************************************************************************
Open an ASCII formatted data
******************************************************************************/
template<class T>
int TVImg<T>::openImgASCII(const char *path, id<float> id)
{
  rect.setll(0,0);
  if( readASCII(path, id) == TCL_OK ) {
    rect.setur(width, height);
    return TCL_OK;
  }
  return TCL_ERROR;
}

/*open an image(tiff) file*/
template<class T>
int TVImg<T>::openImg(void* fp, int optin, int optout)
{  
  rect.setll(0,0);
  if( openTiff(fp, optin, optout) == TCL_OK ){
    rect.setur(width, height);
    return TCL_OK;
  }
  return TCL_ERROR;
}

typedef void (*Fit2DFunc)(float,float,float*,int);

template<class T>
void TVImg<T>::fit2d(Box* box1, Box* box2, int pn, int subopt, int step){
  /*
    Fit the light background using the two box method
    box1: the forbiden area for fitting
    box2: the area to be fitted
    step: the step for 
    pn  : order of the polynomials
  */
  Fit2DFunc fit2dfunc = funcs_n;
  YFRect<int> rect1, rect2, rect3;
  int ma=(pn+1)*(pn+2)/2;
  if(subopt & 1024) {
    ma = pn+1;
    fit2dfunc = funcs_1;
  }
  if(descending == 0){
    rect1.set( box1->xc[0], height-1-box1->yc[0], box1->xc[1], height-1-box1->yc[1] );
    rect2.set( box2->xc[0], height-1-box2->yc[0], box2->xc[1], height-1-box2->yc[1] );
  } else {
    rect1.set( box1->xc[0], box1->yc[0], box1->xc[1], box1->yc[1] );
    rect2.set( box2->xc[0], box2->yc[0], box2->xc[1], box2->yc[1] );
  }
  rect1.shift(rect.ll());
  rect2.shift(rect.ll());
  float* a=(float *)malloc(sizeof(float)*ma);
  if(a==NULL) printf("can not allocate enought memory for a\n");

  rect3=rect2;
  if(verticalfit == 1){
    for(rect3.urx() = rect2.llx(); ; rect3.llx() = rect3.urx() ){
      if( (rect2.urx()-rect3.urx() )<2*step ){
	rect3.urx() = rect2.urx();
	Fit2D_poly(data, width, a,  rect1, rect3, ma, fit2dfunc, pn, subopt);
	return;
      }
      else{
	rect3.urx() += step;
	Fit2D_poly(data, width, a,  rect1, rect3, ma, fit2dfunc, pn, subopt);
      }
    }
  } else {
    for(rect3.ury() = rect2.lly(); ; rect3.lly() = rect3.ury() ){
      if( (rect2.ury()-rect3.ury() )<2*step ){
	rect3.ury() = rect2.ury();
	Fit2D_poly(data, width, a,  rect1, rect3, ma, fit2dfunc, pn, subopt);
	return;
      }
      else{
	rect3.ury() += step;
	Fit2D_poly(data, width, a,  rect1, rect3, ma, fit2dfunc, pn, subopt);
      }
    }
  }
  free(a);
}

template<class T>
void _stat(T *dp, int strip, int w, int h, double *report){
  double sum = 0;
  double mean;
  double sse = 0;
  T max = *dp, min = *dp, t;
  for(int i = 0; i < h; i++)
    for(int j = 0; j < w; j++){
      t = dp[j+i*strip];
      sum += t;
      if(t>max) max = t;
      else if(t<min) min = t;
    }
  report[0] = sum;
  mean = report[1] = sum/w/h;
  report[2] = min;
  report[3] = max;

  for(int i = 0; i < h; i++)
    for(int j = 0; j < w; j++){
      sse += sqr(dp[j+i*strip]-mean);
    }
  cout<<"sse "<<sse<<endl;
  report[4] = sse/w/h;
}

template<class T>
void TVImg<T>::stat(Box* box, double *report){
  /*
    Compute some statistics 
    [pre]
          box: region for computing the statistics
          report: an array to store the result, for now, size of 5
    [post]
          report[0] : sum
          report[1] : average
          report[2] : minimum
          report[3] : maximum  
          report[4] : variance
  */
  YFRect<int> tr;
  if(descending == 0){
    tr.set( box->xc[0], height - 1 - box->yc[0], box->xc[1], height - 1 - box->yc[1] );
  } else {
    tr.set( box->xc[0], box->yc[0], box->xc[1], box->yc[1] );
  }
  tr.shift(rect.ll());
  _stat(data+tr.lly()*width+tr.llx(), width, tr.width(), tr.height(), report);
}

template<class T>
void TVImg<T>::report(FILE *fp){
  fprintf(fp, "width: %d height: %d size: %d=%d\n", width, height, size, width*height);
}

template<class T>
void TVImg<T>::bgadjust(Data *dp){
  /*
    [pre]: dp points to a Data object that has been fitted
    [post]: the backaground adjustment will be made
  */
  for(size_t k=0; k<dp->qs.size(); k++){
    QzSlice* qzsp = &(dp->qs[k]);
    int y = (int)qzsp->qz;
    y = height-1-y;

    for(int x=0; x<width; x++){
      data[x+y*width] += qzsp->bias;
    }
  }
}

template<class T>
void TVImg<T>::scadjust( double *sc){
  /*
    [pre]: the scaling factors from external file
    [post]: the scale adjustment will be made
  */
  for(int y=0; y<height; y++){
    for(int x=0; x<width; x++){
      data[x+y*width] *= sc[height-1-y];
    }
  }
}

template<class T>
int TVImg<T>::exportTiff(char *file){
  TIFF *tifout;
  uint16 bitspersample=16;
  if( !(tifout=TIFFOpen(file,"w")) ) {
    printf("error in openning file\n");
    return TCL_ERROR;
  }
  bitspersample=sizeof(T)*8;
  TIFFSetField(tifout, TIFFTAG_IMAGEWIDTH, width);
  TIFFSetField(tifout, TIFFTAG_IMAGELENGTH, height);
  TIFFSetField(tifout, TIFFTAG_BITSPERSAMPLE, bitspersample);
  //  TIFFSetField(tifout, TIFFTAG_DOCUMENTNAME, "can I put something here");
  TIFFSetField(tifout, TIFFTAG_SUBFILETYPE, 0);
  TIFFSetField(tifout, TIFFTAG_COMPRESSION, 1);
  TIFFSetField(tifout, TIFFTAG_PHOTOMETRIC, 1);
  TIFFSetField(tifout, TIFFTAG_SAMPLESPERPIXEL, 1);
  TIFFSetField(tifout, TIFFTAG_ROWSPERSTRIP, height);
  TIFFSetField(tifout, TIFFTAG_MINSAMPLEVALUE, 0);
  TIFFSetField(tifout, TIFFTAG_MAXSAMPLEVALUE, 65535);
  TIFFSetField(tifout, TIFFTAG_XRESOLUTION, 292.571F);
  TIFFSetField(tifout, TIFFTAG_YRESOLUTION, 292.571F);
  TIFFSetField(tifout, TIFFTAG_RESOLUTIONUNIT, 2);
  TIFFSetField(tifout, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
  // for (int row = 0; row < height; row++) {
  //TIFFWriteScanline(tifout, 0, data+row*width, row, 0);
  //}
  TIFFWriteRawStrip(tifout, 0, data, size*sizeof(T));
  TIFFClose(tifout);
  return TCL_OK;
}


template<class T>
void TVImg<T>::plotCrossSection(int x, int y, int opt){
  /* update the plot in the plot window if requested */
  if(plotflag){
    int asz = yfmax(width, height);
    int dsz;
    if(opt==0){
      dsz = _xplot(data, xvPtr->valueArr, yvPtr->valueArr, width, height, y-rect.lly(), rect.llx(), gData.xswath);
    } else {
      dsz = _yplot(data, xvPtr->valueArr, yvPtr->valueArr, width, height, x-rect.llx(), rect.lly(), gData.yswath);
    }
    Blt_ResetVector(xvPtr, xvPtr->valueArr, dsz, asz, TCL_DYNAMIC);
    Blt_ResetVector(yvPtr, yvPtr->valueArr, dsz, asz, TCL_DYNAMIC);
  }
}

template<class T>
void TVImg<T>::ppmBlock(ColorT* pixelPtr, ColorSetup* csetupPtr) {
  /* draw on pixmap pointed by pixelPtr */
  _ppmline(data, pixelPtr, size, csetupPtr);
}




template<class T>
void TVImg<T>::cleanup(Tcl_Interp* interp){
  /* release resource allocated */  
  YFTcl_EvalEx(hostInterp, 128, "catch {$xyplot element delete %s}", name);
  free(data); data=NULL;
  if(xvPtr!=NULL) Blt_DeleteVector(xvPtr);
  if(yvPtr!=NULL) Blt_DeleteVector(yvPtr);
}

template<class T>
int TVImg<T>::openTiff(void *fp, int optin, int optout){
  TIFF *tif=(TIFF*)fp;
  size_t scanlinesize;
  uint32 imgwidth,imglength;
  uint16 bitspersample;
  int typein = optin;
  int typeout = optout;


  // gather information from the tiff header
  TIFFGetFieldDefaulted(tif, TIFFTAG_IMAGEWIDTH, &imgwidth); 
  TIFFGetFieldDefaulted(tif, TIFFTAG_IMAGELENGTH, &imglength); 
  TIFFGetFieldDefaulted(tif, TIFFTAG_BITSPERSAMPLE, &bitspersample);
  //if( bitspersample/8 != sizeof(T) ) return TCL_ERROR;
  scanlinesize=TIFFScanlineSize(tif);
  cout<<imgwidth<<' '<<imglength<<' '<< bitspersample<<' '<<scanlinesize<<endl;
  setDim(imgwidth, imglength);

  // reallocate memory
  free(data);
  data=(T*)malloc(sizeof(T)*size);
	//data=(int*)malloc(sizeof(int)*size);

	// specify data type; currently all the data are processed as float 
	// in order to keep precision after manipulation on data
	if (typeout == yfFLOAT) {
		dp.floatp = (float*)data;
		type = yfFLOAT;
	} else if (typeout == yfINT16) {
		dp.int16p = (short*)data;	
		type = yfINT16;
	} else if (typeout == yfUINT16) {
		dp.uint16p = (unsigned short*)data;	
		type = yfUINT16;
	} else {
		dp.int32p = (int *)data;
		type = yfINT32;
	}

	// yfFLOAT for images exported by NFIT and TVIEW
	// yfUINT16 for raw data collected by FLICAM CCD at G1 station
  if(typein < 0){ 
    for(uint32 i=0; i<imglength; i++)
      TIFFReadScanline(tif, (char*)data+i*scanlinesize, i, 0);
  }else{
    void *a = malloc(scanlinesize);
    for(uint32 i=0; i<imglength; i++){
      TIFFReadScanline(tif, a, i, 0);
      switch(typein){
      case yfINT16:{   tconvert(imgwidth, (short*)a, 1, data+i*imgwidth, 1); break;}
      case yfFLOAT:{   tconvert(imgwidth, (float*)a, 1, data+i*imgwidth, 1); break;}
      case yfUINT16:{  tconvert(imgwidth, (unsigned short*)a, 1, data+i*imgwidth, 1); break;}
      }
    }
    free(a);
  }

  for(unsigned int i=0; i<callbackFuncPtr.size(); i++)
    callbackFuncPtr[i](NULL, this);

  return TCL_OK;
}


template<class T> 
inline void _ppmline(T *dataPtr, ColorT* pixelPtr, int sz, ColorSetup* csetupPtr){
  /* helper
     draw on pixmap
   */
  float s;
  float sat = csetupPtr->saturation;
  float min = csetupPtr->min;
  float max = csetupPtr->max;
  float bot = csetupPtr->bottom;
  //float fil = csetupPtr->filter;
  ColorT *ct = csetupPtr->colorTable;
  ColorT nc = csetupPtr->negcolor;
  ColorT sc = csetupPtr->satcolor;
  ColorT fc = csetupPtr->filcolor;

  if(csetupPtr->mode==0){
    s = 256 / (max-min);
    for(int k=0; k<sz; k++, dataPtr++){ 
      int dex;
      dex=(int)( (*dataPtr - min) * s );
      if(dex>255)       pixelPtr[k] = (*dataPtr > sat) ? sc : ct[255];
      else if(dex<0)	pixelPtr[k] = (*dataPtr < bot) ? nc : ct[0];
      else              pixelPtr[k] = ct[dex];
      if(*dataPtr ==-32000)	pixelPtr[k] = fc;
    }
  } else if(csetupPtr->mode==1){
    s = 256 / log(max-min);
    for(int k=0; k<sz; k++, dataPtr++){ 
      int dex;
      if( (*dataPtr) <= min ) {
	dex = -1;
      } else {
	dex=(int)( log(*dataPtr - min) * s);
      }
      if(dex>255)       pixelPtr[k] = (*dataPtr > sat) ? sc : ct[255];
      else if(dex<0)	pixelPtr[k] = (*dataPtr < bot) ? nc : ct[0];
      else              pixelPtr[k] = ct[dex];
      if(*dataPtr ==-32000)	pixelPtr[k] = fc;
    }
  }
}

template<class T> 
inline int _xplot(T *ptr, double *xp, double *yp, int width, int height, int y, int xstart, int xswath){
  /* helper
     prepare for plot x crossections
   */
  if (descending == 0) y = height - y -1;
  if (y<xswath || y >= height-xswath) return 0;
  ptr+=width*y;
  int i=0;
  int jlow  = yfmax((int)xlow - xstart, 0);
  int jhigh = yfmin((int)xhigh - xstart, width);
  for(int j=jlow; j<jhigh; j++, i++){
    yp[i]=0;
    for(int k=-xswath;k<=xswath;k++){
      yp[i]+=ptr[j+k*width];
    }
    yp[i]/=(1+2*xswath);
    xp[i]=xstart+j;
  }
  return jhigh-jlow;
}

template<class T> inline 
int _yplot(T *ptr, double *xp, double *yp, int width, int height, int x, int ystart, int yswath){
  /* helper
     prepare for plot y crossections
   */
  if (x<yswath ||x >= width-yswath) return 0;
  ptr+=x;
  int i=0;
  int jlow  = yfmax((int)xlow - ystart, 0);
  int jhigh = yfmin((int)xhigh - ystart, height);
  if(descending == 0) {
    for(int j=jlow; j<jhigh; j++, i++){
      yp[i]=0;
      for(int k=-yswath;k<=yswath;k++){
	yp[i]+=ptr[k+(height-j-1)*width];
      }
      yp[i]/=(1+2*yswath);
      xp[i]=ystart+j;
    }
  } else {
    for(int j=jlow; j<jhigh; j++, i++){
      yp[i]=0;
      for(int k=-yswath;k<=yswath;k++){
	yp[i]+=ptr[k+j*width];
      }
      yp[i]/=(1+2*yswath);
      xp[i]=ystart+j;
    }
  }
  return jhigh-jlow;
}


template<class T>
TVImg<T>* TVImg<T>::crop(const YFRect<int> rect){
  /* not finished */
  TVImg<T>* img=new TVImg<T>;
  return img;
}

template<class T>
YFImgBase& TVImg<T>::operator+=(double f){
  /* add constant */
  yf_addConst(data, f, size);
  return *this;
}

template<class T>
YFImgBase& TVImg<T>::operator*=(double f){
  /* multiply by a scalar */
  yf_multiConst(data, f, size);
  return *this;
}


template<class T>
YFImgBase& TVImg<T>::operator+=(YFImgBase& a){
  /* pointwise addition */
  YFRect<int> tmprect=rect;
  YFRect<int> recta=a.coords();
  YFRect<int> rectan=rect.overlap(a.rect);
  int x1=rectan.llx();
  int y1=rectan.lly();
  int x2=rectan.urx();
  int y2=rectan.ury();
  int offset1=width*(tmprect.ury()-y2)+(x1-tmprect.llx());
  int offset2=a.width*(recta.ury()-y2)+(x1-recta.llx());
  cout<<offset1<<" "<<offset2<<endl;
  switch(type){
  case yfFLOAT:
    if(a.type==yfFLOAT) nk_addArr(dp.floatp+offset1, width, a.dp.floatp+offset2, a.width, x2-x1, y2-y1);
     else if(a.type==yfINT16) nk_addArr(dp.floatp+offset1, width, a.dp.int16p+offset2, a.width, x2-x1, y2-y1);
    break;
  case yfINT16:
    if(a.type==yfFLOAT) nk_addArr(dp.int16p+offset1, width, a.dp.floatp+offset2, a.width, x2-x1, y2-y1);
     else if(a.type==yfINT16) nk_addArr(dp.int16p+offset1, width, a.dp.int16p+offset2, a.width, x2-x1, y2-y1);
    break;
  }
  rect=tmprect;
  return *this;
}

template<class T>
YFImgBase& TVImg<T>::operator-=(YFImgBase& a){
  /* pointwise subtraction */
  YFRect<int> tmprect=rect;
  YFRect<int> recta=a.coords();
  YFRect<int> rectan=rect.overlap(a.rect);
  int x1=rectan.llx();
  int y1=rectan.lly();
  int x2=rectan.urx();
  int y2=rectan.ury();
  int offset1=width*(tmprect.ury()-y2)+(x1-tmprect.llx());
  int offset2=a.width*(recta.ury()-y2)+(x1-recta.llx());
  cout<<offset1<<" "<<offset2<<endl;
  switch(type){
  case yfFLOAT:
    if(a.type==yfFLOAT) nk_subArr(dp.floatp+offset1, width, a.dp.floatp+offset2, a.width, x2-x1, y2-y1);
     else if(a.type==yfINT16) nk_subArr(dp.floatp+offset1, width, a.dp.int16p+offset2, a.width, x2-x1, y2-y1);
    break;
  case yfINT16:
    if(a.type==yfFLOAT) nk_subArr(dp.int16p+offset1, width, a.dp.floatp+offset2, a.width, x2-x1, y2-y1);
     else if(a.type==yfINT16) nk_subArr(dp.int16p+offset1, width, a.dp.int16p+offset2, a.width, x2-x1, y2-y1);
    break;
  }
  rect=tmprect;
  return *this;
}

template<class T>
YFImgBase& TVImg<T>::operator*=(YFImgBase& a){
  /* pointwise multiplication */
  YFRect<int> tmprect=rect;
  YFRect<int> recta=a.coords();
  YFRect<int> rectan=rect.overlap(a.rect);
  int x1=rectan.llx();
  int y1=rectan.lly();
  int x2=rectan.urx();
  int y2=rectan.ury();
  int offset1=width*(tmprect.ury()-y2)+(x1-tmprect.llx());
  int offset2=a.width*(recta.ury()-y2)+(x1-recta.llx());
  cout<<offset1<<" "<<offset2<<endl;
  switch(type){
  case yfFLOAT:
    if(a.type==yfFLOAT) nk_multiArr(dp.floatp+offset1, width, a.dp.floatp+offset2, a.width, x2-x1, y2-y1);
     else if(a.type==yfINT16) nk_multiArr(dp.floatp+offset1, width, a.dp.int16p+offset2, a.width, x2-x1, y2-y1);
    break;
  case yfINT16:
    if(a.type==yfFLOAT) nk_multiArr(dp.int16p+offset1, width, a.dp.floatp+offset2, a.width, x2-x1, y2-y1);
     else if(a.type==yfINT16) nk_multiArr(dp.int16p+offset1, width, a.dp.int16p+offset2, a.width, x2-x1, y2-y1);
    break;
  }
  rect=tmprect;
  return *this;
}

template<class T>
YFImgBase& TVImg<T>::operator/=(YFImgBase& a){
  /* pointwise division */
  YFRect<int> tmprect=rect;
  YFRect<int> recta=a.coords();
  YFRect<int> rectan=rect.overlap(a.rect);
  int x1=rectan.llx();
  int y1=rectan.lly();
  int x2=rectan.urx();
  int y2=rectan.ury();
  int offset1=width*(tmprect.ury()-y2)+(x1-tmprect.llx());
  int offset2=a.width*(recta.ury()-y2)+(x1-recta.llx());
  cout<<offset1<<" "<<offset2<<endl;
  switch(type){
  case yfFLOAT:
    if(a.type==yfFLOAT) nk_divArr(dp.floatp+offset1, width, a.dp.floatp+offset2, a.width, x2-x1, y2-y1);
     else if(a.type==yfINT16) nk_divArr(dp.floatp+offset1, width, a.dp.int16p+offset2, a.width, x2-x1, y2-y1);
    break;
  case yfINT16:
    if(a.type==yfFLOAT) nk_divArr(dp.int16p+offset1, width, a.dp.floatp+offset2, a.width, x2-x1, y2-y1);
     else if(a.type==yfINT16) nk_divArr(dp.int16p+offset1, width, a.dp.int16p+offset2, a.width, x2-x1, y2-y1);
    break;
  }
  rect=tmprect;
  return *this;
}

template<class T>
YFImgBase& TVImg<T>::rotate(double rad, int nw, int nh, double cx, double cy, double ncx, double ncy){
  yf_rotateImg(&data, rad, width, height, nw, nh, cx, cy, ncx, ncy);
  size = nw * nh;
  width = nw;
  height = nh;
  if( sizeof(T) == sizeof(float) ){
    dp.floatp = (float*)data;
    type = yfFLOAT;
  } else {
    dp.int16p = (short*)data;
    type = yfINT16;
  }
  if(xvPtr != NULL )  Blt_ResizeVector (xvPtr, width);
  if(yvPtr != NULL )  Blt_ResizeVector (yvPtr, height);
  for(unsigned int i=0; i<callbackFuncPtr.size(); i++)
    callbackFuncPtr[i](NULL, this);
  return *this;
}

template<class T>
YFImgBase& TVImg<T>::nc_rotate90(int direction, int nw, int nh, double cx, double cy){
  nc_rotateImg(&data, direction, nw, nh, cx, cy);
  size = nw * nh;
  width = nw;
  height = nh;
  if( sizeof(T) == sizeof(float) ){
    dp.floatp = (float*)data;
    type = yfFLOAT;
  } else {
    dp.int16p = (short*)data;
    type = yfINT16;
  }
  if(xvPtr != NULL )  Blt_ResizeVector (xvPtr, width);
  if(yvPtr != NULL )  Blt_ResizeVector (yvPtr, height);
  for(unsigned int i=0; i<callbackFuncPtr.size(); i++)
    callbackFuncPtr[i](NULL, this);
  return *this;
}

/* affine transformation of the image */
template <class T>
YFImgBase& TVImg<T>::affine(double rad, double scale, double tx, double ty)
{
  yf_affineImg(&data, rad, scale, width, height, width, height, tx, ty);
  if( sizeof(T) == sizeof(float) ){
    dp.floatp = (float*)data;
    type = yfFLOAT;
  } else {
    dp.int16p = (short*)data;
    type = yfINT16;
  }
  for(unsigned int i=0; i<callbackFuncPtr.size(); i++)
    callbackFuncPtr[i](NULL, this);
  return *this;
}



template<class T>
void TVImg<T>::flipx(){
  yf_flipx(data, width, height);
}

template<class T>
void TVImg<T>::flipy(){
  yf_flipy(data, width, height);
}

template<class T>
void TVImg<T>::shift(YFDuple<int> s){
  nk_shift(data, s, width, height);
}

template<class T>
void TVImg<T>::transpose(){
  yf_transpose(data, width, height);
  yf_swap(width, height);
  for(unsigned int i=0; i<callbackFuncPtr.size(); i++)
    callbackFuncPtr[i](NULL, this);
}

template<class T>
float TVImg<T>::pixData(int x, int y){
  /* return value at (x,y) */
  return data[x+y*width];
}

template<class T>
void TVImg<T>::resize(int w, int h){
  setDim(w, h);
  data = (T*)realloc(data, size*sizeof(T));
  if( sizeof(T) == sizeof(float) ){
    dp.floatp = (float*)data;
    type = yfFLOAT;
  } else {
    dp.int16p = (short*)data;
    type = yfINT16;
  }
  rect.setll(0,0);
  rect.setur(width, height);
}


#endif
