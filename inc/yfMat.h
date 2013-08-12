#ifndef _YFMAT_HEADER_
#define _YFMAT_HEADER_

/*
YFMatView

immature trial on creating views.
a view means a rectangular area of
an image
*/

struct YFMatView{
  DataPtr dp;
  int offset;
  int strip;
  //int step; 1
  int w, h;
  int type;
  void medianFilterInplace();
  void medianFilter();
  void peakCenter();
  void convfilter();
  void FFTsmooth(double*);
};

#endif
