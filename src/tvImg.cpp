#include <cstdlib>
#include <cstring>
#include <tcl.h>
#include <tk.h>
#include <float.h>

#include "tvImg.h"
#include "tvLinFit2d.h"

void YFImgBase::updateCoords(){
  YFDuple<int> wh(width, height);
  rect.setur(rect.ll()+wh);
}

void YFImgBase::setInterp(Tcl_Interp* interp){
  hostInterp = interp;
}


void YFImgBase::UpdatePlot(int opt){
  /* delete or create a curve in the plot window */
  if(opt==0){
    YFTcl_EvalEx(hostInterp, 512, "catch {$xyplot element delete %s}",name);
  } else if (opt==1) {
    YFTcl_EvalEx(hostInterp, 1024,
		 "eval $xyplot element create %s -xdata %s -ydata %s [gPopt %s] [cPopt %s]" // \\{ %s\\}
		 , name, xvname, yvname, name, name);
  }
}


void YFImgBase::plotHide(int opt){
  /* this toggles the plot */
  if(opt==0){
    if(plotflag==1){
      UpdatePlot(0);
      if(xvPtr!=NULL) Blt_DeleteVector(xvPtr);
      if(yvPtr!=NULL) Blt_DeleteVector(yvPtr);
      xvPtr=yvPtr=NULL;
    }
    plotflag=0;
  } else if(opt==1){
    if(plotflag==0){
      int asz = yfmax(width, height);
      sprintf(xvname, "xv%ld", (long unsigned int)this);
      Blt_CreateVector(hostInterp, xvname, asz, &xvPtr);
      sprintf(yvname, "yv%ld", (long unsigned int)this);
      Blt_CreateVector(hostInterp, yvname, asz, &yvPtr);
    }
    UpdatePlot(1);
    plotflag=1;
    plotCrossSection(csix, csiy, gData.plotopt);
  } else if(opt==2){
    if(plotflag==1) plotHide(0);
    else if(plotflag==0) plotHide(1);
  }
}

void YFImgBase::setname(char* n){
  name=n;
}


char* YFImgBase::getname(){
  return name;
}

void YFImgBase::show(Tk_PhotoHandle photoHd, ColorSetup* cstp){
  /*
    photoHd : image handle of the image window
    cstp:     the color scheme to be used
    [post]: the image on the screen will be updated
  */

  Tk_PhotoImageBlock block;
  if(photoHd==NULL) return;

  // get the pixmap of the PhotoHandle
  Tk_PhotoSetSize(photoHd, width, height);
  Tk_PhotoGetImage(photoHd,&block);

  // draw this image on the pixmap
  ppmBlock((ColorT*)block.pixelPtr, cstp);

  // notify tk that the image has been changed
  Tk_ImageChanged( *((Tk_ImageMaster*)photoHd), 0, 0, width, height, width, height);

  // invoke a tcl/tk line to show the image
  YFTcl_EvalEx(hostInterp, 512, "$imgWN marker configure bg -coords {%d %d %d %d}",
	       rect.llx(), rect.lly(), rect.urx(), rect.ury() );
}


void funcs_n(dfprec x, dfprec y, dfprec *afunc, int n){
  int i,k;
  x*=1e-3f;
  y*=1e-3f;
  *afunc=1;
  for(i=0; i<n; i++, afunc++) *(afunc+1) = (*afunc) * x;
  for(k=n+1; k>1; k--){
    for(i=1; i<k; i++){
      afunc++;
      *afunc = *(afunc-k) * y;
    }
  }
}

void funcs_1(dfprec x, dfprec y, dfprec *afunc, int n){
  int i;
  x*=1e-3f;
  y*=1e-3f;
  *afunc=1;
  for(i=0; i<n; i++, afunc++) *(afunc+1) = (*afunc) * x;
}

double YFImgBase::getvalsafe(int a, int b){
  /* return value at (a,b) with boundary checking */
  a -= rect.llx();
  b -= rect.lly();
  if(b>=0 && b<height && a>=0 && a<width) return getval(a,b);
  return DBL_MIN;
}

void ColorSetup::buildColorTable( ColorStyle opt ){
  /* given a style, build a color table*/
  free(colorTable);
  colorTable=(ColorT*)malloc(sizeof(ColorT)*256);
  switch(opt){
  case GRAYSCALE1:
    for(int j=0;j<256;j++) colorTable[j].set(j,j,j,255);
    break;
  case GRAYSCALE2:
    for(int j=0;j<256;j++) colorTable[255-j].set(j,j,j,255);
  }
  negcolor.set(255,0,0,255);
  satcolor.set(0,255,0,255);
  filcolor.set(0,0,255,255);
  printf("color table built\n");
}

int YFImgBase::getSize(){
  if (plotflag==1) {
    return xvPtr->numValues;
  } else return 0;
}

void YFImgBase::getValue(int j, double *x, double *y){
  *x=xvPtr->valueArr[j];
  *y=yvPtr->valueArr[j];
}



