#include <cstdlib>
#include <cstring>
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <vector>
#include <tcl.h>
#include <tk.h>
#include <blt.h>
#include <tiff.h>
#include <tiffio.h>
#include <semaphore.h>
#include <pthread.h>
#include <math.h>
#include <stdio.h>
#include <ctype.h>
#include <time.h>
#include <unistd.h>

#include "toadcmd.h"
#include "tvImg.h"
#include "toad.h"
#include "toadmisc.h"
#include "dataset.h"
#include "modelcalculator.h"
#include "funsupport.h"
#include "tvds.h"
#include "funclmdif.h"
#include "mydll.h"
#include "nfit.h"
#include "globalVariables.h"
#include "utable.h"
#include "Para.h"


using std::cout; using std::endl;

extern "C" int Blt_Init(Tcl_Interp *interp);

int FLICAM = 0;
extern int linkxy;
extern int nocz, noscale, kiyomask;
extern double NKfilter;
int descending = 0;
int verticalfit = 0;
double tmpsigma = 1e-7;
extern double dupe;

extern double aFactor;
extern double backgroundSigmaSquare;
extern int updateSigma;
extern int nthreads;

extern double g_spStrFctAbserr;
extern double g_spStrFctRelerr;
extern double g_spStrFctMindx;
extern double g_spStrFctMaxdx;
extern double g_spsumnAbserr;
extern double g_spsumnRelerr;
extern double g_spsumnMindx;
extern double g_spsumnMaxdx;
extern double g_spHrAbserr;
extern double g_spHrRelerr;
extern double g_spHrMindx;
extern double g_spHrMaxdx;

extern double g_omega;
const char *Varname[] = {
  "Kc", "B", "Lr", "Mz", "D", "mosaic", "edisp", "bFWHM",
  "s", "bc2b", "wavelength", "pixelSize","qxzero","nindex",
  "T","bsig", (char *) NULL
};

char *dset;
int recordflag;
time_t elapsedMins = 0, elapsedSecs = 0, startTime;

double tcl_chi;
double ln[32],step[4]={1,1,1,1};
int boolx[32],xpin[4];
double ln1[32],ln0[32];
double sx[3];
int sxpin[4]={1,1,1,1};
double xisquare;
double NKparams[18];
SetUp setup;
Blt_Vector *xVP[4], *plotxVP, *plotyVP;


Tcl_Interp *NKinterp;
UpdateStruct *update; //Structure to deliver data to updatelinks_async
/*Handlers to safely perform updates of data in the Tcl script from threads
  other than the one handling the Tcl event loop
*/
Tcl_AsyncHandler updateHandler;
Tcl_AsyncHandler resetHandler;
Tcl_AsyncHandler timerHandler;

int colorMIN, colorMAX, colorBOT, colorSAT, colorFIL;

double STresult;

Tcl_HashTable datasetHT;
Tcl_HashTable modelccHT;
Tcl_HashTable paraHT;

Para g_ParaStruct;


GlobalData gData;
BoxVec boxvec;
YFImgBase *imgshown=NULL;

const char *YFtypeOptions[] = {"uint16", "int16", "float", (char *) NULL};

void textcolor(FILE* fp, int attr, int fg, int bg){
  fprintf(fp, "%c[%d;%d;%dm", 0x1B, attr, fg + 30, bg + 40);
}

void textcolor(ostream& out, int attr, int fg, int bg){
  char command[13];
  sprintf(command, "%c[%d;%d;%dm", 0x1B, attr, fg + 30, bg + 40);
  out<<command;
}


void checkImgUpdate(YFImgBase *img){
  /* check if the screen should be updated */
  if(img==imgshown) img->show(gData.tv_photoHd, &(gData.csetup));
}

int TV_BoxV(ClientData clientData, Tcl_Interp *interp,int objc, Tcl_Obj *const objv[]){
  /* handle the 'boxv' command */
  static const char *boxOptions[] = {
    "2nd", "1st", "del", "halo", "3rd", "coords", "set",
    (char *) NULL
  };
  enum options {
    BOX_2ND, BOX_1ST, BOX_DEL, BOX_HALO, BOX_3RD, BOX_COORDS, BOX_SET
  };
  int index;
  if (Tcl_GetIndexFromObj(interp, objv[1], boxOptions, "option", 0, &index) != TCL_OK) return TCL_ERROR;
  switch ((enum options) index) {
  case BOX_2ND:{
    double x, y;
    Tcl_GetDoubleFromObj(interp, objv[2], &x);
    Tcl_GetDoubleFromObj(interp, objv[3], &y);
    boxvec.step2((int)x,(int)y);
    break;
  }
  case BOX_1ST:{
    double x, y;
    Tcl_GetDoubleFromObj(interp, objv[2], &x);
    Tcl_GetDoubleFromObj(interp, objv[3], &y);
    if(objc<5) boxvec.step1((int)x,(int)y,NULL);
    else boxvec.step1(x,y, Tcl_GetString(objv[4]));
    break;
  }
  case BOX_3RD:{
    boxvec.step3();
    break;
  }
  case BOX_DEL:{
    boxvec.del(Tcl_GetString(objv[2]));
    break;
  }
  case BOX_HALO:{
    double x, y, x1, y1;
    Tcl_GetDoubleFromObj(interp, objv[2], &x);
    Tcl_GetDoubleFromObj(interp, objv[3], &x1);
    Tcl_GetDoubleFromObj(interp, objv[4], &y);
    Tcl_GetDoubleFromObj(interp, objv[5], &y1);
    boxvec.sethalo(fabs(x1-x), fabs(y-y1));
    //cout<<fabs(x-x1)<<' '<<fabs(y-y1)<<endl;
    break;
  }
  case BOX_COORDS:{
    Box* bptr;
    bptr = boxvec.boxPtr( Tcl_GetString( objv[2] ) );
    if(! bptr){
      return TCL_ERROR;
    }
    int x0 = bptr->xc[0];
    int x1 = bptr->xc[1];
    int y0 = bptr->yc[0];
    int y1 = bptr->yc[1];
    if(x0 > x1) tswap( x0 , x1);
    if(y0 > y1) tswap( y0 , y1);
    YFTcl_SetResult(interp, 128, "%d %d %d %d", x0, y0, x1, y1);
    YFTcl_EvalEx(interp, 128, "set boxx0 %d; set boxy0 %d; set boxx1 %d; set boxy1 %d", x0, y0, x1, y1);
    break;
  }
  case BOX_SET:{
    Box* bptr;
    bptr = boxvec.boxPtr( Tcl_GetString( objv[2] ) );
    int x , x1, y, y1;
/* old code: input is x0 y0 x1 y1
    Tcl_GetIntFromObj(interp, objv[3], &x);
    Tcl_GetIntFromObj(interp, objv[4], &y);
    Tcl_GetIntFromObj(interp, objv[5], &x1);
    Tcl_GetIntFromObj(interp, objv[6], &y1);
*/
// new code: input is x0 x1 y0 y1
    Tcl_GetIntFromObj(interp, objv[3], &x);
    Tcl_GetIntFromObj(interp, objv[4], &x1);
    Tcl_GetIntFromObj(interp, objv[5], &y);
    Tcl_GetIntFromObj(interp, objv[6], &y1);

    if(x > x1) tswap( x , x1);
    if(y > y1) tswap( y , y1);
    if(bptr == NULL)
      boxvec.newbox(x, y, x1, y1, Tcl_GetString( objv[2] ) );
    else{
      bptr->xc[0] = x;
      bptr->xc[1] = x1;
      bptr->yc[0] = y;
      bptr->yc[1] = y1;
      boxvec.updateBox(bptr);
    }
    break;
  }
  }
  return TCL_OK;
}

int tv_checkImgValidity_( int imgNum ){
  // for debugging
  return ( imgNum >= (int)gData.imgPtrVec.size() || imgNum < 0) ? TCL_ERROR : TCL_OK;
}

/*
  The image pointers are stored in a hash table and also in an array
*/
int tv_getExistingImgNum(YFImgBase* imgPtr){
  /* obtain the array index of the image pointer */
  for(unsigned int i=0; i<gData.imgPtrVec.size(); i++)
    if(gData.imgPtrVec[i]==imgPtr) return i;
  return -1;
}

int TV_ColorSetup(ClientData clientData, Tcl_Interp *interp,int objc, Tcl_Obj *const objv[]){
  /* handles the 'colorsetup' command */
  const char *Options[] = {
    "minmax", "mode", "bottom", "saturation", "filter", "style", "disp", "update",
    (char *) NULL
  };
  enum options {
    CLR_MINMAX, CLR_MOD, CLR_BOT, CLR_SAT, CLR_FIL, CLR_STY, CLR_DISP, CLR_UPDATE
  };
  int index;
  if (Tcl_GetIndexFromObj(interp, objv[1], Options, "option", 0, &index) != TCL_OK) return TCL_ERROR;

  switch ((enum options) index) {
  case CLR_MINMAX:{
    double tmp;
    if(objc>2) Tcl_GetDoubleFromObj(interp, objv[2], &tmp);
    gData.csetup.min = tmp;
    if(objc>3) Tcl_GetDoubleFromObj(interp, objv[3], &tmp);
    gData.csetup.max = tmp;
    break;}
  case CLR_SAT:{
    double sat;
    Tcl_GetDoubleFromObj(interp, objv[2], &sat);
    gData.csetup.saturation = sat;
    break;}
  case CLR_BOT:{
    double bot;
    Tcl_GetDoubleFromObj(interp, objv[2], &bot);
    gData.csetup.bottom = bot;
    break;}
  case CLR_FIL:{
    double fil;
    Tcl_GetDoubleFromObj(interp, objv[2], &fil);
    gData.csetup.filter = fil;
    break;}
  case CLR_MOD:{
    const char *Options[] = { "linear", "log", (char *) NULL };
    enum options { LIN, LOG };
    if (Tcl_GetIndexFromObj(interp, objv[2], Options, "option", 0, &gData.csetup.mode) != TCL_OK) return TCL_ERROR;
    break;}
  case CLR_STY:{
    const char *Options[] = { "gray1", "gray2", (char *) NULL };
    int style;
    if (Tcl_GetIndexFromObj(interp, objv[2], Options, "option", 0, &style) != TCL_OK) return TCL_ERROR;
    gData.csetup.buildColorTable( (ColorStyle)style );
    break;}
  case CLR_DISP:{
    cout<<"min: "<<gData.csetup.min<<endl;
    cout<<"max: "<<gData.csetup.max<<endl;
    cout<<"sat: "<<gData.csetup.saturation<<endl;
    cout<<"bot: "<<gData.csetup.bottom<<endl;
    cout<<"fil: "<<gData.csetup.filter<<endl;
    cout<<"mod: "<<gData.csetup.mode<<endl;
    break;}
  case CLR_UPDATE:{
    char cmd[15];
    colorSAT=(int)gData.csetup.saturation;
    sprintf(cmd,"colorSAT");
    Tcl_UpdateLinkedVar(interp, cmd);
    colorMAX=(int)gData.csetup.max;
    sprintf(cmd,"colorMAX");
    Tcl_UpdateLinkedVar(interp, cmd);
    colorMIN=(int)gData.csetup.min;
    sprintf(cmd,"colorMIN");
    Tcl_UpdateLinkedVar(interp, cmd);
    colorBOT=(int)gData.csetup.bottom;
    sprintf(cmd,"colorBOT");
    Tcl_UpdateLinkedVar(interp, cmd);
    colorFIL=(int)gData.csetup.filter;
    sprintf(cmd,"colorFIL");
    Tcl_UpdateLinkedVar(interp, cmd);
    break;}
  }
  return TCL_OK;
}

int csix, csiy;

int TV_PlotCS(ClientData clientData, Tcl_Interp *interp,int objc, Tcl_Obj *const objv[]){
  /* plot cross-sections */
  double x,y;

  if(objc>2){
    Tcl_GetDoubleFromObj(interp, objv[1], &x);
    Tcl_GetDoubleFromObj(interp, objv[2], &y);
    csix = (int)x; csiy = (int)y;

    YFTcl_EvalEx(interp, 128, "set tvar_xpos %d; set tvar_ypos %d", csix, csiy );
    if(imgshown!=NULL){
      YFTcl_EvalEx(interp, 128, "set tvar_inte %g", imgshown->getvalsafe(csix, csiy) );
    }
  }
  for(unsigned int i=0; i<gData.imgPtrVec.size(); i++)
    gData.imgPtrVec[i]->plotCrossSection(csix, csiy, gData.plotopt);
  return TCL_OK;
}

int TV_PlotHide(ClientData clientData, Tcl_Interp *interp,int objc, Tcl_Obj *const objv[]){
  int i, opt;
  Tcl_GetIntFromObj(interp, objv[1], &i);
  Tcl_GetIntFromObj(interp, objv[2], &opt);
  gData.imgPtrVec[i]->plotHide(opt);
  return TCL_OK;
}

int NK_qzrexport(ClientData clientData, Tcl_Interp *interp,int objc, Tcl_Obj *const objv[]){
  FILE *fp;
  fp=fopen(Tcl_GetString(objv[1]),"w");
  if(fp==NULL) return TCL_OK;
  int size, maxsize=0;
  double x,y;
  for (unsigned int i=0; i<gData.imgPtrVec.size(); i++) {
    size=gData.imgPtrVec[i]->getSize();
    if (size>maxsize) maxsize=size;
//    cout<<size<<endl;
  }
  for (int j=0; j<maxsize; j++) {
    for (unsigned int i=0; i<gData.imgPtrVec.size(); i++) {
      size=gData.imgPtrVec[i]->getSize();
      if (size==0) continue;
      if (size>j) {
        gData.imgPtrVec[i]->getValue(j,&x,&y);
        fprintf(fp,"%g\t%g\t", x, y);
      } else fprintf(fp,"\t\t");
    }
    fprintf(fp,"\n");
  }
  fclose(fp);
  return TCL_OK;
}

int readDataFromASCII(ClientData clientData, Tcl_Interp *interp,int objc, 
                      Tcl_Obj *const objv[])
{
  YFImgBase* img;
  const char *fileRootName = Tcl_GetString(objv[1]);
  const char *absFilePath = Tcl_GetString(objv[2]);
    
	TVImg<float> *newimg = new TVImg<float>;
	img = newimg;
	
	img->openImgASCII(absFilePath);
	img->setInterp(interp);

  int newflag = 1;
  Tcl_HashEntry *entry;
  entry = Tcl_CreateHashEntry(&gData.imgNT, fileRootName, &newflag);
  entry->clientData = img;
  img->setname(entry->key.string);
  gData.imgPtrVec.push_back(img);
  
  return TCL_OK;
}

int TV_ImgCmd(ClientData clientData, Tcl_Interp *interp,int objc, Tcl_Obj *const objv[]){
  /* handles 'img' command */
  const char *imgOptions[] = {
    "load", "delete", "report", "flip",
    "transpose", "rotate", "coords",
    "show", "convert", "expr", "shifttif","shiftfit",
    "export", "configure", "names", "toggleplot", "fitlbg",
    "medianfilter", "convfilter", "stat", "affine", "rename","fftsmooth", "address",
    "bgadjust","scadjust",
    (char *) NULL
  };

  enum options {
    IMG_LOAD, IMG_DELETE, IMG_REPORT, IMG_FLIP,
    IMG_TRANSPOSE, IMG_ROTATE, IMG_COORDS,
    IMG_SHOW, IMG_CONVERT, IMG_EXPR, IMG_SHIFTTIF, IMG_SHIFTFIT,
    IMG_EXPORT, IMG_CONFIGURE, IMG_NAMES, IMG_TOGGLEPLOT, IMG_FITLBG,
    IMG_MEDIANF, IMG_CONVF,  IMG_STAT, IMG_AFFINE, IMG_RENAME, IMG_FFTSMOOTH, IMG_ADDRESS,
    IMG_BGADJUST, IMG_SCADJUST
  };
  const char *directOpts[]={
    "x", "y", NULL};
  enum direction { X, Y };
  int index;
  YFImgBase* img;
  Tcl_HashEntry *entryPtr;

  if (Tcl_GetIndexFromObj(interp, objv[1], imgOptions, "option", 0,
			  &index) != TCL_OK) return TCL_ERROR;

  switch((enum options) index) {
		case IMG_NAMES: {
		  Tcl_HashSearch search;
		  Tcl_Obj* r = Tcl_NewStringObj(NULL,0);
		  entryPtr=Tcl_FirstHashEntry(&gData.imgNT, &search);
		  while(entryPtr!=NULL){
		    Tcl_AppendStringsToObj(r, entryPtr->key.string, NULL);
		    Tcl_AppendStringsToObj(r, " ", NULL);
		    entryPtr=Tcl_NextHashEntry(&search);
		  }
		  Tcl_SetObjResult(interp, r);
		  return TCL_OK;
		}
		default:{ }
  }

  entryPtr=_getExistingPtr(interp, &gData.imgNT, objv[2], (void**)&img);
  if( img == NULL && (enum options) index != IMG_LOAD) {
    cout<<"hash table error"<<endl;
    cout<<Tcl_HashStats(&gData.imgNT)<<endl;
    return TCL_ERROR;
  }
  switch ((enum options) index) {
  case IMG_LOAD:{
    TIFF *tif;
    uint16 bitspersample;
    int newflag=0;
    int typein=-1;
    int typeout=-1;

    if((tif=TIFFOpen(Tcl_GetString(objv[3]), "r"))==NULL) return TCL_ERROR;
    if( objc > 4 )
      if (Tcl_GetIndexFromObj(interp, objv[4], YFtypeOptions, "typeOpt", 0,
			      &typeout) != TCL_OK) return TCL_ERROR;
    if( objc > 5 )
      if (Tcl_GetIndexFromObj(interp, objv[5], YFtypeOptions, "typeOpt", 0,
			      &typein) != TCL_OK) return TCL_ERROR;

    TIFFGetFieldDefaulted(tif, TIFFTAG_BITSPERSAMPLE, &bitspersample);

	  //Reminder: enum DataType {yfUINT16, yfINT16, yfFLOAT, yfINT32};
	  //FLICAM is the ccd at G1 station, which has the dynamic range
    if(bitspersample == 16) {
    	cout << "FLICAM = " << FLICAM << endl;
    	if (FLICAM != 0 && FLICAM != 1) return TCL_ERROR;
    	if (FLICAM) {
    		if (typeout < 0) {
    			typeout = yfFLOAT;
    			cout << "typeout = yfFLOAT" << endl;
    		}
    		if (typein < 0) {
    			typein = yfUINT16;
    			cout << "typein = yfUINT16" << endl;
    		}
    	} else {
    		if(typeout<0) {
    			typeout = yfINT16;
    			cout << "typeout = yfINT16" << endl;
    		}
      		if(typein<0) {
      			typein = yfUINT16;
      			cout << "typein = yfUINT16" << endl;
      		}
      		if(typein != yfINT16 && typein != yfUINT16) return TCL_ERROR;
    	}
    } else if (bitspersample == 32) {
      if(typeout<0) {
      	typeout = yfFLOAT;
      	cout << "typeout = yfFLOAT" << endl;
      }
      if(typein<0) {
      	typein = yfFLOAT;
      	cout << "typein = yfFLOAT" << endl;
      }
      if(typein != yfFLOAT) return TCL_ERROR;
    }
    if( img == NULL ) {
      newflag=1;
      switch(typeout) {
				case yfUINT16:{
					TVImg<unsigned short> *newimg = new TVImg<unsigned short>;
					img = newimg;
					break;}
				case yfINT16:{
					TVImg<short> *newimg = new TVImg<short>;
					img = newimg;
					break;}
				case yfINT32:{
					TVImg<int> *newimg = new TVImg<int>;
					img = newimg;
					break;}
				case yfFLOAT:{
					TVImg<float> *newimg = new TVImg<float>;
					img = newimg;
					break;}
				default:{
					cout<<"I can only handle 16(short) or 32(float) bits data"<<endl;
					return TCL_ERROR;
					break;}
      }
    }
    if( (img->openImg(tif, typein, typeout) ) != TCL_OK ) {
      cout<<"error in loading file"<<endl;
      return TCL_ERROR;
    } else {
      img->setInterp(interp);
    }
    if(newflag){
      Tcl_HashEntry *entry=Tcl_CreateHashEntry(&gData.imgNT, Tcl_GetString(objv[2]), &newflag);
      entry->clientData=img;
      img->setname(entry->key.string);
      gData.imgPtrVec.push_back(img);
    }
    break;}
  case IMG_RENAME:{
    int newflag;
    Tcl_DeleteHashEntry(entryPtr);
    Tcl_HashEntry *entry=Tcl_CreateHashEntry(&gData.imgNT, Tcl_GetString(objv[3]), &newflag);
    entry->clientData=img;
    img->setname(entry->key.string);
    break;}
  case IMG_ROTATE:{
    double rad;
    double cx, cy;
    int w = img->imgwidth();
    int h = img->imgheight();

    if(objc>5) {
      Tcl_GetDoubleFromObj(interp, objv[4], &cx);
      Tcl_GetDoubleFromObj(interp, objv[5], &cy);
    } else {
      cx=w/2.0;
      cy=h/2.0;
    }

    if(strcmp(Tcl_GetString(objv[3]),"inf")==0){
      img->nc_rotate90(1, w, h, cx, cy);
    } else if (strcmp(Tcl_GetString(objv[3]),"-inf")==0){
      img->nc_rotate90(-1, w, h, cx, cy);
    } else if( Tcl_GetDoubleFromObj(interp, objv[3], &rad) == TCL_OK ){
      img->rotate(rad, w, h, cx, cy, cx, cy);
    } else if( Tcl_GetDoubleFromObj(interp, objv[3], &rad) != TCL_OK ) return TCL_ERROR;
    break;}

  case IMG_AFFINE:{
    double rad, scale, tx, ty;
    if( objc < 7) return TCL_ERROR;
    if( Tcl_GetDoubleFromObj(interp, objv[3], &rad) != TCL_OK ) return TCL_ERROR;
    if( Tcl_GetDoubleFromObj(interp, objv[4], &scale) != TCL_OK ) return TCL_ERROR;
    if( Tcl_GetDoubleFromObj(interp, objv[5], &tx) != TCL_OK ) return TCL_ERROR;
    if( Tcl_GetDoubleFromObj(interp, objv[6], &ty) != TCL_OK ) return TCL_ERROR;

    img->affine(rad, scale,  tx, ty);

    break;}

  case IMG_REPORT:
    img->report(stdout);
    break;

  case IMG_COORDS:{
    YFRect<int> rect = img->coords();
    char *str=Tcl_Alloc(64);
    sprintf(str, "%d %d %d %d", rect.llx(), rect.lly(), rect.urx(), rect.ury());
    Tcl_SetResult(interp, str, TCL_DYNAMIC);
    break;}

  case IMG_ADDRESS:{
    char *str=Tcl_Alloc(64);
    sprintf(str, "%ld",  (long unsigned int)img);
    Tcl_SetResult(interp, str, TCL_DYNAMIC);
    break;}

  case IMG_FLIP:{
    int direct;
    if (Tcl_GetIndexFromObj(interp, objv[3], directOpts, "directOpt", 0,
			    &direct) != TCL_OK) return TCL_ERROR;
    switch( (enum direction) direct ){
    case X:
      img->flipx();break;
    case Y:
      img->flipy();break;
    }
    break;}

  case IMG_SHOW:
    imgshown = img;
    //cout<<"imgshown: "<<imgshown<<endl;
    checkImgUpdate(img); break;

  case IMG_TRANSPOSE:
    img->transpose();
    break;

  case IMG_DELETE:{
    for(int i=2; i<objc; i++){
      entryPtr=_getExistingPtr(interp, &gData.imgNT, objv[i], (void**)&img);
      if( img == NULL ) continue;
      img->cleanup(interp);
      Tcl_DeleteHashEntry(entryPtr);
      gData.imgPtrVec.erase(gData.imgPtrVec.begin()+tv_getExistingImgNum(img));
      if(img == imgshown) imgshown=NULL;
      delete img;
    }
    break;}

  case IMG_SHIFTTIF:{
  /*
  This command shifts each image pixel by x and y inside the area of the original picture
  */
    YFDuple<int> rv;
    Tcl_GetIntFromObj(interp, objv[3], &rv.x);
    Tcl_GetIntFromObj(interp, objv[4], &rv.y);
    img->shift(rv);
    //img->shiftposition(rv);
    break;}

  case IMG_SHIFTFIT:{
  /*
  This command shifts the whole area of the image, but it leaves the pixel positions the same;
  so the left lower corner will be still [0;0] when it is going to be fit
  */
    YFDuple<int> rv;
    Tcl_GetIntFromObj(interp, objv[3], &rv.x);
    Tcl_GetIntFromObj(interp, objv[4], &rv.y);
    //img->shift(rv);
    img->shiftposition(rv);
    break;}

  case IMG_EXPORT:{
    img->exportTiff(Tcl_GetString(objv[3]));
    break;}
  case IMG_CONFIGURE:{
    //    img->configurePlot(objv+3, objc-3);
    break;}
  case IMG_TOGGLEPLOT:{
    int togopt=2;
    if(objc>3) Tcl_GetIntFromObj(interp, objv[3], &togopt);
    img->plotHide(togopt);
    break;}
  case IMG_FITLBG:{
    Box *b1, *b2;
    int order=2;
    int opt = 4;
    int step = 1;
    if(objc<5) return TCL_ERROR;
    b1 = boxvec.boxPtr( Tcl_GetString(objv[3]) );
    b2 = boxvec.boxPtr( Tcl_GetString(objv[4]) );
    if(b1==NULL || b2==NULL) return TCL_ERROR;
    if(objc>5) Tcl_GetIntFromObj(interp, objv[5], &order);
    if(objc>6) Tcl_GetIntFromObj(interp, objv[6], &opt);
    if(objc>7) Tcl_GetIntFromObj(interp, objv[7], &step);
    if(order > 20) return TCL_OK;
    img->fit2d( b1, b2, order, opt, step);
    break;}
  case IMG_STAT:{
    Box *b;
    char *res = (char*)malloc(128*sizeof(char));
    double report[5];
    b = boxvec.boxPtr( Tcl_GetString(objv[3]) );
    if(!b){
      free(res);
      return TCL_ERROR;
    }
    img->stat( b, report);
    sprintf(res, "%g %g %g %g %g", report[0], report[1], report[2], report[3], report[4]);
    STresult=report[4];
    cout<<res<<endl;
    char cmd[10];
    sprintf(cmd,"STresult");
    Tcl_UpdateLinkedVar(interp, cmd);
    break;}
  case IMG_MEDIANF:{
    break;}
  case IMG_CONVF:{
    break;}
  case IMG_FFTSMOOTH:{
    break;}
  case IMG_EXPR:{
    char *op = Tcl_GetString(objv[3]);
    char *next = Tcl_GetString(objv[4]);
    char *tmp;
    double scalar;

    scalar = strtod(next, &tmp);
    if( tmp != next ){
      cout<<"scalar = "<<scalar<<endl;
      if( ! strcmp (  op , "+=" ) ){
	(*img) += scalar;
      } else if( ! strcmp ( op , "*=" ) ){
	(*img) *= scalar;
      } else if( ! strcmp ( op , "-=" ) ){
	(*img) += (-scalar);
      } else if( ! strcmp ( op , "/=" ) ){
	(*img) *= (1.0/scalar);
      } else {
	cout<<"unknown operations"<<endl;
      }
    } else {
      YFImgBase * e;
      cout<<"two img operations"<<endl;
      entryPtr=_getExistingPtr(interp, &gData.imgNT, objv[4], (void**)&e);
      if(entryPtr == NULL) return TCL_ERROR;
      if( ! strcmp (  op , "+=" ) ){
	(*img) += (*e);
      } else if( ! strcmp ( op , "*=" ) ){
	(*img) *= (*e);
      } else if( ! strcmp ( op , "-=" ) ){
	(*img) -= (*e);
      } else if( ! strcmp ( op , "/=" ) ){
	(*img) /= (*e);
      } else {
	cout<<"unknown operations"<<endl;
      }
    }
    break;}
  case IMG_BGADJUST:{
  //subtratcs cz parameter from fitted image
    extern Tcl_HashTable datasetHT;
    Data *dp;
    entryPtr=_getExistingPtr(interp, &datasetHT, objv[3], (void**)&dp);
    img->bgadjust(dp);
    break;}
  case IMG_SCADJUST:{
  //multiplies image with scalling factors from external file
  /*
  use to simulate data: 1. structure factor from toad (use noscale=1)
  			2. scalling factors from df
			3. multiply in toad (img scadjust "imgname" "ext.file")
  */
    int p;
    double sc[2049],tmpsc;
    for (p=0; p<2049;p++) sc[p]=1;
    FILE *frm;
    if ((frm=fopen(Tcl_GetString(objv[3]), "r"))==NULL) return TCL_ERROR;
    while (fscanf(frm,"%d %lf %*[^\n]\n" ,&p,&tmpsc) == 2) {sc[p]=tmpsc;}
    fclose(frm);
    img->scadjust(sc);
    break;}
  default:{ }
  }
  return TCL_OK;
}


/****************************************************************************************
dataset op name ....
****************************************************************************************/

int YF_dataset(ClientData clientData, Tcl_Interp *interp,
               int objc, Tcl_Obj *const objv[]) {
  const char *Options[] = {
    "load", "delete", "info", "filter",
    (char *) NULL
  };
  enum options {
    DATA_LOAD, DATA_DEL, DATA_INFO, DATA_FILTER
  };
  int index;
  Tcl_HashEntry *entryPtr;
  Data *dp;


  if (Tcl_GetIndexFromObj(interp, objv[1], Options, "option", 0, &index) != TCL_OK) return TCL_ERROR;

  entryPtr=_getExistingPtr(interp, &datasetHT, objv[2], (void**)&dp);
  if( dp == NULL && (enum options) index != DATA_LOAD){
    cout<<"hash table error"<<endl;
    cout<<Tcl_HashStats(&datasetHT)<<endl;
    return TCL_ERROR;
  }

  switch ((enum options) index) {
  case DATA_LOAD:{
    /*
      dataset new name filename basecnt qrstart qrend qzstart qzend .. .. .. ..
      0       1   2    3        4         5       6      7      8
    */
    int newflag = 0;
    Instruction *ins = new Instruction;
    if( dp == NULL){
      newflag = 1;
      dp = new Data;
    }
    int temp;
    double dtemp;
    int err=0;
    double range=0.0;


    Tcl_GetDoubleFromObj(interp, objv[4], &dtemp);
    ins->basecnt = dtemp ;

    Tcl_GetIntFromObj(interp, objv[5], &temp);
    ins->qrstart = temp ;
    range = temp;
    Tcl_GetIntFromObj(interp, objv[6], &temp);
    ins->qrend   = temp ;
    range -= temp;
    if ( range > 0 ) {err=1;}
    for(int i = 7; i < objc; i+=2){
      range = 0.0;
      Tcl_GetIntFromObj(interp, objv[i], &temp);
      ins->qzstart.push_back( temp );
      range = temp;
      Tcl_GetIntFromObj(interp, objv[i+1], &temp);
      ins->qzend.push_back( temp );
      range -= temp;
      if ( range > 0 ) {err=1;}
    }

    if ( err == 1 ){
      cout<<"wrong order of arguments"<<endl;
      Tcl_EvalEx(interp,"set startflag 0",-1,TCL_EVAL_GLOBAL);
      return TCL_ERROR;
    }

    YFImgBase *img;
    entryPtr=_getExistingPtr(interp, &gData.imgNT, objv[3], (void**)&img);

    dp->readin( img, ins );
    dp->flag = false;

    TVImg<float> *fitimg;
    entryPtr=_getExistingPtr(interp, &gData.imgNT, (char *) "FIT", (void**)&fitimg);
    if(fitimg==NULL){
      int newflag;
      fitimg = new TVImg<float>;
      Tcl_HashEntry *entry = Tcl_CreateHashEntry(&gData.imgNT, "FIT", &newflag);
      entry->clientData = fitimg;
      fitimg->setname(entry->key.string);
      gData.imgPtrVec.push_back(fitimg);
      fitimg->setInterp(interp);
    }
    dp->setimg(fitimg);

    if(newflag){
      Tcl_HashEntry *entry = Tcl_CreateHashEntry(&datasetHT, Tcl_GetString(objv[2]), &newflag);
      entry->clientData = dp;
    }
    break;}
  case DATA_DEL:{
    Tcl_DeleteHashEntry(entryPtr);
    delete dp;
    break;}
  case DATA_INFO:{
    dp->info();
    break;}

  case DATA_FILTER:{
    YFImgBase *img;
    entryPtr=_getExistingPtr(interp, &gData.imgNT, objv[3], (void**)&img);

    TVImg<float> *filterimg;
    entryPtr=_getExistingPtr(interp, &gData.imgNT, (char *) "FILTER", (void**)&filterimg);
    if(filterimg==NULL){
      int newflag;
      filterimg = new TVImg<float>;
      Tcl_HashEntry *entry = Tcl_CreateHashEntry(&gData.imgNT, "FILTER", &newflag);
      entry->clientData = filterimg;
      filterimg->setname(entry->key.string);
      gData.imgPtrVec.push_back(filterimg);
      filterimg->setInterp(interp);
    }
    dp->filter(filterimg);
    break;}
  }
  return TCL_OK;
}

/******************************************************************************
modelcc op name
******************************************************************************/
int YF_modelcc(ClientData clientData, Tcl_Interp *interp,int objc, Tcl_Obj *const objv[])
{
  const char *Options[] = {
    "new", "delete","setutabpar","wutab","tutab", "retol",
    (char *) NULL
  };
  enum options {
    MODCC_NEW, MODCC_DEL, MODCC_SETUTABPAR, MODCC_WUTAB, MODCC_TUTAB, MODCC_RETOL
  };
  int index;
  Tcl_HashEntry *entryPtr;
  ModelCalculator *mc;

  if (Tcl_GetIndexFromObj(interp, objv[1], Options, "option", 0, &index) != TCL_OK) return TCL_ERROR;

  entryPtr=_getExistingPtr(interp, &modelccHT, objv[2], (void**)&mc);
  if( mc == NULL && (enum options) index != MODCC_NEW){
    cout<<"hash table error"<<endl;
    cout<<Tcl_HashStats(&modelccHT)<<endl;
    return TCL_ERROR;
  }

  switch ((enum options) index) {
  case MODCC_NEW: {
    int newflag=1;
    if( mc != NULL) break;
    mc = new ModelCalculator;

    mc->setModelParameter(0.01, "edisp");
    mc->read_in_utable((char *)"dat/utab_nfit12.15.dat");

    Tcl_HashEntry *entry = Tcl_CreateHashEntry(&modelccHT, Tcl_GetString(objv[2]), &newflag);
    entry->clientData = mc;
    } break;
  case MODCC_DEL: {
    Tcl_DeleteHashEntry(entryPtr);
    delete mc;
    } break;
  case MODCC_SETUTABPAR:
    break;
  case MODCC_WUTAB: {
    Utable utab;
    utab.writeUtableFile("new_utab.dat", 600, 2000);
    break; }
  case MODCC_RETOL:
    mc->resetTOL();
    break;
  case MODCC_TUTAB:
    break;
  default:
    break;
  }
  return TCL_OK;
}


/******************************************************************************
updateParaObject
******************************************************************************/
int updateParaObject(ClientData clientData, Tcl_Interp *interp,
                     int objc, Tcl_Obj *const objv[])
{
  //int index;
  double value;
  if (objc != 3) return TCL_ERROR;
  //if (Tcl_GetIndexFromObj(interp, objv[1], Varname, "varname", 0, &index) != TCL_OK)
  //  return TCL_ERROR;
  string paraName(Tcl_GetString(objv[1]));
  Tcl_GetDoubleFromObj(interp, objv[2], &value);
  g_ParaStruct.setValue(value, paraName);
  return TCL_OK;
}

/******************************************************************************
setNFIT
******************************************************************************/
int setNFIT(ClientData clientData, Tcl_Interp *interp, int objc, Tcl_Obj *const objv[])
{
  int numFreeParams = objc - 1;
  g_ParaStruct.setNfit(interp, objv + 1, numFreeParams);
  return TCL_OK;
}

/****************************************************************************************
para op name
****************************************************************************************/
/*
int YF_paraset(ClientData clientData, Tcl_Interp *interp,int objc, Tcl_Obj *const objv[])
{
  const char *Options[] = {
    "set", "update", "delete", "nfit", "disp", "dfit", "qz", "qr",
    (char *) NULL
  };
  enum options {
    PARA_SET, PARA_UPDATE, PARA_DEL, PARA_NFIT, PARA_DISP, PARA_DFIT, PARA_QZ, PARA_QR
  };
  int index;
  Tcl_HashEntry *entryPtr;
  Para *p;

  if( objc < 3) return TCL_ERROR;
  if(Tcl_GetIndexFromObj(interp, objv[1], Options, "option", 0, &index) != TCL_OK) return TCL_ERROR;

  entryPtr=_getExistingPtr(interp, &paraHT, objv[2], (void**)&p);
  if( p == NULL && (enum options) index != PARA_SET){
    cout<<(char)0x1B<<"[31;1mhash table error[0m"<<endl;
    return TCL_ERROR;
  }

  switch ((enum options) index) {
  case PARA_SET:
    int newflag=0;
    if( p == NULL) {
      p = new Para;
      newflag = 1;
    }

    double *paraptr = (double*)p;

    if( objc < 4 ) return TCL_ERROR;
    if( isalpha( *Tcl_GetString(objv[3]) ) ) {
      double a;
      if( objc < 5 ) return TCL_ERROR;
      if (Tcl_GetIndexFromObj(interp, objv[3], Varname, "varname", 0, &index) != TCL_OK) return TCL_ERROR;
      if (Tcl_GetDoubleFromObj(interp, objv[4], &a) != TCL_OK) return TCL_ERROR;
      if (index==18) {
	      a=1.18/YFPI*a*paraptr[9]*paraptr[13]/paraptr[14];
	      index=8;
      }
      paraptr[index] = a;
    } else {
      for(int i=3; i<objc; i++) Tcl_GetDoubleFromObj(interp, objv[i], paraptr+i-3 );
    }
    if(newflag == 1) {
      Tcl_HashEntry *entry = Tcl_CreateHashEntry(&paraHT, Tcl_GetString(objv[2]), &newflag);
      p->nfit = 0;
      entry->clientData = p;
    }
    break;
  case PARA_UPDATE:
    char cmd[20];
    for(int i=0; i<18; i++) {
      NKparams[i]=((double*)p)[i];
      sprintf(cmd,"NKparams(%d)",i);
      Tcl_UpdateLinkedVar(interp, cmd);
    }
    int i=0;
    NKparams[i]=int(1e17*((double*)p)[i]+0.5)*1e-17;
    sprintf(cmd,"NKparams(%d)",i);
    Tcl_UpdateLinkedVar(interp, cmd);
    i=1;
    NKparams[i]=int(1e-9*((double*)p)[i]+0.5)*1e9;
    sprintf(cmd,"NKparams(%d)",i);
    Tcl_UpdateLinkedVar(interp, cmd);
    break;
  case PARA_DEL:
    Tcl_DeleteHashEntry(entryPtr);
    delete p;
    break;

  case PARA_NFIT:
    //objv+3 because the first arguments are "paraset nfit p"
    //objc-3 is equal to the number of free parameters
    p->setNfit(interp, objv+3, objc-3);
    break;

  case PARA_DISP:{
    for(int i=0; i<18; i++)
      cout<<i<<' '<<Varname[i]<<"\t:\t"<<((double*)p)[i]<<endl;
    cout<<"\nnfit = "<<p->nfit<<" {";
    for(int i=0; i<p->nfit; i++){
      cout<<p->idx[i]<<' ';
    }
    cout<<'}'<<endl;
    break;}

  case PARA_DFIT:{
    int num = (objc-2)/2;
    double *pp = new double[num];
    int *pn = new int[num];
    for(int i=0; i<num; i++){
      Tcl_GetIntFromObj(interp, objv[3+i*2], pn+i);
      Tcl_GetDoubleFromObj(interp, objv[3+i*2+1], pp+i);
    }
    p->dfit(num, pn, pp);
    break;}
  case PARA_QZ:{
    int pix;
    Tcl_GetIntFromObj(interp, objv[3], &pix);
    cout<<p->setup.getqz(pix)<<endl;
    break;}
  case PARA_QR:{
    int pix;
    Tcl_GetIntFromObj(interp, objv[3], &pix);
    cout<<p->setup.getqr(pix)<<endl;
    break;}
  }
  return TCL_OK;
}
*/

/* multi-threading utilities */
struct NfitThreadPars {
  int iter;
  Data *dp;
  ModelCalculator *mc;
  Para *p;
  FILE *fp;
};

/* error - wrapper for perror used for bad syscalls */
void error(char *msg) {
  perror(msg);
  exit(1);
}

/* initialize semaphore sem to value */
/* pshared=0 if thread, pshared=1 if process */
void Sem_init(sem_t *sem, int pshared, unsigned int value){
  if (sem_init(sem, pshared, value) < 0) error((char *) "Sem_init");
}

/* P operation on semaphore sem */
void P(sem_t *sem){
  if (sem_wait(sem)) error((char *) "P");
}

/* V operation on semaphore sem */
void V(sem_t *sem){
  if (sem_post(sem)) error((char *) "V");
}

/*
 * Pthread_create - wrapper for pthread_create
 */
void Pthread_create(pthread_t *thread, pthread_attr_t *attr,
		    void *(*start_routine)(void *), void *arg){
  if(pthread_create(thread, attr, start_routine, arg)) error((char *) "ERROR creating thread");
}
/*
 * Pthread_detach - wrapper for pthread_detach
 */
void Pthread_detach(pthread_t tid){
  if(pthread_detach(tid)) error((char *) "ERROR detaching thread");
}

/*This alters data shared with the Tcl script. For thread safety, this method
  should be invoked through the thread-safe resetFlags function.
*/
int resetFlags_async(ClientData clientData, Tcl_Interp *interp, int code)
{
  interp = NKinterp;
  if(recordflag) Tcl_EvalEx(interp,"recordfit finished",-1,TCL_EVAL_GLOBAL);
  Tcl_EvalEx(interp,"set stopflag 0",-1,TCL_EVAL_GLOBAL);
  Tcl_EvalEx(interp,"set startflag 0",-1,TCL_EVAL_GLOBAL);
  return code;
}

void resetFlags(){
  //Wait until it is safe to do an asynchronous modification.
  Tcl_AsyncMark(resetHandler);
}

/*Update variables tracking the minutes and seconds elapsed since the beginning
  of the most recent fit.
*/
void updateTimer(){
  time_t diff = time(NULL) - startTime;
  elapsedMins = diff / 60;
  elapsedSecs = diff % 60;
  Tcl_AsyncMark(timerHandler);
}

/*This function updates the minutes and seconds since the beginning of the most
  recent fit within the Tcl script. As this function is expected to be used
  a thread seperate from the one in which the Tcl event loop is running, the
  thread-safe updateTimer function should be used to invoke this function.
*/
int updateTimer_async(ClientData clientData, Tcl_Interp *interp, int code){
  Tcl_UpdateLinkedVar(NKinterp, "elapsedMins");
  Tcl_UpdateLinkedVar(NKinterp, "elapsedSecs");
  return code;
}

/*As long as a fit is running, continue to update the minutes and seconds
  elapsed since the start of the most recent fitting run.
*/
void* runTimer(void* arg){
  while( startflag ){
    sleep(1);
    updateTimer();
  }
  return NULL;
}

void updatelinks(double chisq, char *chain);

/****************************************************************************
Either pthread_join(3) or pthread_detach() should be called for each
thread that an application creates, so that system resources for the
thread can be released.  (But note that the resources of all threads
are freed when the process terminates.)

Pthread_create(&tid, NULL, thread_lmdif, (void *)s);
binds a newly created thread called tid with
void *thread_lmdif(void *s) routine. The argument, void *s, or
(NfitTrheadPars *)s contains all the info regarding a fit.

FuncLmdif object essentially atcs as a wrapper between ModelCalculator object
and the fitting routine from MINPACK. It keeps a track of free parameters
for a current fit and updates the model parameters for each iteration.
*****************************************************************************/
void *thread_lmdif(void *s)
{
  Pthread_detach(pthread_self());
  NfitThreadPars *ntp = (NfitThreadPars *)s;
  FuncLmdif *func = new FuncLmdif;
  func->setData(ntp->dp);
  func->setMC(ntp->mc);
  func->setPara(ntp->p);
  Lmdif *lmdif = new  Lmdif;
  char *bestChain;
  double bestChisq;
  pthread_t timerThread; // for the Elapsed Time window

  cout << "iter: " << ntp->iter << endl;
  cout << "nthreads: " << nthreads << endl;
  /*
  void Lmdif::setup(LmdifFunc f, int _m, int _n, int _maxfev, double _ftol,
         double _gtol, double _xtol, double _factor, double *_epsfcn, void *_client)

  FuncLmdif::WrapperFunclmdif: (I_e-I_c) vector wrapper
  func->npoint(): total number of data points
  ntp->p->nfit: number of free parameters
  ntp-iter: maximum iterations
  ftol
  gtol
  xtol
  LAMBDA_LM: lambda of L-M algorithm
  ntp->p->epsfcn: precision
  func: (I_e-I_c) vector
  */
  const double FTOL = 0;
  const double GTOL = 0;
  const double XTOL = 0;
  const double LAMBDA_LM = 0.1;
  lmdif->setup(FuncLmdif::WrapperFunclmdif, func->npoint(), ntp->p->nfit, ntp->iter,
               FTOL, GTOL, XTOL, LAMBDA_LM, ntp->p->epsfcn, func);
  lmdif->init(ntp->p->xp, ntp->p->nfit);
  startTime = time(NULL);
  pthread_create(&timerThread, NULL, runTimer, NULL);

  // Begins a fit with ntp->iter number of iterations
  lmdif->fit();

  // covar_cminpack() results in segfault if total iteration is less than 2
  // derivative gets calculated for iter > 1
	if (ntp->iter > 1) lmdif->covar_cminpack();

  /*
	Tcl_Interp *interp;
	interp = NKinterp;
	Para *p;
	Tcl_HashEntry *entryPtr;
	entryPtr=_getExistingPtr(interp, &paraHT, (char *) "p", (void**)&p);
	*/

  resetFlags();
	bestChain = func->getBestChain();
	// Update the hash table instead of NKparams. Type casting is necessary.
	//bestChisq = func->recoverBestParams((double*)p);
	bestChisq = func->recoverBestParams(&g_ParaStruct);
  printf("%s Xr: %g\n", bestChain, bestChisq);
  updatelinks(bestChisq, bestChain);
  delete lmdif;
  delete func;
  fprintf(ntp->fp, "finished\n"); fflush(ntp->fp);
  delete ntp;
  return NULL;
}

/****************************************************************************
http://man7.org/linux/man-pages/man3/pthread_create.3.html

The pthread_create() function starts a new thread in the calling
process. The new thread starts execution by invoking
start_routine(); arg is passed as the sole argument of
start_routine().

The new thread terminates in one of the following ways:
* It calls pthread_exit(3), specifying an exit status value that is
  available to another thread in the same process that calls
  pthread_join(3).

* It returns from start_routine().  This is equivalent to calling
  pthread_exit(3) with the value supplied in the return statement.

* It is canceled (see pthread_cancel(3)).

* Any of the threads in the process calls exit(3), or the main thread
  performs a return from main().  This causes the termination of all
  threads in the process.

Before returning, a successful call to pthread_create() stores the ID
of the new thread in the buffer pointed to by thread; this identifier
is used to refer to the thread in subsequent calls to other pthreads
functions.
*****************************************************************************/
int YF_fitdata(ClientData clientData, Tcl_Interp *interp, int objc, Tcl_Obj *const objv[])
{
  NfitThreadPars *s;
  Data *dp;
  ModelCalculator *mc;
  //Para *p;
  int iter=10;
  pthread_t tid;

  // get the pointer to point to corresponding objects
  Tcl_HashEntry *entryPtr;
  entryPtr=_getExistingPtr(interp, &datasetHT, objv[1], (void**)&dp);
  if( entryPtr == NULL || dp == NULL ) return TCL_ERROR;
  entryPtr=_getExistingPtr(interp, &modelccHT, objv[2], (void**)&mc);
  if( entryPtr == NULL || mc == NULL ) return TCL_ERROR;
  //entryPtr=_getExistingPtr(interp, &paraHT, objv[3], (void**)&p);
  //if( entryPtr == NULL || p == NULL ) return TCL_ERROR;
  Tcl_GetIntFromObj( interp, objv[4], &iter );

  s = new NfitThreadPars;
  s->iter = ( iter > 100 ? 100 : iter );
  s->dp = dp;
  mc->paraSet(&g_ParaStruct); // sync the model with the parameter set.
  mc->setOmega(g_omega); // set the angle of incidence, omega
  s->mc = mc;
  s->p = &g_ParaStruct;
  s->fp = stdout;

  Tcl_EvalEx(interp, "set stopflag 0", -1, TCL_EVAL_GLOBAL);

  Pthread_create(&tid, NULL, thread_lmdif, (void *)s);

  return TCL_OK;
}


/****************************************************************************************
This updates linked parameters in param window. For thread safety, it should
only be invoked through the thread-safe updatelinks function when invoked
from the fitting thread as opposed to the Tcl event loop thread.
****************************************************************************************/
int updatelinks_async(ClientData clientData, Tcl_Interp *interp, int code)
{
  double xisq = update->xisq;
  //char *chain = update->chain;
  interp = NKinterp;
  Para *p;
  _getExistingPtr(interp, &paraHT, (char *) "p", (void**)&p);
  char cmd[20];
  /*
  for(int i=0;i<18;i++){
    NKparams[i]=((double*)p)[i];
    sprintf(cmd,"NKparams(%d)",i);
    Tcl_UpdateLinkedVar(interp, cmd);
  }
  //YFTcl_EvalEx(interp,256,"updateLambdaEta %g %g %g %g",NKparams[0],NKparams[1],NKparams[5],NKparams[17]);
  int i=0;
  NKparams[i]=int(1e17*((double*)p)[i]+0.5)*1e-17;
  sprintf(cmd,"NKparams(%d)",i);
  Tcl_UpdateLinkedVar(interp, cmd);
  i=1;
  NKparams[i]=int(1e-9*((double*)p)[i]+0.5)*1e9;
  sprintf(cmd,"NKparams(%d)",i);
  Tcl_UpdateLinkedVar(interp, cmd);
  */

  // empty the string stream, convert double to string, and then call Tcl_SetVar2 function
  std::ostringstream strs;
  for (int i = 0; i < Para::s_numParams; i++) {
    strs.str("");
    strs << g_ParaStruct.getValue(i);
    std::string newValue = strs.str();
    Tcl_SetVar2(interp, "paramArray", Varname[i], newValue.c_str(), TCL_GLOBAL_ONLY);
  }

  xisquare=xisq;
  sprintf(cmd, "xisquare");
  Tcl_UpdateLinkedVar(interp, cmd);
  //char temp[256];
  //sprintf(temp,"\"%s\"",chain);
  //if(recordflag) YFTcl_EvalEx(interp,256,"recordfit %s",temp);
  return code;
}

extern int Init_DSvectors(Tcl_Interp *interp);

extern "C" int Toad_Init(Tcl_Interp *interp){
  /*
    this function
    create coustomize commands in tcl
    link C variables with tcl variable
    call the starting script 'toad.tcl'
  */
  if (Tcl_Init(interp) != TCL_OK) {
    printf("%s\n", interp->result);
    return TCL_ERROR;
  }
  if (Tk_Init(interp) != TCL_OK) {
    printf("tk not loaded");
    return TCL_ERROR;
  }
  if (Blt_Init(interp) !=TCL_OK) {
    return TCL_ERROR;
  }

  Tcl_InitHashTable(&datasetHT, TCL_STRING_KEYS);
  Tcl_InitHashTable(&modelccHT, TCL_STRING_KEYS);
  Tcl_InitHashTable(&paraHT, TCL_STRING_KEYS);

  Tcl_CreateObjCommand(interp,"DSfitting", YF_DSfitting, (ClientData) NULL, (Tcl_CmdDeleteProc *) NULL);

  Tcl_CreateObjCommand(interp,"tv::img", TV_ImgCmd, (ClientData) NULL, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateObjCommand(interp,"tv::colorsetup", TV_ColorSetup, (ClientData) NULL, (Tcl_CmdDeleteProc *) NULL);

  Tcl_LinkVar(interp,"colorBOT", (char*)&colorBOT, TCL_LINK_INT);
  Tcl_LinkVar(interp,"colorMIN", (char*)&colorMIN, TCL_LINK_INT);
  Tcl_LinkVar(interp,"colorMAX", (char*)&colorMAX, TCL_LINK_INT);
  Tcl_LinkVar(interp,"colorSAT", (char*)&colorSAT, TCL_LINK_INT);
  Tcl_LinkVar(interp,"colorFIL", (char*)&colorFIL, TCL_LINK_INT);

  Tcl_CreateObjCommand(interp,"plotcs", TV_PlotCS, (ClientData) NULL, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateObjCommand(interp,"plothide", TV_PlotHide, (ClientData) NULL, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateObjCommand(interp,"boxv", TV_BoxV, (ClientData) NULL, (Tcl_CmdDeleteProc *) NULL);

  Tcl_CreateObjCommand(interp,"setNFIT", setNFIT, (ClientData) NULL, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateObjCommand(interp,"updateParaObject", updateParaObject, (ClientData) NULL, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateObjCommand(interp,"readDataFromASCII", readDataFromASCII, (ClientData) NULL, (Tcl_CmdDeleteProc *) NULL);

  Tcl_CreateObjCommand(interp,"fitdata", YF_fitdata, (ClientData) NULL, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateObjCommand(interp,"nfitdata", YF_nfitdata, (ClientData) NULL, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateObjCommand(interp,"dataset", YF_dataset, (ClientData) NULL, (Tcl_CmdDeleteProc *) NULL);
  //Tcl_CreateObjCommand(interp,"paraset", YF_paraset, (ClientData) NULL, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateObjCommand(interp,"modelcc", YF_modelcc, (ClientData) NULL, (Tcl_CmdDeleteProc *) NULL);


  Tcl_CreateObjCommand(interp,"qzrexport", NK_qzrexport, (ClientData) NULL, (Tcl_CmdDeleteProc *) NULL);


  TIFFSetWarningHandler(NULL);

  Tcl_LinkVar(interp, "plotopt", (char*)&gData.plotopt, TCL_LINK_INT);
  Tcl_LinkVar(interp, "startflag", (char*)&startflag, TCL_LINK_INT);
  Tcl_LinkVar(interp, "stopflag", (char*)&stopflag, TCL_LINK_INT);
  Tcl_LinkVar(interp, "recordflag", (char*)&recordflag, TCL_LINK_INT);
  Tcl_LinkVar(interp, "tvar_xswath", (char*)&gData.xswath, TCL_LINK_INT);
  Tcl_LinkVar(interp, "tvar_xswath", (char*)&gData.yswath, TCL_LINK_INT);

  Tcl_LinkVar(interp, "xlow", (char*)&xlow, TCL_LINK_DOUBLE);
  Tcl_LinkVar(interp, "xhigh", (char*)&xhigh, TCL_LINK_DOUBLE);
  Tcl_LinkVar(interp, "nocz", (char*)&nocz, TCL_LINK_INT);
  Tcl_LinkVar(interp, "noscale", (char*)&noscale, TCL_LINK_INT);
  Tcl_LinkVar(interp, "kiyomask", (char*)&kiyomask, TCL_LINK_INT);
  Tcl_LinkVar(interp, "FLICAM", (char*)&FLICAM, TCL_LINK_INT);

  Tcl_LinkVar(interp, "NKfilter", (char*)&NKfilter, TCL_LINK_DOUBLE);
  Tcl_LinkVar(interp, "descending", (char*)&descending, TCL_LINK_INT);

  Tcl_LinkVar(interp, "dataset", (char*)&dset, TCL_LINK_STRING);

  Tcl_LinkVar(interp, "verticalfit", (char*)&verticalfit, TCL_LINK_INT);

  Tcl_LinkVar(interp, "tmpsigma", (char*)&tmpsigma, TCL_LINK_DOUBLE);
  Tcl_LinkVar(interp, "dupe", (char*)&dupe, TCL_LINK_DOUBLE);

  Tcl_LinkVar(interp, "aFactor", (char*)&aFactor, TCL_LINK_DOUBLE);
  Tcl_LinkVar(interp, "backgroundSigmaSquare", (char*)&backgroundSigmaSquare, TCL_LINK_DOUBLE);
  Tcl_LinkVar(interp, "updateSigma", (char*)&updateSigma, TCL_LINK_INT);
  Tcl_LinkVar(interp, "nthreads", (char*)&nthreads, TCL_LINK_INT);
  Tcl_LinkVar(interp, "elapsedMins", (char*)&elapsedMins, TCL_LINK_INT);
  Tcl_LinkVar(interp, "elapsedSecs", (char*)&elapsedSecs, TCL_LINK_INT);

	Tcl_LinkVar(interp, "omega", (char*)&g_omega, TCL_LINK_DOUBLE);

  Tcl_EvalEx(interp,"image create photo realimg0",-1,TCL_EVAL_GLOBAL);
  gData.tv_photoHd=Tk_FindPhoto(interp, "realimg0");

  char varname[20];
  sprintf(varname,"S");
  Tcl_LinkVar(interp,varname,(char *)&(setup.s), TCL_LINK_DOUBLE);
  sprintf(varname,"R");
  Tcl_LinkVar(interp,varname,(char *)&(setup.R), TCL_LINK_DOUBLE);
  sprintf(varname,"x0");
  Tcl_LinkVar(interp,varname,(char *)&(setup.x0), TCL_LINK_DOUBLE);
  sprintf(varname,"lambda");
  Tcl_LinkVar(interp,varname,(char *)&(setup.lambda), TCL_LINK_DOUBLE);
  sprintf(varname,"D");
  Tcl_LinkVar(interp,varname,(char *)&(setup.d), TCL_LINK_DOUBLE);
  sprintf(varname,"alpha");
  Tcl_LinkVar(interp,varname,(char *)&(setup.alpha), TCL_LINK_DOUBLE);
  sprintf(varname,"pz");
  Tcl_LinkVar(interp,varname,(char *)&(setup.pz), TCL_LINK_DOUBLE);
  sprintf(varname,"delts");
  Tcl_LinkVar(interp,varname,(char *)&(setup.delts), TCL_LINK_DOUBLE);
  //added by NC for ref. index, 3.08.04
  sprintf(varname,"refn");
  Tcl_LinkVar(interp,varname,(char *)&(setup.refn), TCL_LINK_DOUBLE);
  //end of addition NC 3.08.04
  sprintf(varname,"chi");
  Tcl_LinkVar(interp,varname,(char *)&(tcl_chi), TCL_LINK_DOUBLE);
  for(int i=0;i<32;i++){
    sprintf(varname,"boolx(%d)",i);
    Tcl_LinkVar(interp,varname,(char *)&(boolx[i]), TCL_LINK_INT);
  }
  for(int i=0;i<4;i++){
    sprintf(varname,"xpin(%d)",i);
    Tcl_LinkVar(interp,varname,(char *)&(xpin[i]), TCL_LINK_INT);
    sprintf(varname,"step(%d)",i);
    Tcl_LinkVar(interp,varname,(char *)&(step[i]), TCL_LINK_DOUBLE);
    sprintf(varname,"sxpin(%d)",i);
    Tcl_LinkVar(interp,varname,(char *)&(sxpin[i]), TCL_LINK_INT);
  }
  for(int i=0;i<18;i++){
    sprintf(varname,"NKparams(%d)",i);
    Tcl_LinkVar(interp,varname,(char *)&(NKparams[i]), TCL_LINK_DOUBLE);
  }
  sprintf(varname,"xisquare");
  Tcl_LinkVar(interp,varname,(char *)&(xisquare), TCL_LINK_DOUBLE);

  sprintf(varname,"STresult");
  Tcl_LinkVar(interp,varname,(char *)&(STresult), TCL_LINK_DOUBLE);

  Init_DSvectors(interp);

  NKinterp=interp;
  update = new UpdateStruct;
  /*Set up handlers for asynchronous updates to the Tcl script invoked by
    threads other than the one running the Tcl event loop.
  */
  updateHandler = Tcl_AsyncCreate(updatelinks_async, NULL);
  resetHandler = Tcl_AsyncCreate(resetFlags_async, NULL);
  timerHandler = Tcl_AsyncCreate(updateTimer_async, NULL);

  Tcl_EvalEx(interp,"source tcl/toad.tcl",-1,TCL_EVAL_GLOBAL);
  boxvec.init(interp);
  return TCL_OK;
}

/*Thread safe means of updating variables associated with a fitting run in the
  Tcl script from the fitting thread.
*/
void updatelinks(double xisq, char *chain){
  update->xisq = xisq;
  update->chain = chain;
  Tcl_AsyncMark(updateHandler);
}

//function controlling the refinement (filtering the outstanding points)
double refine(double theor, double exper){
  double tmp;
  //tmp=fabs(theor-exper);
  if (theor != 0) {tmp=fabs(theor-exper)/theor;}
  else {tmp=fabs(exper);}
  return tmp;
}

/*
void kiyoMask ()
{
  if (kiyomask){
    // find lower left corner of the fitting box
    int x0 = dp->qs[0].sdata[0].qx;
    int y0 = dp->qs[0].qz;
    cout << "x0 = " << x0 << "   y0 = " << y0 << endl;
    //cout << "-----Data points that will be ignored in the chi squared evaluation-----" << endl;

    int lowx, highx, tempy;
    //std::string line;
    std::ifstream infile("kiyomask.txt");

    while (infile >> tempy >> lowx >> highx){
    //cout << "qz = " << dp->qs[tempy-y0].qz << "   qx = ";
    for (int i = lowx; i <=highx; i++){
                dp->qs[tempy-y0].sdata[i-x0].sigma = 1000000000;
                //cout << dp->qs[tempy-y0].sdata[i-x0].qx << " ";
            }
        }

  }
}
*/

int YF_nfitdata(ClientData clientData, Tcl_Interp *interp,
int objc, Tcl_Obj *const objv[]){
  Data *dp;
  ModelCalculator *mc;
  Para *p;
  int iter=10;

  Tcl_HashEntry *entryPtr;
  entryPtr=_getExistingPtr(interp, &datasetHT, objv[1], (void**)&dp);
  if( entryPtr == NULL || dp == NULL ) return TCL_ERROR;
  entryPtr=_getExistingPtr(interp, &modelccHT, objv[2], (void**)&mc);
  if( entryPtr == NULL || mc == NULL ) return TCL_ERROR;
  entryPtr=_getExistingPtr(interp, &paraHT, objv[3], (void**)&p);
  if( entryPtr == NULL || p == NULL ) return TCL_ERROR;
  Tcl_GetIntFromObj( interp, objv[4], &iter );

  iter = ( iter>100 ? 100 : iter );
  mc->paraSet(p);
  FuncLmdif *func = new FuncLmdif;
  func->setData(dp);
  func->setMC(mc);
  func->setPara(p);

  Lmdif *lmdif = new  Lmdif;
  cout<<"iter: "<<iter<<endl;
  lmdif->setup(FuncLmdif::WrapperFunclmdif, func->npoint(), p->nfit, iter, 0, 0, 0, 0.1, p->epsfcn, func);

  lmdif->init(p->xp, p->nfit);
  lmdif->fit();

  delete lmdif;
  delete func;
  fprintf(stdout, "***********finished"); fflush(stdout);

  return TCL_OK;
}

