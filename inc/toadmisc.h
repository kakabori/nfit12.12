#ifndef _TOADMISC_HEADER_
#define _TOADMISC_HEADER_

#include <cstdlib>
#include <cstring>
#include <cstdio>
#include <tcl.h>
#include <tk.h>
#include <blt.h>
#include <tiff.h>
#include <tiffio.h>
#include <sstream>
#include <iostream>
#include <stdlib.h>
#include <vector>
#include "mydll.h"

using namespace std;

/* simple array operations */
template<class T, class A>
inline void yf_addConst(T* data, A a, int sz){
  for(int i=0; i<sz; i++) data[i]=(T)(data[i]+a);
}

template<class T, class A>
inline void yf_multiConst(T* data, A a, int sz){
  for(int i=0; i<sz; i++) data[i]=(T)(data[i]*a);
}

template<class T, class A>
inline void yf_addArr(T* data, A* a, int sz){
  for(int i=0; i<sz; i++) data[i]=(T)(data[i]+a[i]);
}

template<class T, class A>
inline void yf_multiArr(T* data, A* a, int sz){
  for(int i=0; i<sz; i++) data[i]=(T)(data[i]*a[i]);
}

template<class T, class A>
inline void yf_subArr(T* data, A* a, int sz){
  for(int i=0; i<sz; i++) data[i]=(T)(data[i]-a[i]);
}

template<class T, class A>
inline void nk_addArr(T* ptrT, int wT, A* ptrA, int wA, int w, int h){
  T *pT=ptrT;
  A *pA=ptrA;
  for(int i=0; i<h; i++){
    for(int k=0; k<w; k++){
      pT[k]+=(T)(pA[k]);
    }
    pT+=wT;
    pA+=wA;
  }
}
template<class T, class A>
inline void nk_subArr(T* ptrT, int wT, A* ptrA, int wA, int w, int h){
  T *pT=ptrT;
  A *pA=ptrA;
  for(int i=0; i<h; i++){
    for(int k=0; k<w; k++){
      pT[k]-=(T)(pA[k]);
    }
    pT+=wT;
    pA+=wA;
  }
}
template<class T, class A>
inline void nk_multiArr(T* ptrT, int wT, A* ptrA, int wA, int w, int h){
  T *pT=ptrT;
  A *pA=ptrA;
  for(int i=0; i<h; i++){
    for(int k=0; k<w; k++){
      pT[k]=(T)(pT[k]*pA[k]);
    }
    pT+=wT;
    pA+=wA;
  }
}
template<class T, class A>
inline void nk_divArr(T* ptrT, int wT, A* ptrA, int wA, int w, int h){
  T *pT=ptrT;
  A *pA=ptrA;
  for(int i=0; i<h; i++){
    for(int k=0; k<w; k++){
      pT[k]=(T)(pT[k]/pA[k]);
    }
    pT+=wT;
    pA+=wA;
  }
}

template<class T, class A>
inline void yf_divArr(T* data, A* a, int sz){
  for(int i=0; i<sz; i++) data[i]=(T)(data[i]/a[i]);
}

template<class T> void yf_rotateImg(T **ptrptr, double rad,
				    int w, int h,
				    int nw, int nh,
				    double cx, double cy,
				    double ncx, double ncy){
  /* rotate an image
     rad: angle in radian
     w  : width of the image
     h  : height of the image
     nw : new width of the image
     nh : new height of the image
     (cx,cy) : origin of the coordinates for the old image
     (ncx,ncy): origin of the coordinates for the new image
   */
  double cosT=cos(rad);
  double sinT=sin(rad);
  T *nptr;
  T *ptr=(*ptrptr);
  nptr=(T *)calloc(sizeof(T), nw*nh);

  for(int i=0;i<nw;i++){
    for(int k=0;k<nh;k++){
      int intx, inty;
      double x, y;
      double tempx, tempy;
      int dex;

      tempx=i-ncx;
      tempy=k-ncy;
      x=cosT*tempx + sinT*tempy + cx;
      y=cosT*tempy - sinT*tempx + cy;
      x=(x-(intx=(int)floor(x)));
      y=(y-(inty=(int)floor(y)));
      if(intx<0 || inty<0 || intx>w-2 || inty>h-2) *(nptr+k*nw+i)=0;
      else{
	dex=inty*w+intx;
	*(nptr+k*nw+i)=(T)(*(ptr+dex)*(1-x)*(1-y)+
			   *(ptr+dex+1)*x*(1-y)+
			   *(ptr+dex+w)*(1-x)*y+
			   *(ptr+dex+w+1)*x*y);
      }
    }
  }

  free(ptr);
  *ptrptr=nptr;
}
//added by NC 4.21.03
template<class T> void nc_rotateImg(T **ptrptr,int direction, int width, int height, double ox, double oy){
	//printf("i'm really doing rotate90 now\n");

	//copied from _rotate, revised cosT & sinT 's values
	//cosT in fact can be omitted, and sinT = 1|-1
  int ip,kp,dex;
  double rx,ry,x,y;
  double cosT=0;//double cosT=1.0/sqrt(1+tanT*tanT); 
  double sinT=direction;//double sinT=tanT*cosT;
  T *nptr;
  T *ptr=(*ptrptr);
  nptr=(T *)malloc(sizeof(T)*width*height);
  for(int i=0;i<width;i++){
    for(int k=0;k<height;k++){
      rx=i-ox;
      ry=k-oy;
      x=cosT*rx+sinT*ry+ox;
      y=cosT*ry-sinT*rx+oy;
      x=(x-(ip=(int)floor(x)));
      y=(y-(kp=(int)floor(y)));
      if(ip<0 || kp<0 || ip>width-2 || kp>height-2) *(nptr+k*width+i)=1;
      else{
	dex=kp*width+ip;
	*(nptr+k*width+i)=(T)(*(ptr+dex)*(1-x)*(1-y)+
			      *(ptr+dex+1)*x*(1-y)+
			      *(ptr+dex+width)*(1-x)*y+
			      *(ptr+dex+width+1)*x*y);
      }
    }
  }
  free(ptr);
  *ptrptr=nptr;
}
//end of NC

template<class T> void yf_affineImg(T **ptrptr, double rad, double scale,
				    int w, int h,
				    int nw, int nh,
				    double tx, double ty){
  /* scaling , rotation and translation 
     obviously, this supersede the 'yf_rotateImg'
   */
  double cosT=cos(rad) * scale;
  double sinT=sin(rad) * scale;
  T *nptr;
  T *ptr=(*ptrptr);
  nptr=(T *)calloc(sizeof(T), nw*nh);

  for(int i=0;i<nw;i++){
    for(int k=0;k<nh;k++){
      int intx, inty;
      double x, y;
      int dex;

      x=cosT * i + sinT * k + tx;
      y=cosT * k - sinT * i + ty;
      x=(x-(intx=(int)floor(x)));
      y=(y-(inty=(int)floor(y)));
      if(intx<0 || inty<0 || intx>w-2 || inty>h-2) *(nptr+k*nw+i)=0;
      else{
	dex=inty*w+intx;
	*(nptr+k*nw+i)=(T)(*(ptr+dex)*(1-x)*(1-y)+
			   *(ptr+dex+1)*x*(1-y)+
			   *(ptr+dex+w)*(1-x)*y+
			   *(ptr+dex+w+1)*x*y);
      }
    }
  }

  free(ptr);
  *ptrptr=nptr;
}

template<class T> inline void yf_swap(T& a, T& b){
  T temp;
  temp=a;
  a=b;
  b=temp;
}

template<class T> void yf_flipx(T *ptr, int width, int height){
  /* flip images horizontally */
  int halfw=width/2;
  for(int i=0; i<height; i++)
    for(int k=0; k<halfw; k++)
      yf_swap(ptr[width*i+k], ptr[width*i+width-k-1]);
}

template<class T> void yf_flipy(T *ptr, int width, int height){
  /* flip images vertically */
  int halfh=height/2;
  for(int i=0;i<width;i++)
    for(int k=0;k<halfh;k++)
      yf_swap(ptr[width*k+i], ptr[width*(height-k-1)+i]);
}

template<class T> void yf_transpose(T *ptr, int width, int height){
  /* transpose image */
  T* tp=(T*)malloc(sizeof(T)*width*height);
  memcpy(tp, ptr, sizeof(T)*width*height);
  for(int i=0;i<width;i++)
    for(int k=0;k<height;k++)
      ptr[height*i+k]= tp[width*k+i];
  free(tp);
}

template<class A, class B> B* yf_convert(A *a, B *b, int sz){
  //covert from type A to type B
  if(sizeof(B)>sizeof(A)){
    b=(B *)realloc(a,sizeof(B)*sz);
    if(b==NULL){
      printf("reallocation fails\n");
      return NULL;
    }
    a=(A *)b;
    for(int i=sz-1;i>=0;i--) b[i]=(B)a[i];
  } else {
    b=(B *)a;
    for(int i=0;i<sz;i++) b[i]=(B)a[i];
    b=(B *)realloc(b,sizeof(B)*sz);
  }
  return b;
}

template<class T> void nk_shift(T *ptr, YFDuple<int> s, int width, int height){
  s.y=-1*s.y;
  if ((height-abs(s.y)>0) || (width-abs(s.x)>0)){
    if (s.y>0) {
	memmove(ptr+s.y*width, ptr, width*sizeof(T)*(height-s.y));
	memset(ptr,0,width*sizeof(T)*s.y);
    } else {
	memmove(ptr,ptr-s.y*width, width*sizeof(T)*(height+s.y));
	memset(ptr+width*(height+s.y),0,width*sizeof(T)*(-1*s.y));
    }
    T *tmptr;
    for(int i=0;i<height;i++){
      tmptr=ptr+i*width;
      if(s.x>0) {
	  memmove(tmptr+s.x,tmptr,sizeof(T)*(width-s.x));
	  memset(tmptr,0,sizeof(T)*s.x);
      } else {
	  memmove(tmptr,tmptr-s.x,sizeof(T)*(width+s.x));
	  memset(tmptr+(width+s.x),0,sizeof(T)*(-1*s.x));
      }
    }
  } else cout<<"cannot shift"<<endl;
}

Tcl_HashEntry* _getExistingPtr(Tcl_Interp *interp, Tcl_HashTable *htptr, Tcl_Obj *obj, void** ptrptr);
Tcl_HashEntry* _getExistingPtr(Tcl_Interp *interp, Tcl_HashTable *htptr, char * name , void** ptrptr);
void YFTcl_EvalEx(Tcl_Interp* interp, size_t sz, const char*fmt, ...);
void polint(const double *xa, const double *ya,int n,double x,double *y,double *dy);
void YFTcl_SetResult(Tcl_Interp* interp, size_t sz, const char*fmt, ...);


enum TermAttr { TERMATTR_RESET	=	0,
		TERMATTR_BRIGHT =	1,
		TERMATTR_DIM	=	2,
		TERMATTR_UNDERLINE =	3,
		TERMATTR_BLINK	=	4,
		TERMATTR_REVERSE=	7,
		TERMATTR_HIDDEN	=	8
};

enum TermColor {TERMCOLOR_BLACK 	=	0,
		TERMCOLOR_RED	=	1,
		TERMCOLOR_GREEN	=	2,
		TERMCOLOR_YELLOW	=	3,
		TERMCOLOR_BLUE	=	4,
		TERMCOLOR_MAGENTA	=	5,
		TERMCOLOR_CYAN	=	6,
		TERMCOLOR_WHITE	=	7
};

void textcolor(FILE* fp, int attr, int fg, int bg);
void textcolor(ostream& out, int attr, int fg, int bg);

#endif
