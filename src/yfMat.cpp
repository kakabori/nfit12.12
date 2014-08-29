/*****************************
This is something I did for fun
I will document this part later
******************************/

#include "tvImg.h"
#include "yfMat.h"

#define s2(a,b) {if (b<a) {t=a; a=b; b=t;}}
#define mn3(a,b,c) s2(a,b); s2(a,c);
#define mx3(a,b,c) s2(b,c); s2(a,c);
#define mnmx3(a,b,c) mx3(a,b,c); s2(a,b);
#define mnmx4(a,b,c,d) s2(a,b); s2(c,d); s2(a,c); s2(b,d);
#define mnmx5(a,b,c,d,e) s2(a,b); s2(c,d); mn3(a,c,e); mx3(b,d,e);
#define mnmx6(a,b,c,d,e,f) s2(a,d); s2(b,e); s2(c,f);\
                            mn3(a,b,c); mx3(d,e,f);
template<class T>
T med3x3(T *b1, T *b2, T *b3)
/*
 * Find median on a 3x3 input box of integers.
 * b1, b2, b3 are pointers to the left-hand edge of three
 * parallel scan-lines to form a 3x3 spatial median.
 * Rewriting b2 and b3 as b1 yields code which forms median
 * on input presented as a linear array of nine elements.
 */
{
  register T t, r1, r2, r3, r4, r5, r6;
  r1 = *b1++; r2 = *b1++; r3 = *b1++;
  r4 = *b2++; r5 = *b2++; r6 = *b2++;
  mnmx6(r1, r2, r3, r4, r5, r6);
  r1 = *b3++;
  mnmx5(r1, r2, r3, r4, r5);
  r1 = *b3++;
  mnmx4(r1, r2, r3, r4);
  r1 = *b3++;
  mnmx3(r1, r2, r3);
  return(r2);
}

template<class T>
void _medianfilter_inplace(T *dp, int strip, int w, int h){
  T *b1, *b2, *b3;
  b1 = dp;
  b2 = b1 + strip;
  b3 = b2 + strip;
  dp += w - 2;
  for(int k=2; k<h; k++){
    for(; b1 < dp; b1++, b2++, b3++)
      b2[1] = med3x3(b1, b2, b3);

    b1 = b1 + strip - w + 2;
    b2 = b1 + strip;
    b3 = b2 + strip;
    dp = dp + strip;
  }
}

template<class T>
void _medianfilter(T *dp, int strip, int w, int h){
  T *b1, *b2, *b3, *tp1, *tp2;
  b1 = dp;
  b2 = b1 + strip;
  b3 = b2 + strip;
  dp += w - 2;
  tp1 = (T*)malloc(sizeof(T)*(w-2));
  tp2 = (T*)malloc(sizeof(T)*(w-2));

  for(int i=0; b1 < dp; i++, b1++, b2++, b3++)
    tp1[i] = med3x3(b1, b2, b3);

  for(int k=3; k<h; k++){
    b1 = b1 + strip - w + 2;
    b2 = b1 + strip;
    b3 = b2 + strip;
    dp = dp + strip;
    for(int i=0; b1 < dp; i++, b1++, b2++, b3++)
      tp2[i] = med3x3(b1, b2, b3);
    
    memcpy( dp-w+3, tp1, sizeof(T)*(w-2) );
    b3 = tp1; tp1 = tp2; tp2 = b3;
  }
  memcpy( b2-w+3, tp1, sizeof(T)*(w-2) );
  free(tp1);
  free(tp2);
}








/* reorgnize the way data saved so that FFT can be performed */
template<class T> void _fft_convert(T *ptr, double *d, int width, int height, int jump, int oriwidth){
  for(int i=0;i<height;i++){
    for(int j=0;j<width;j++){
      d[i*jump+j]=ptr[i*oriwidth+j];
    }
  }
}

/* helper function for FFT , kernel smoothing */
static double*  _FFT(void *dataPtr, int oriwidth, YFRect<int> relative, int type,  double *hs, int lowcut){
  int jump,i,ij,j;
  int lnum;
  int width=relative.width();
  int length=relative.height();
  double *a,x;
  fftw_complex *A;
  rfftwnd_plan p, pinv;

  lnum=width*length;
  a=(double *)calloc((long)length*2*(width/2+1),sizeof(double));
  jump=2*(width/2+1);
  switch(type){
  case yfINT16:
    _fft_convert((short *)dataPtr+oriwidth*relative.lly()+relative.llx(), a, width, length, jump, oriwidth);
    break;
  case yfFLOAT:
    _fft_convert((float *)dataPtr+oriwidth*relative.lly()+relative.llx(), a, width, length, jump, oriwidth);
    break;
  }

  p    = rfftw2d_create_plan(length, width, FFTW_REAL_TO_COMPLEX,FFTW_ESTIMATE|FFTW_IN_PLACE);
  pinv = rfftw2d_create_plan(length, width, FFTW_COMPLEX_TO_REAL,FFTW_ESTIMATE|FFTW_IN_PLACE);
  if(a==NULL) printf("no space\n");
  rfftwnd_one_real_to_complex(p, (fftw_real *) a, NULL);
  A = (fftw_complex*) a;
  for(i = 0; i < length; ++i){
    for (j = 0; j < width/2+1; ++j) {
      x=exp(-2*(hs[0]*j/width*3.1416)*(hs[0]*j/width*3.1416));
      if(i>0&&i<length/2) x*=exp(-2*(hs[1]*i/length*3.1416)*(hs[1]*i/length*3.1416));
      else if(i>0) x*=exp(-2*(hs[1]*(length-i)/length*3.1416)*(hs[1]*(length-i)/length*3.1416));
      if(abs(length-i)<lowcut || i<lowcut || j<lowcut) x=0;
      ij = i*(width/2+1) + j;
      A[ij].re *=x; //			A[ij].re /=x;
      A[ij].im *=x;//			A[ij].im /=x;
    }
  }
  rfftwnd_one_complex_to_real(pinv, A, NULL);
  rfftwnd_destroy_plan(p);rfftwnd_destroy_plan(pinv);
  return a;
}




void YFMatView::medianFilterInplace(){
  switch(type){
  case yfUINT16: {
    _medianfilter_inplace( dp.uint16p + offset, strip, w, h);
    break;}
  case yfINT16: {
    _medianfilter_inplace( dp.int16p + offset, strip, w, h);
    break;}
  case yfFLOAT: {
    _medianfilter_inplace( dp.floatp + offset, strip, w, h);
    break;}
  }
}
void YFMatView::medianFilter(){
  switch(type){
  case yfUINT16: {
    _medianfilter( dp.uint16p + offset, strip, w, h);
    break;}
  case yfINT16: {
    _medianfilter( dp.int16p + offset, strip, w, h);
    break;}
  case yfFLOAT: {
    _medianfilter( dp.floatp + offset, strip, w, h);
    break;}
  }
}


template<class T>
void _peakCenter(T *dp, int strip, int w, int h){
  
}

template<class T>
inline float _convstep(T *p, int strip, float *f, int sz){
  float res = 0;
  for(int i=0; i<sz; i++){
    for(int j=0; j<sz; j++){
      res += (*p++) * (*f++);
    }
    p = p + strip - sz;
  }
  return res;
}

/*template<class T>
void _convfilter(T *dp, T *ndp, int strip, int nstrip, int w, int h, float *f, int sz){
  int nw = w - sz + 1;
  int nh = h - sz + 1;
  int jump = strip - nw;
  int njump = nstrip -nw;
  ndp = ndp + sz/2*(nstrip+1);

  for(int i=0; i<nh; i++){
    for(int k=0; k<nw; k++){
      *ndp++ = (T)_convstep(dp++, strip, f, sz);
    }
    ndp = ndp + njump;
    dp  =  dp + jump;
  }
}
*/

template<class T>
void _convfilter(T *dp, int strip, int w, int h, float *f, int sz){
  int nw = w - sz + 1;
  int nh = h - sz + 1;
  int jump = strip - nw;
  cout<<strip<<' '<<nw<<' '<<nh<<' '<<sz<<endl;
  T *p =dp + (sz/2)*(strip+1);
  T *ndp = (T*)malloc(sizeof(T)*nw*nh);
  for(int i=0; i<nh; i++){
    for(int k=0; k<nw; k++){
      *ndp++ = (T)_convstep(dp++, strip, f, sz);
    }
    dp  =  dp + jump;
  }
  ndp -= nw*nh;
  for(int i=0; i<nh; i++){
    memcpy(p, ndp, sizeof(T)*nw);
    p = p + strip;
    ndp = ndp + nw;
  }
  free(ndp-nw*nh);
}

float LoG1[]={
   0,  0, -3, -2, -2, -2, -3,  0,  0,
   0, -2, -3, -5, -5, -5, -3, -2,  0,
  -3, -3, -5, -3,  0, -3, -5, -3, -3,
  -2, -5, -3, 12, 24, 12, -3, -5, -2,
  -2, -5,  0, 24, 40, 24,  0, -5, -2,
  -2, -5, -3, 12, 24, 12, -3, -5, -2,
  -3, -3, -5, -3,  0, -3, -5, -3, -3,
   0, -2, -3, -5, -5, -5, -3, -2,  0,
   0,  0, -3, -2, -2, -2, -3,  0,  0};
float LoG[]={
   0, -1, -1, -2, -2, -2, -1, -1,  0,
  -1, -2, -4, -5, -5, -5, -4, -2, -1,
  -1, -4, -5, -3,  0, -3, -5, -4, -1,
  -2, -5, -3, 12, 24, 12, -3, -5, -2,
  -2, -5,  0, 24, 40, 24,  0, -5, -2,
  -2, -5, -3, 12, 24, 12, -3, -5, -2,
  -1, -4, -5, -3,  0, -3, -5, -4, -1,
  -1, -2, -4, -5, -5, -5, -4, -2, -1,
   0, -1, -1, -2, -2, -2, -1, -1,  0};

void YFMatView::convfilter(){
  float *f = LoG;
  int sz = 9;
  switch(type){
  case yfUINT16: {
    _convfilter( dp.uint16p + offset, strip, w, h, f, sz);
    break;}
  case yfINT16: {
    _convfilter( dp.int16p + offset, strip, w, h, f, sz);
    break;}
  case yfFLOAT: {
    _convfilter( dp.floatp + offset, strip, w, h, f, sz);
    break;}
  }
}

void YFMatView::peakCenter(){
}


/* this function provide kernel smoothing for now
 the result is save as floating point */
void YFMatView::FFTsmooth(double *hs){
  int jump=2*(w/2+1), lnum=w*h, i, j;
  double *a = NULL;
  YFRect<int> pos;
  pos.set(0, 0, w, h);
  switch(type){
  case yfINT16: {
    short *sp = dp.int16p;
    a=_FFT(dp.int16p+offset, w, pos, type, hs, -1);
    for(i=0;i<h;i++){
      for(j=0;j<w;j++){
	sp[i*w+j]=(short)(a[i*jump+j]/lnum);
      }
    }
    break;}
  case yfFLOAT: {
    float *sp = dp.floatp;
    a=_FFT(dp.floatp+offset, w, pos, type, hs, -1);
    for(i=0;i<h;i++){
      for(j=0;j<w;j++){
	sp[i*w+j]=(float)(a[i*jump+j]/lnum);
      }
    }
    break;}
  }
  free(a);
}
/* this function is to test the effect of high frequency cutoff
 in addition to kernel smoothing
void YFImage2D::filter(double *hs, YFRectangle rect, int lowcut){
  int jump, lnum, i, j, rectwidth, rectheight, offset;
  double *a;

  rect.shift(position.lleft);
  rectwidth=rect.width();
  rectheight=rect.height();
  jump=2*(rectwidth/2+1);
  lnum=rectwidth*rectheight;
  offset=rect.offset(width);


  a=_FFT(dataPtr,width, rect, type, hs, lowcut);
  if(type==0){
    short *sp;
    sp=(short *)dataPtr;
    sp+=offset;
    for(i=0;i<rectheight;i++){
      for(j=0;j<rectwidth;j++){
	sp[i*width+j]=(short)(a[i*jump+j]/lnum);
      }
    }
  }else if(type==2){
    df_prec *sp;
    sp=(df_prec *)dataPtr;
    sp+=offset;
    for(i=0;i<rectheight;i++){
      for(j=0;j<rectwidth;j++){
	sp[i*width+j]=(df_prec)(a[i*jump+j]/lnum);
      }
    }
  }
  free(a);
  }*/
