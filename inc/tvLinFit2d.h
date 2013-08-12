#ifndef _TVLINFIT2D_
#define _TVLINFIT2D_

#include <cstdlib>
#include <cstdio>
#include <cmath>

#define dfprec float

void svdfit(int ndata,dfprec *a,int ma,dfprec *u,dfprec *v,dfprec *w, dfprec *b);
extern double tmpsigma;

//ma: num of linear parameter
//(*funcs)(x, y, result)

template<class T> void  Fit2D_poly(T *y, int imgwidth, dfprec *a,
				   YFRect<int> rect1, YFRect<int> rect2,
				   int ma, void (*funcs)(dfprec,dfprec,dfprec*,int), int pn, int subopt){
  /*
    This fit area defined by two boxes using polynomials
    too much to explain the SVD algorithm (see numerical recipe)
  */

/*
subopt controls background-subtraction options as well as helps to speed up fitting
transformed into binary mode it represents 11 switching variables
e.g.: "10010000100" = 4 + 128 + 1024 = 1156
*/
//  enum{subtractR2Overlap=1,subtractR1Overlap=2,subtractR2=4,subtractR1=8,subtractOverlap=16,
//       replaceByApprox=128,funcs_1=1024};
  int nout, j, ndata;
  int n,w,h, overlap_flag=1;
  dfprec tmp,*u, *b, *v, *wsp, *afunc;
  YFRect<int> overlap;
  dfprec sum, sum2;
  //printf("%d %d %d %d \t%d %d %d %d\n", rect1.llx(), rect1.lly(), rect1.urx(), rect1.ury(), rect2.llx(), rect2.lly(), rect2.urx(), rect2.ury());
  overlap_flag=rect1.overlap(rect2, overlap);
  
  if(overlap_flag) ndata=rect2.area()-overlap.area();
  else ndata=rect2.area();

  u=(dfprec *)malloc(sizeof(dfprec)*ma*ndata);
  v=(dfprec *)malloc(sizeof(dfprec)*ma*ma);
  b=(dfprec *)malloc(sizeof(dfprec)*ndata);
  wsp=(dfprec *)malloc(sizeof(dfprec)*ndata);
  afunc=(dfprec *)malloc(sizeof(dfprec)*ma);

  if(u==NULL || v==NULL || b==NULL || wsp==NULL || afunc==NULL){
    printf("can not allocate enought memory in Fit2D_poly\n");
    return;
  }

  n=0;
  for(h=rect2.lly(); h<rect2.ury(); h++){//initialize arrays for fitting, iterate throngh pixels corresponding box2-box
    for(w=rect2.llx(); w<rect2.urx(); w++){
      if(overlap_flag && w<overlap.urx() && w>=overlap.llx() && h>=overlap.lly() && h<overlap.ury()) {w=overlap.urx()-1; continue;}
      (*funcs)(w, h, afunc,pn);
      for (j=0;j<ma;j++) (u+j*ndata)[n]=afunc[j];
      b[n]=y[h*imgwidth+w];
      n++;
    }
  }
  if(n!=ndata) printf("need debugging, report this please\n");
  svdfit(ndata,a,ma,u,v,wsp,b);
  sum2=0;
  for(h=rect2.lly(); h<rect2.ury(); h++){
    for(w=rect2.llx(); w<rect2.urx(); w++){
      if(overlap_flag && w<overlap.urx() && w>=overlap.llx() && h>=overlap.lly() && h<overlap.ury()) {w=overlap.urx()-1; continue;}
      (*funcs)(w, h, afunc,pn);
      for (sum=0.0,j=0;j<ma;j++) sum += a[j]*afunc[j];//sum is the value of ploy
      sum2+=fabs(y[h*imgwidth+w]-sum);//sum2 is total deviation
    }
  }
  sum2 /= ndata/2;
  n=0;nout=0;
  for(h=rect2.lly(); h<rect2.ury(); h++){//initialize arrays for fitting, iterate throngh pixels corresponding box2-box
    for(w=rect2.llx(); w<rect2.urx(); w++){
      if(overlap_flag && w<overlap.urx() && w>=overlap.llx() && h>=overlap.lly() && h<overlap.ury()) {w=overlap.urx()-1; continue;}
      (*funcs)(w, h, afunc,pn);
      for (sum=0.0,j=0;j<ma;j++) sum += a[j]*afunc[j];
      if( fabs(y[h*imgwidth+w]-sum) > sum2 ) {tmp=tmpsigma;nout++;} else tmp=1;
      for (j=0;j<ma;j++) (u+j*ndata)[n]=afunc[j]*tmp;
      b[n]=y[h*imgwidth+w]*tmp;
      n++;
    }
  }
  fprintf(stderr,"\r# nout: %d",nout);
  svdfit(ndata,a,ma,u,v,wsp,b);

  if(subopt&1){
    for(h=rect2.lly(); h<rect2.ury(); h++){//subtract the quadratic background: rect2-overlap
      for(w=rect2.llx(); w<rect2.urx(); w++){
	if(overlap_flag && w<overlap.urx() && w>=overlap.llx() && h>=overlap.lly() && h<overlap.ury()) {w=overlap.urx()-1; continue;}
	(*funcs)(w, h, afunc,pn);
	for (sum=0.0,j=0;j<ma;j++) sum += a[j]*afunc[j];
	y[h*imgwidth+w] =(T)(y[h*imgwidth+w]-sum);
      }
    }
  }
  if(subopt&2){
    for(h=rect1.lly(); h<rect1.ury(); h++){//subtract the quadratic background: rect1-overlap
      for(w=rect1.llx(); w<rect1.urx(); w++){
	if(overlap_flag && w<overlap.urx() && w>=overlap.llx() && h>=overlap.lly() && h<overlap.ury()) {w=overlap.urx()-1; continue;}
	(*funcs)(w, h, afunc,pn);
	for (sum=0.0,j=0;j<ma;j++) sum += a[j]*afunc[j];
	if(subopt & 128) y[h*imgwidth+w] =(T)(sum);
	else y[h*imgwidth+w] =(T)(y[h*imgwidth+w]-sum);
      }
    }
  }
  if(subopt&4){
    for(h=rect2.lly(); h<rect2.ury(); h++){//subtract the quadratic background: rect 2
      for(w=rect2.llx(); w<rect2.urx(); w++){
	(*funcs)(w, h, afunc,pn);
	for (sum=0.0,j=0;j<ma;j++) sum += a[j]*afunc[j];
	if(subopt & 128) y[h*imgwidth+w] =(T)(sum);
	else 	y[h*imgwidth+w] =(T)(y[h*imgwidth+w]-sum);
      }
    }
  }
  if(subopt&8){
    for(h=rect1.lly(); h<rect1.ury(); h++){//subtract the quadratic background: rect1
      for(w=rect1.llx(); w<rect1.urx(); w++){
	(*funcs)(w, h, afunc,pn);
	for (sum=0.0,j=0;j<ma;j++) sum += a[j]*afunc[j];
	if(subopt & 128) y[h*imgwidth+w] =(T)(sum);
	else 	y[h*imgwidth+w] =(T)(y[h*imgwidth+w]-sum);
      }
    }
  }
  if((subopt&16) && overlap_flag){
    printf("overlap\n");
    for(h=overlap.lly(); h<overlap.ury(); h++){//subtract the quadratic background: overlap
      for(w=overlap.llx(); w<overlap.urx(); w++){
	(*funcs)(w, h, afunc,pn);
	for (sum=0.0,j=0;j<ma;j++) sum += a[j]*afunc[j];
	if(subopt & 128) y[h*imgwidth+w] =(T)(sum);
	else 	y[h*imgwidth+w] =(T)(y[h*imgwidth+w]-sum);
      }
    }
  }
     
  free(afunc);
  free(u);
  free(v);
  free(b);
  free(wsp);
}

#endif
