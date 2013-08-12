#include <math.h>
#include <stdlib.h>
#include <limits.h>
#include <float.h>
#include <tcl.h>
#include <tk.h>
#include <blt.h>

#include "boxchain.h"

void BoxVec::updateBox(Box* bptr){
  /*
    pre: bptr points to initialized or updated Box object
    post: $imgWN updates the box marker
  */
  YFTcl_EvalEx(itp,
	       512,
	       "$imgWN marker configure %s -coords {%d %d %d %d %d %d %d %d %d %d} ; $imgWN marker configure _%s -coords {%d %d}",
	       bptr->name,
	       bptr->xc[0], bptr->yc[0],
	       bptr->xc[0], bptr->yc[1],
	       bptr->xc[1], bptr->yc[1],
	       bptr->xc[1], bptr->yc[0],
	       bptr->xc[0], bptr->yc[0],
 	       bptr->name,
	       bptr->minxc(), bptr->minyc() );
}



void BoxVec::init(Tcl_Interp* interp){
  /*
    pre: valid Tcl interpreter
    post: bvhash (Hash table for box marker name searching) initialized
  */
  itp=interp;
  Tcl_InitHashTable(&bvhash, TCL_STRING_KEYS);
  yhalo=xhalo=25;
  curb=NULL;
  nextname=1;
}

/*
  step1 : determine the purpose for the first mouse click (resize, move, new)
  step2 : keep track of the mouse movement
  step3 : finalize when mouse button released
*/


void BoxVec::step1(double dbx, double dby, char*nm){
  /*
    pre: valid Tcl interpreter
    post: bvhash (Hash table for box marker name searching) initialized
  */
  Box *bp;
  Box *minbp=NULL, *moveb=NULL;
  double mindxy=DBL_MAX;
  double mind =DBL_MAX;
  Tcl_HashSearch search;
  Tcl_HashEntry *entryPtr=Tcl_FirstHashEntry(&bvhash, &search);

  while(entryPtr!=NULL){
    double dx, dy;
    bp=(Box*)(entryPtr->clientData);
    dx=bp->nearx(dbx);
    dy=bp->neary(dby)*xhalo/yhalo;
    if( dx < mindxy && fabs(dby-bp->midyc()) < bp->height()/2.0 ) {mindxy=dx; moveb=bp;}
    if( dy < mindxy && fabs(dbx-bp->midxc()) < bp->width()/2.0  ) {mindxy=dy; moveb=bp;}
    if( (dx=sqrt(dx*dx+dy*dy*xhalo*xhalo/yhalo/yhalo)) < mind) {mind=dx; minbp=bp;}
    entryPtr=Tcl_NextHashEntry(&search);
  }
  if( mind<xhalo ){
    option=RESIZEBOX;
    curb=minbp;
    curb->nearcorner(dbx,dby,xi,yi);
    return;
  } else if( mindxy<xhalo){
    option=MOVEBOX;
    curb=moveb;
  } else {
    option=NEWBOX;
  }
  pinx=(int)dbx;
  piny=(int)dby;
}

void BoxVec::step2(int x, int y){
  /* (x,y) is the current cursor position */
  switch(option){
  case RESIZEBOX:
    if(curb==NULL) return;
    curb->xc[xi]=x;
    curb->yc[yi]=y;
    break;
  case MOVEBOX:{
    if(curb==NULL) return;
    int dx=x-pinx, dy=y-piny;
    pinx=x;
    piny=y;
    curb->xc[0]+=dx;
    curb->xc[1]+=dx;
    curb->yc[0]+=dy;
    curb->yc[1]+=dy;
    break;
  }
  case NEWBOX:
    newbox(pinx,piny,pinx,piny,NULL);
    option=RESIZEBOX;
    break;
  }
  YFTcl_EvalEx(itp, 256,
	       "$imgWN marker configure %s -coords {%d %d %d %d %d %d %d %d %d %d}",
	       curb->name,
	       curb->xc[0], curb->yc[0],
	       curb->xc[0], curb->yc[1],
	       curb->xc[1], curb->yc[1],
	       curb->xc[1], curb->yc[0],
	       curb->xc[0], curb->yc[0]
	       );
}

void BoxVec::step3(){
  if(curb!=NULL)
  YFTcl_EvalEx(itp, 128,
	       "$imgWN marker configure _%s -coords {%d %d}",
	       curb->name,
	       curb->minxc(), curb->minyc()
	       );
}

void BoxVec::newbox(int x, int y, int x1, int y1, char* nm){
  /***
      pre: given the two diagnal points (x,y) and (x1,y1), and the box name 'nm'
      post: create a box on $imgWN window
  ****/

  int newflag;
  char* name=nm; /* ?? bad operation, will change */
  if(name==NULL){
    int temp;
    name=(char*)malloc(32);
    if(recyclename.size()==0) { temp=nextname; nextname++;}
    else { temp=recyclename[recyclename.size()-1]; recyclename.pop_back();}
    sprintf(name,"b%d",temp);
  }

  curb=(Box*)malloc(sizeof(Box)+strlen(name)) ;
  Tcl_HashEntry *entry=Tcl_CreateHashEntry(&bvhash, name, &newflag);
  strcpy(curb->name, name);
  curb->set(x, y, x1, y1);
  entry->clientData=curb;
  xi=yi=1;
  YFTcl_EvalEx(itp, 256,
	       "eval $imgWN marker create line -name %s $bxlopt -coords \\{%d %d %d %d %d %d %d %d %d %d\\}\neval $imgWN marker create text -name _%s -text %s -coords \\{%d %d\\} $bxtopt",
	       curb->name,
	       curb->xc[0], curb->yc[0],
	       curb->xc[0], curb->yc[1],
	       curb->xc[1], curb->yc[1],
	       curb->xc[1], curb->yc[0],
	       curb->xc[0], curb->yc[0],
	       curb->name,curb->name,
	       curb->minxc(), curb->minyc()
	       );
  if(nm==NULL) free(name);
}

void BoxVec::del(char* nm){
  /***
      delete the box named by 'nm'
  ****/

  /* get the box pointer from the hash table */
  Tcl_HashEntry *entryPtr=Tcl_FindHashEntry(&bvhash, nm);
  Box* bp=(Box*)(entryPtr->clientData);

  /* delete the box from the screen */
  YFTcl_EvalEx(itp, 128,
	       "$imgWN marker delete %s; $imgWN marker delete _%s",
	       bp->name, bp->name);

  /* clean up resources */
  int id=atoi(nm+1);
  recyclename.push_back(id); /* ? better recycling scheme */
  free(entryPtr->clientData);
  Tcl_DeleteHashEntry(entryPtr);
  curb=NULL;
}

Box* BoxVec::boxPtr(char* nm){
  /****
       given the name of the box 'nm',
       return the pointer to the box object
  ****/
  Tcl_HashEntry *entryPtr=Tcl_FindHashEntry(&bvhash, nm);
  if( entryPtr == NULL ) return NULL;
  else return (Box*)(entryPtr->clientData);
}
