// use DSscale for similar triangle
// use DSfit2 to fit 2 sample given S

#include <cmath>
#include <iostream>
#include <tk.h>
#include <blt.h>

#include "toadcmd.h"

using namespace std;

extern Blt_Vector *xVP[4], *plotxVP, *plotyVP;
extern double tcl_chi;
extern SetUp setup;
/* the following table shows what each x means
   x[0]	s
   x[1]	x0
   x[2]	alpha
   x[3]	d
   x[4]	R
   x[5]	lambda
   x[6] pz
   x[7] delts
*/
extern double ln[32],step[4];
extern int boolx[32],xpin[4];
extern double ln1[32],ln0[32];
extern double sx[3];
extern int sxpin[4];


static double getl(double theta,double *x){
  return (x[4]*cos(theta)+tan(2*theta)*(x[0]+x[4]*sin(theta)))/(cos(x[2])-tan(2*theta)*sin(x[2]))+x[1];
}
static double fit2chi(double *x){
  double temp,chi=0,*setupx=(double *)&setup,thetai;
  for(int i=0;i<32;i++){
	if(boolx[i]){
	thetai=asin(i*setupx[5]/2/setupx[3]);
        //thetai = acos(0.999998 * cos(thetai) );
          thetai = acos(setup.refn*cos(thetai) );//changed by NC for ref. index 3.8.04
	if(i==0) thetai = 0;
	} else continue;
    temp=(setupx[4]*cos(thetai)+tan(2*thetai)*(x[0]+setupx[4]*sin(thetai)))/(cos(setupx[2])-tan(2*thetai)*sin(setupx[2]))-ln0[i]+x[1];
    chi+=temp*temp;
    temp=(setupx[4]*cos(thetai)+tan(2*thetai)*(x[0]+setupx[7]+setupx[4]*sin(thetai)))/(cos(setupx[2])-tan(2*thetai)*sin(setupx[2]))-ln1[i]+x[2];
    chi+=temp*temp;
    //		cout<<i<<endl;
  }
  return chi;
}
static double get_theta(double L,double S, int n){
  double R=setup.R;
  double x0=0,Lstar;
  int i;
  for(i=0;i<n;i++){ // n=20 will be enough to get precise value of x0
    Lstar=L/sqrt(1-x0*x0);
    x0=(sqrt(S*S+2*(Lstar-R)*Lstar)-S)/2/Lstar;
  }
  return asin(x0);
}
static double scalechi(double *x){
  double temp,chi=0,*setupx=(double *)&setup,thetai;
  for(int i=0;i<32;i++){
    if(boolx[i])  thetai=get_theta(ln0[i]-x[1],x[0],20); else continue;
    //		temp=(setupx[4]*cos(thetai)+tan(2*thetai)*(x[0]+setupx[4]*sin(thetai)))/(cos(setupx[2])-tan(2*thetai)*sin(setupx[2]))-ln0[i]+x[1];
    //		chi+=temp*temp;
    temp=(setupx[4]*cos(thetai)+tan(2*thetai)*(x[0]+setupx[7]+setupx[4]*sin(thetai)))/(cos(setupx[2])-tan(2*thetai)*sin(setupx[2]))-ln1[i]+x[2];
    chi+=temp*temp;
  }
  return chi;
}
static double chis(double *x){
  double chi=0,temp, thetai;
  int i;
  for(i=0;i<32;i++){
	if(boolx[i]){
		thetai=asin(i*x[5]/2/x[3]);
	        //thetai = acos(0.999998 * cos(thetai) );
        	  thetai = acos(setup.refn*cos(thetai) );//changed by NC for ref. index 3.8.04
		if(i==0) thetai = 0;
	} else continue;
    temp=(x[4]*cos(thetai)+tan(2*thetai)*(x[0]+x[4]*sin(thetai)))/(cos(x[2])-tan(2*thetai)*sin(x[2]))-ln[i]+x[1];
    chi+=temp*temp;
  }
  //	cout<<chi<<endl;
  return chi;
}

static double findmin(double *x,double *step, int *xpin,int varnum, int iter){
  int j,i,k,adjust,limit,limit2;
  double abr,err;
  err=chis(x);
  //	cout<<x[0]<<' '<<x[1]<<' '<<x[2]<<' '<<x[3]<<endl;
  //	cout<<xpin[0]<<' '<<xpin[1]<<' '<<xpin[2]<<' '<<xpin[3]<<endl;
  //	cout<<step[0]<<' '<<step[1]<<' '<<step[2]<<' '<<step[3]<<endl;
  for(j=0,limit=0;j<iter;j++){
    do{
      adjust=0;
      for(i=0;i<varnum;i++){
	if(!xpin[i]) continue;
	x[i]+=step[i];
	if((abr=chis(x))<err) {err=abr;adjust++;}
	else{x[i]-=step[i];step[i]=-step[i];}
	limit2=0;
	do{
	  x[i]+=step[i];
	  if((abr=chis(x))<err) {err=abr;adjust++;}
	  else {
	    x[i]-=step[i];
	    break;
	  }
	  if(limit2++>6000) {tcl_chi=err;return(err);}
	}while(1);
      }
      if(limit++>100) {tcl_chi=err;return(err);}
    }while(adjust>3);
    for(k=0;k<varnum;k++) step[k]/=3.;
  }
  //	cout<<limit<<' '<<err<<endl;
  tcl_chi=err;
  return(err);
}

typedef double (*CHIFunc)(double *);

static double findminG(double *x,double *step, int *xpin,int varnum, CHIFunc chifunc, int iter){
  int j,i,k,adjust,limit,limit2;
  double abr,err;
  err=(*chifunc)(x);
  for(j=0,limit=0;j<iter;j++){
    do{
      adjust=0;
      for(i=0;i<varnum;i++){
	if(!xpin[i]) continue;
	x[i]+=step[i];
	if((abr=(*chifunc)(x))<err) {err=abr;adjust++;}
	else{x[i]-=step[i];step[i]=-step[i];}
	limit2=0;
	do{
	  x[i]+=step[i];
	  if((abr=(*chifunc)(x))<err) {err=abr;adjust++;}
	  else {
	    x[i]-=step[i];
	    break;
	  }
	  if(limit2++>6000) {return(err);}
	}while(1);
      }
      if(limit++>100) {return(err);}
    }while(adjust>3);
    for(k=0;k<varnum;k++) step[k]/=3.;
  }
  return(err);
}

void UpdateVarDS(Tcl_Interp *interp){
  char varname[64];
  sprintf(varname,"S");
  Tcl_UpdateLinkedVar(interp,varname);
  sprintf(varname,"x0");
  Tcl_UpdateLinkedVar(interp,varname);
  sprintf(varname,"D");
  Tcl_UpdateLinkedVar(interp,varname);
  sprintf(varname,"alpha");
  Tcl_UpdateLinkedVar(interp,varname);
  sprintf(varname,"chi");
  Tcl_UpdateLinkedVar(interp,varname);
}

void UpdateVPDS(int spnum, SetUp *setupP){
  int i,k=0;
  double *x=plotxVP->valueArr,*y=plotyVP->valueArr;
  double *expt=xVP[0]->valueArr,*theo=xVP[1]->valueArr;
  if(spnum==2) expt=xVP[2]->valueArr,theo=xVP[3]->valueArr;
  for(i=0;i<32;i++){
	if(i==0) theo[i]=getl(0,(double *)setupP)/setup.pz;
	//else     theo[i] = getl(acos(0.999998*cos(asin(i*setup.lambda/2/setup.d))),(double *)setupP)/setup.pz;
        else     theo[i] = getl(acos(setup.refn*cos(asin(i*setup.lambda/2/setup.d))),(double *)setupP)/setup.pz;
    	//changed by NC for ref. index 3.8.04
    if(boolx[i]) {
      x[k]=i;
      y[k]=(expt[i]-theo[i]);
      k++;
    }
  }
  Blt_ResetVector(plotxVP,x, k, 32, TCL_DYNAMIC);
  Blt_ResetVector(plotyVP,y, k, 32, TCL_DYNAMIC);
  //	Blt_ResetVector(xVP[1],theo, 32, 32, TCL_STATIC);setup.pz**setup.pz
}
int YF_DSfitting(ClientData clientData, Tcl_Interp *interp,int objc, Tcl_Obj *const objv[]){
  double *x=xVP[0]->valueArr; //defalt: use xVP[0] for x, spnum updates it
  //	double *step=xVP[3]->valueArr;
  int iter,spnum;
  Tcl_GetIntFromObj(interp, objv[1], &iter);
  Tcl_GetIntFromObj(interp, objv[2], &spnum);
  if(spnum==2) x=xVP[2]->valueArr;
  for(int i=0;i<32;i++) ln[i]=setup.pz*x[i];
  findmin((double *)&setup,step,xpin,4,iter);
  UpdateVarDS(interp);
  UpdateVPDS(spnum,&setup);
  return TCL_OK;
}

int YF_DSdev(ClientData clientData, Tcl_Interp *interp,int objc, Tcl_Obj *const objv[]){
  int spnum;
  SetUp tsetup=setup;
  Tcl_GetIntFromObj(interp, objv[1], &spnum);
  if(spnum==2){
    tsetup.s=sx[0]+setup.delts;
    tsetup.x0=sx[2];
  }else{
    tsetup.s=sx[0];
    tsetup.x0=sx[1];
  }
  UpdateVPDS(spnum,&tsetup);
  return TCL_OK;
}
int YF_DSfit2(ClientData clientData, Tcl_Interp *interp,int objc, Tcl_Obj *const objv[]){
  int iter;
  Tcl_GetIntFromObj(interp, objv[1], &iter);
  if(objc>2) Tcl_GetDoubleFromObj(interp, objv[2], sx);
  if(objc>3) Tcl_GetDoubleFromObj(interp, objv[3], sx+1);
  if(objc>4) Tcl_GetDoubleFromObj(interp, objv[4], sx+2);
  for(int i=0;i<32;i++){
    ln0[i]=setup.pz*xVP[0]->valueArr[i];
    ln1[i]=setup.pz*xVP[2]->valueArr[i];
  }
  step[0]=step[1]=step[2]=1;
  if(clientData){
    cout<<findminG(sx,step,sxpin,3,fit2chi, iter)<<endl;
  }else{
    cout<<findminG(sx,step,sxpin,3,scalechi, iter)<<endl;
  }
  cout<<sx[0]<<' '<<sx[1]<<' '<<sx[2]<<endl;
  return TCL_OK;
}

int Init_DSvectors(Tcl_Interp *interp)
{
  char vectorname[8];
  int i;
  sprintf(vectorname,"plotxDS");
  if (Blt_CreateVector(interp, vectorname, 32, &plotxVP) != TCL_OK) return TCL_ERROR;
  sprintf(vectorname,"plotyDS");
  if (Blt_CreateVector(interp, vectorname, 32, &plotyVP) != TCL_OK) return TCL_ERROR;
  for (i=0; i<4; i++)
    {
      sprintf(vectorname,"xv%d",i);
      if (Blt_CreateVector(interp, vectorname, 32, &xVP[i]) != TCL_OK) return TCL_ERROR;
    }
  return TCL_OK;
}
