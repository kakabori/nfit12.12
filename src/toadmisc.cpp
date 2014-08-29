#include "toadmisc.h"


void YFTcl_SetResult(Tcl_Interp* interp, size_t sz, const char*fmt, ...){
  /* helps to set tcl result */
  va_list ap;
  va_start(ap, fmt);
  char* cmd=(char*)malloc(sz);
  vsprintf(cmd, fmt, ap);
  Tcl_SetResult(interp, cmd, TCL_DYNAMIC);
  va_end(ap);
}

/* helps to execute tcl command */
void YFTcl_EvalEx(Tcl_Interp* interp, size_t sz, const char*fmt, ...)
{
  va_list ap;
  va_start(ap, fmt);
  char* cmd=(char*)malloc(sz);
  vsprintf(cmd, fmt, ap);
  Tcl_EvalEx(interp, cmd, -1, TCL_EVAL_GLOBAL);
  free(cmd);
  va_end(ap);
}

double polint4t(double xmx0,double *ya){
  return (ya[0]*(3-xmx0)+ya[3]*xmx0)*(xmx0-1)*(xmx0-2)/6+
    (ya[2]*(1-xmx0)+ya[1]*(xmx0-2))*(xmx0-3)*xmx0/2;
}
double polint3_log(double xmxa,float *ya){
  return xmxa*50*((100*xmxa-2)*(ya[0]+ya[2]-2*ya[1])+ya[2]-ya[0])+ya[0];
}
double polint30_log(double x,float *ya){
  return (0.9772372212*ya[0]*(x-1)-42.93136751*ya[1]*x)*(x-1.023292992)+41.95413029*ya[2]*(x-1)*x;
}

void polint(const double *xa, const double *ya,int n,double x,double *y,double *dy){
  /* ploynomial interpolation code from numerical recipe */
  int i,m,ns=0;
  double den,dif,dift,ho,hp,w;
  double *c=(double *)malloc(n*sizeof(double)),*d=(double *)malloc(n*sizeof(double));
  dif=fabs(x-xa[0]);
  for (i=0;i<n;i++) {
    if ( (dift=fabs(x-xa[i])) < dif) {
      ns=i;
      dif=dift;
    }
    c[i]=ya[i];
    d[i]=ya[i];
  }
  *y=ya[ns];
  for (m=1;m<n;m++) {
    for (i=0;i<n-m;i++) {
      ho=xa[i]-x;
      hp=xa[i+m]-x;
      w=c[i+1]-d[i];
      if ( (den=ho-hp) == 0.0){
	printf("Error in routine POLINT");
      }
      den=w/den;
      d[i]=hp*den;
      c[i]=ho*den;
    }
    *y += (*dy=(2*ns<(n-m) ? c[ns] : d[--ns]));
  }
  free(d);
  free(c);
}

/****************************************************************************************
Tcl_HashEntry * Tcl_FindHashEntry(Tcl_HashTable * table, const char* key);
locates the entry corresponding to a particular key
returns NULL if the key doesn't exit
Tcl_HashEntry is a pointer to hash table entry

****************************************************************************************/

/* help to find object in a hash table */
Tcl_HashEntry* _getExistingPtr(Tcl_Interp *interp, Tcl_HashTable *htptr, Tcl_Obj *obj, void** ptrptr)
{
  return _getExistingPtr(interp, htptr, Tcl_GetString(obj), ptrptr);
}

Tcl_HashEntry* _getExistingPtr(Tcl_Interp *interp, Tcl_HashTable *htptr, char* name, void** ptrptr)
{
  Tcl_HashEntry *entryPtr = Tcl_FindHashEntry(htptr, name);
  if(entryPtr != NULL){
    *ptrptr = entryPtr->clientData;
  } else {
    *ptrptr = NULL;
  }
  return entryPtr;
}
