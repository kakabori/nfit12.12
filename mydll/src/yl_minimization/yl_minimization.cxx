#include <mydll.h>

DfpminCtrl::DfpminCtrl(){
    eps = 3.0e-8;
    tolx = 1e-7;
    stpmx = 100;
    gtol = 1e-7;
    alf = 1e-4;
    maxiter = 200;
}


YLSUBPLEXctrl::YLSUBPLEXctrl(){
    alpha = 1.;
    beta = .5;
    gamma = 2.;
    delta = .5;
    psi = .25;
    omega = .1;
    irepl = 0;
    ifxsw = 1;
    bonus = 1.;
    nfstop = 0;
    tol = 1e-7;
    maxnfe = 1000;
}


void dfpmin(double* p, int n,int *iter, double *fret,
	    YLMinimizationFunc func,
	    YLMinimizationDFunc dfunc, void* client, DfpminCtrl *ctrl);


int YLMinimization::fit(){
  switch(method){
  case DFPVM:{
    DfpminCtrl *dfpctrl = (DfpminCtrl*)ctrl;
    iter = dfpctrl->maxiter;
    dfpmin(x,
	   dfpctrl->nvar,
	   &iter,
	   &val,
	   dfpctrl->func,
	   dfpctrl->dfunc,
	   client,
	   dfpctrl);
    break;}
  case SUBPLEX:{
    double* work;
    double *scale;
    int iflag=-1;
    int *iwork;
    int mdcont=0, totnfe=100 ;
    double scl=0.1;
    double tol=0.1;

    YLSUBPLEXctrl *subctrl = (YLSUBPLEXctrl*)ctrl;
    n = subctrl->nvar;
    scale = (double*)calloc(n, sizeof(double) );
    work  = (double*)calloc(2*n + subctrl->nsmax*(subctrl->nsmax+4) + 1, sizeof(double) );
    iwork = (int*)calloc( n + (int)(n/subctrl->nsmin), sizeof(int) );
    scale[0] = -fabs(scl);

    do{
      yl_subplex_subplx(subctrl->func,
			 subctrl->nvar,
			 tol,
			 totnfe,
			 mdcont,
			 scale,
			 x,
			 &val,
			 &iter,
			 work,
			 iwork, 
			 &iflag,
			 ctrl);
      printf("%d\t%g\t%g\t%d\t%d\n", totnfe, tol, val, iter, iflag);
      if (iflag == -1) {
	if (totnfe >= subctrl->maxnfe) goto Result;
	totnfe += 100;
      } else if (iflag == 0) {
	if (tol <= subctrl->tol) goto Result;
	tol *= 0.1;
      } else goto Result;
         mdcont = 1;
    } while(1);
  Result:
    printf("******  optimization of fun  ********\n");
    printf("iflag = %d, nfe = %d, fx = %g\n",iflag,iter,val);
    printf("x =\n");
    for (int i = 0; i < subctrl->nvar; ++i) printf("\t%g",x[i]);
    printf("\n");
    free(scale);
    free(work);
    free(iwork);
    break;
  }
  default: ;
  }
  return 0;
}


int YLMinimization::setup(void *_ctrl, int _meth){
  ctrl = _ctrl;
  method = _meth;
  switch(method){
  case DFPVM:{
    DfpminCtrl *dfpctrl = (DfpminCtrl*)ctrl;
    tvecr(&x, dfpctrl->nvar);
    break;}
  case SUBPLEX:{
    YLSUBPLEXctrl *subctrl = (YLSUBPLEXctrl*)ctrl;
    tvecr(&x, subctrl->nvar);
    break;
  }
  default: ;
  }
  return 0;
}
