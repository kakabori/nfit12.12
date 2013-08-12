#include <stdio.h>
#include <mydll.h>

static int nfunc,ndfunc;

double func(double x[], void* v_1, void* v_2)
{
        double x1p2sqr=sqr(2.0+x[0]);

        nfunc++;
        return 1*(10.0*
                sqr(sqr(x[1])*(3.0-x[0])-sqr(x[0])*(3.0+x[0]))+
                x1p2sqr/(1.0+x1p2sqr));
}

void dfunc(double x[],double df[], void* v_1, void* v_2)
{
        double x1sqr=sqr(x[0]),x2sqr=sqr(x[1]),x1p2=x[0]+2.0;
        double x1p2sqr=sqr(x1p2);

        ndfunc++;
        df[0]=20.0*(x2sqr*(3.0-x[0])-x1sqr*(3.0+x[0]))*(-x2sqr-6.0*x[0]-3.0*x1sqr)+
                2.0*x1p2/(1.0+x1p2sqr)-2.0*x1p2*x1p2sqr/sqr((1.0+x1p2sqr));
        df[1]=40.0*(x2sqr*(3.0-x[0])-x1sqr*(3.0+x[0]))*x[1]*(3.0-x[0]);
}


double fun_(double *x, void *v1, void *v2)
{
    double sum=0;
    int i;
	/* fstar = 0 at xstar = (0,...,0) */
    for (i = 0; i < 10; ++i) {
		sum += -x[i] * x[i]/(i+1)/(i+1);
    }
    return sum-2*x[0]+4*x[2]/9;
} 


int demo()
{
    double work[66];
    int i, n, iflag=-1;
    double scale[10], x[10];
    int nfinc, iwork[15];
    YLSUBPLEXctrl *ctrl = (YLSUBPLEXctrl *)malloc(sizeof(YLSUBPLEXctrl));
    static int mode;
    static double fx, tolfac;
    static int maxnfe, mdsing, mdcont, mduser, nf1, nf2;
    static int nfe;
    static double scl;
    static double tol, tol1, tol2;

/* This program uses subplx to minimize the function fun. */

/* constants */
/*     See subplx comments for storage requirements. */


    printf("********************************************\n");
    printf("******  subplx minimization of fun  ********\n");
    printf("********************************************\n");

/* For descriptions of subplx arguments see subplx comments. */
	n=10;
	n=2;
/* The following two read statements determine when subplx */
/* is interrupted so the user can examine intermediate */
/* results.  subplx can be interrupted and then continued as */
/* if no interrupt had occured when an optimization tolerance */
/* is satisfied and/or when the maximum number of objective */
/* function evaluations is reached. */

/* Variables that define the sequence of tolerances. */
/* If tol = 0 is used, subplx will optimize to the limits */
/* of machine precision. */
/* See subplx comments for description of tol. */

   	tol1=0.1;
	tol2=1e-4;
	tolfac=0.1;

/* Variables that define the sequence of maximum number of */
/* function evaluations. */
/* See subplx comments for description of maxnfe. */

   	nf1=100;
	nf2=1000;
	nfinc=100;

/* Set initial stepsizes for optimization. */
/* See subplx comments for description of scale. */

	scl=0.1;
    scale[0] = -fabs(scl);

/* Set starting point for optimization. */
/* See subplx comments for description of x. */

    for (i = 0; i < n; ++i) {
		x[i]=1;
    }

    printf("n =  %d\n", n);
    printf("tol1,tol2,tolfac = %g %g %g\n",tol1,tol2,tolfac);

    tol = tol1;

    printf("nf1,nf2,nfinc = %d %d %d\n", nf1,nf2,nfinc);

    maxnfe = nf1;

    printf("scale = %g\n", *scale);

    printf("x0 =\n");
    for (i = 0; i < n; ++i) {
      printf("\t%g", x[i]);
    }
    printf("\n");
    /* Print output headers. */
    printf("maxnfe\ttol\tfx\tnfe\tiflag\n");
    printf("%d\t%g\t%g\t%d\t%d\n", maxnfe, tol, fx, nfe, iflag);
    
    /* Set subplx's operating mode. */
    /* See subplx comments for description of mode. */
    
    /* First call to subplx so continuation mode is off. */
    mdcont = 0;
    /* Using default options so user options mode is off. */
    mduser = 0;
    /* Using optimization so single-step mode is off. */
    mdsing = 0;

L20:

    mode = (mdsing << 2) + (mduser << 1) + mdcont;
    
    yl_subplex_subplx(func, n, tol, maxnfe, mode, scale, x, &fx, &nfe, work, iwork, 
		       &iflag, ctrl);
    
    /* Print intermediate results. */
    printf("%d\t%g\t%g\t%d\t%d\n", maxnfe, tol, fx, nfe, iflag);
    
    /* Check iflag to see if done or which termination */
    /* test needs to be reset before resuming optimization. */
    
    if (iflag == -1) {
      if (maxnfe >= nf2) goto Result;
      maxnfe += nfinc;
    } else if (iflag == 0) {
      if (tol <= tol2) goto Result;
      tol *= tolfac;
    } else goto Result;

/* Resume optimization in continuation mode. */
    mdcont = 1;
    goto L20;
Result:
    printf("******  optimization of fun  ********\n");
    printf("iflag = %d, nfe = %d, fx = %g\n",iflag,nfe,fx);
    printf("x =\n");
    for (i = 0; i < n; ++i) printf("\t%g",x[i]);
    printf("\n");
    return 0;
}


int main(void)
{
	YLMinimization mini;
	YLSUBPLEXctrl* subctrl = new YLSUBPLEXctrl;
	DfpminCtrl* dfpctrl = new DfpminCtrl;

	dfpctrl->func = func;
	dfpctrl->dfunc = dfunc;
	dfpctrl->nvar = 2;
	dfpctrl->maxiter = 200;

	mini.setup(dfpctrl, DFPVM);
	mini.x[0] = 0.1;
	mini.x[1] = 4.2;

        printf("True minimum is at (-2.0,+-0.89442719)\n");
        nfunc=ndfunc=0;

        printf("Starting vector: (%7.4f,%7.4f)\n",mini.x[0], mini.x[1]);

	mini.fit();

	//        dfpmin(p, 2, 1.0e-4, &iter, &fret, func, dfunc);
        printf("Iterations: %3d\n", mini.iter);
        printf("Func. evals: %3d\n",nfunc);
        printf("Deriv. evals: %3d\n",ndfunc);
        printf("Solution vector: (%9.6f,%9.6f)\n",mini.x[0],mini.x[1]);
        printf("Func. value at solution %14.6g\n",mini.val);


	subctrl->func = func;
	subctrl->nvar = 2;
	subctrl->nsmin = yfmin(2, subctrl->nvar);
	subctrl->nsmax = yfmin(5, subctrl->nvar);
	mini.setup(subctrl, SUBPLEX);
	mini.x[0] = 0.1;
	mini.x[1] = 4.2;
	mini.fit();
        return 0;
}

