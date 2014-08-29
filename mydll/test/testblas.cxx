#include <mydll.h>
#include <iostream>
#include <time.h>

#define DIMI 503
#define DIMK 503
#define DIMJ 503

void calcsum(int n, double* c, double* c2){
  double  sum = 0;
  for(int i=0; i<DIMI; i++) sum+=fabs(c[i]-c2[i]);
  cout<<"sum: "<<sum<<endl;
}

int testblock(){
  double *a, *b, *c, *c2;
  double sum;
  clock_t t1, t2;
  a = (double*)malloc(DIMI*DIMK*sizeof(double));
  b = (double*)malloc(DIMK*DIMJ*sizeof(double));
  c = (double*)malloc(DIMI*DIMJ*sizeof(double));
  c2 = (double*)malloc(DIMI*DIMJ*sizeof(double));
		      
  for(int i=0; i<DIMI*DIMK; i++) a[i] = (rand()%1000);
  for(int i=0; i<DIMK*DIMJ; i++) b[i] = (rand()%1000);

  //yl_blas_dset(DIMI*DIMK, 10001, a, 1);
  //yl_blas_dset(DIMK*DIMJ, 10002, b, 1);

  yl_blas_dset(DIMI*DIMJ, 0, c, 1);
  yl_blas_dset(DIMI*DIMJ, 0, c2, 1);
  t1 = clock();
  yl_cblas_dgemm2(a, DIMK, 1,
		b, DIMJ, 1,
		c2, DIMJ, 1,
		DIMI, DIMK, DIMJ);
  t2 = clock();
  cout<<"b_ijk: "<<(t2-t1)/CLK_TCK<<endl;

  //return 0;

  t1 = clock();
  yl_cblas_dgemm_auto(a, DIMK, 1,
	       b, DIMJ, 1,
	       c, DIMJ, 1,
	       DIMI, DIMK, DIMJ, 1, 0, 0);
  t2 = clock();
  cout<<"bp_ijk: "<<(t2-t1)/CLK_TCK<<endl;
  sum = 0;
  for(int i=0; i<DIMI*DIMJ; i++) sum+=fabs(c[i]-c2[i]);
  cout<<"sum: "<<sum<<endl;
  yl_blas_dset(DIMI*DIMJ, 0, c, 1);
 
  t1 = clock();
  yl_cblas_dgemm_n_ijk(a, DIMK, 1,
		       b, DIMJ, 1,
		       c, DIMJ, 1,
		       DIMI, DIMK, DIMJ,0);
  t2 = clock();
  cout<<"n_ijk: "<<(t2-t1)/CLK_TCK<<endl;
  sum = 0;
  for(int i=0; i<DIMI*DIMJ; i++) sum+=fabs(c[i]-c2[i]);
  cout<<"sum: "<<sum<<endl;
  yl_blas_dset(DIMI*DIMJ, 0, c, 1);
 
  t1 = clock();
  yl_cblas_dgemm_np_ijk(a, DIMK, 1,
		       b, DIMJ, 1,
		       c, DIMJ, 1,
		       DIMI, DIMK, DIMJ);
  t2 = clock();
  cout<<"np_ijk: "<<(t2-t1)/CLK_TCK<<endl;
  sum = 0;
  for(int i=0; i<DIMI*DIMJ; i++) sum+=fabs(c[i]-c2[i]);
  cout<<"sum: "<<sum<<endl;
  
  return 0;
}


int testblock2(char opt, int f){
  double *a, *b, *c, *c2;
  double sum;
  clock_t t1, t2;
  a = (double*)malloc(DIMI*DIMK*sizeof(double));
  b = (double*)malloc(DIMK*DIMJ*sizeof(double));
  c = (double*)malloc(DIMI*DIMJ*sizeof(double));
  c2 = (double*)malloc(DIMI*DIMJ*sizeof(double));
  
  for(int i=0; i<DIMI*DIMK; i++) a[i] = (rand()%1000);
  for(int i=0; i<DIMK*DIMJ; i++) b[i] = (rand()%1000);
  
  
  yl_blas_dset(DIMI*DIMJ, 0, c, 1);
  yl_blas_dset(DIMI*DIMJ, 0, c2, 1);
  
  t1 = clock();
  /* yl_blas_dgemm_np_jik_opt(a, DIMI, 1,
			b, DIMK, 1,
			c, DIMI, 1,
			DIMI, DIMK, DIMJ, opt);*/
  yl_blas_dgemm_auto(a, DIMI, 1,
		     b, DIMK, 1,
		     c, DIMI, 1,
		     DIMI, DIMK, DIMJ, 1, 0, opt);//*/
  t2 = clock();
  cout<<"bp_ijk: "<<(t2-t1)/CLK_TCK<<endl;
  
  if(f) return 0;
  t1 = clock();
  /*
  yl_blas_dgemm_n_ijk(a, DIMI, 1,
		      b, DIMK, 1,
		      c2, DIMI, 1,
		      DIMI, DIMK, DIMJ, opt);*/
  yl_blas_dgemm_np_jik_opt(a, DIMI, 1,
			   b, DIMK, 1,
			   c2, DIMI, 1,
			   DIMI, DIMK, DIMJ, opt);//*/

  t2 = clock();
  cout<<"n_ijk: "<<(t2-t1)/CLK_TCK<<endl;
  sum = 0;
  for(int i=0; i<DIMI*DIMJ; i++) sum+=fabs(c[i]-c2[i]);
  cout<<"sum: "<<sum<<endl;
  free(a); free(b); free(c); free(c2);
  return 0;
}

int testblock3(char opt, int f){
  double *a, *b, *c, *c2;
  double sum;
  clock_t t1, t2;
  a = (double*)malloc(DIMI*DIMK*sizeof(double));
  b = (double*)malloc(DIMK*DIMJ*sizeof(double));
  c = (double*)malloc(DIMI*DIMJ*sizeof(double));
  c2 = (double*)malloc(DIMI*DIMJ*sizeof(double));
  
  for(int i=0; i<DIMI*DIMK; i++) a[i] = (rand()%1000);
  for(int i=0; i<DIMK*DIMJ; i++) b[i] = (rand()%1000);
  
  
  yl_blas_dset(DIMI*DIMJ, 0, c, 1);
  yl_blas_dset(DIMI*DIMJ, 0, c2, 1);
  
  t1 = clock();

  yl_cblas_dgemm_auto(a, DIMK, 1,
		      b, DIMJ, 1,
		      c, DIMJ, 1,
		      DIMI, DIMK, DIMJ, 1, 0, opt);
  t2 = clock();
  cout<<"bp_ijk: "<<(t2-t1)/CLK_TCK<<endl;
  
  if(f) return 0;
  t1 = clock();
   yl_cblas_dgemm_n_ijk(a, DIMK, 1,
		       b, DIMJ, 1,
		       c2, DIMJ, 1,
		       DIMI, DIMK, DIMJ, opt);
  
  t2 = clock();
  cout<<"n_ijk: "<<(t2-t1)/CLK_TCK<<endl;
  sum = 0;
  for(int i=0; i<DIMI*DIMJ; i++) sum+=fabs(c[i]-c2[i]);
  cout<<"sum: "<<sum<<endl;
  
  return 0;
}

void simpledisp(double *c){
  cout<<c[0]<<' '<<c[1]<<' '<<c[2]<<' '<<c[3]<<endl;
  cout<<c[4]<<' '<<c[5]<<' '<<c[6]<<' '<<c[7]<<endl;
  cout<<c[8]<<' '<<c[9]<<' '<<c[10]<<' '<<c[11]<<endl;

}

int testsimple(){
  double a[] = { 1, 2, 3, -1,
		 4, 5, 6, 0,
		 7, 8, 9, 1};
  double b[] = { 0, 1, 0,
		 1, 0, 0};
  double c[12];
  yl_blas_dset(12, 0, c, 1);
  yl_blas_dgemm_n_ijk(a, 4, 1,
		      b, 3, 1,
		      c, 4, 1,
		      4, 3, 3,0);
  simpledisp(c);
  yl_blas_dset(12, 0, c, 1);
  yl_blas_dgemm_np_jik(a, 4, 1,
		      b, 3, 1,
		      c, 4, 1,
		      4, 3, 2);
  simpledisp(c);

  return 0;
}


#define DIMI 1003
#define DIMK 1003


int testdgemv(char opt, int f){
  double *a, *b, *c, *c2;
  double sum;
  clock_t t1, t2;
  a = (double*)malloc(DIMI*DIMK*sizeof(double));
  b = (double*)malloc(DIMK*sizeof(double));
  c = (double*)malloc(DIMI*sizeof(double));
  c2 = (double*)malloc(DIMI*sizeof(double));
  
  for(int i=0; i<DIMI*DIMK; i++) a[i] = (rand()%1000);
  for(int i=0; i<DIMK; i++) b[i] = (rand()%1000);
  
  
  yl_blas_dset(DIMI, 0, c, 1);
  yl_blas_dset(DIMI, 0, c2, 1);
  
  t1 = clock();
  for(int iter=0; iter<100; iter++){
    yl_blas_dgemv_auto(a, DIMI, 1,
		       b, 1,
		       c, 1,
		       DIMI, DIMK, 0.5, 2, opt);
  }
  t2 = clock();
  cout<<"bp_ijk: "<<(t2-t1)/CLK_TCK<<endl;
  
  if(f) return 0;
  t1 = clock();

  for(int iter=0; iter<100; iter++){
    yl_blas_dgemm_auto(a, DIMI, 1,
		       b, 1, 1,
		       c2, 1, 1,
		       DIMI, DIMK, 1, 0.5, 2, opt);
  }
  t2 = clock();
  cout<<"np_ijk: "<<(t2-t1)/CLK_TCK<<endl;
  calcsum(DIMI, c, c2);
  free(a); free(b); free(c); free(c2);
  return 0;
}
















int main(int argc, char**argv){
  //  testsimple();

  int f= atoi(argv[1]);
  /*  
  testblock2(0,f);
  testblock2(1,f);
  testblock2(2,f);
  testblock2(3,f);

  //return 0;
  cout<<"******"<<endl;
  f=1;
  testblock3(0,f);
  testblock3(1,f);
  testblock3(2,f);
  testblock3(3,f);
  //testsimple();*/
  testdgemv(0,f);
  testdgemv(2,f);
  return 0;
}
