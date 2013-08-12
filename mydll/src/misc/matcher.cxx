#include "mydll.h"

double errorWrapper(double *a, void *v1, void *v2 ){
  Matcher *matcher = (Matcher *)v1;
  return matcher->error( a[2], a[0], a[1] );
}

void derrorWrapper(double *a, double *df, void *v1, void *v2 ){
  Matcher *matcher = (Matcher *)v1;
  matcher->derror( a[2], a[0], a[1] , df );
}


Matcher::Matcher(int maxm){ setsize(maxm); m=0;}

void  Matcher::setsize(int maxm){ aut.cautiousResize(4, maxm); b.cautiousResize(1, maxm); }

void Matcher::reset(){ m=0; }

void Matcher::push(double rx, double ry, double x, double y){
  aut.value(0, m) = x;
  aut.value(0, m+1) = y;
  aut.value(1, m) = -y;
  aut.value(1, m+1) = x;  
  aut.value(2, m) = 1;
  aut.value(2, m+1) = 0;
  aut.value(3, m) = 0;
  aut.value(3, m+1) = 1;

  b.rvalue(m) = rx;
  b.rvalue(m+1) = ry;
  m += 2;
}

double Matcher::error(double theta, double dx, double dy){
  cmatrix tx(1,4);
  cmatrix bufa;

  tx.rvalue(0) = cos(theta);
  tx.rvalue(1) = sin(theta);
  tx.rvalue(2) = dx;
  tx.rvalue(3) = dy;
  bufa.gemm( tx, aut.view(0, 4, 0, m) );
  bufa -= b.view(0,1,0,m);
  return bufa.norm22();
}

void Matcher::derror(double theta, double dx, double dy , double *df){
  cmatrix tx(1,4);
  cmatrix bufa;
  double ct = cos(theta);
  double st = sin(theta);

  tx.rvalue(0) = ct;
  tx.rvalue(1) = st;
  tx.rvalue(2) = dx;
  tx.rvalue(3) = dy;
  bufa.gemm( tx, aut.view(0, 4, 0, m) );
  bufa -= b.view(0,1,0,m);


  df[0]=df[1]=df[2]=0;

  for(int i=0; i<m; i+=2){
    df[0] += bufa.rvalue(i);
    df[1] += bufa.rvalue(i+1);
    df[2] += bufa.rvalue(i) * ( -aut.value(0,i) * st - aut.value(0,i+1) * ct ) + 
             bufa.rvalue(i+1) * ( aut.value(0,i)  * ct - aut.value(0,i+1) * st );
  }

  df[0] *= 2;
  df[1] *= 2;
  df[2] *= 2;
  //cout<<df[0]<< ' '<<df[1]<< ' '<<df[2]<< ' '<<endl;
}


void Matcher::svdMatch(cmatrix &x){
  autV = aut.view(0, 4, 0, m); 
  autV.tsvd_nr( w, vt );

  //  cout<<autV<<endl <<w << endl <<vt <<endl<< b.view(0,1,0,m) <<endl;
  autV.svbksb( w, vt, b.view(0,1,0,m), x );
  
  //cout<<x<<endl;
 
  //cmatrix bufa;
  //bufa.gemm(x, aut.view(0, 4, 0, m));
  //bufa -= b.view(0,1,0,m);
  // cout<<"bufa"<< bufa <<endl;
  //cout<<aut.view(0, 4, 0, m)<<endl;
  //cout<<b.view(0,1,0,m)<<endl;
}

double Matcher::dfpMatch(cmatrix &x){
  double theta = atan2( x.rvalue(1) , x.rvalue(0) );

  YLMinimization mini;
  DfpminCtrl* dfpctrl = new DfpminCtrl;

  dfpctrl->func = errorWrapper;
  dfpctrl->dfunc = derrorWrapper;
  dfpctrl->nvar = 3;
  dfpctrl->maxiter = 20;
  //*
  dfpctrl->eps = 0;
  dfpctrl->tolx = 0;
  dfpctrl->gtol = 0;
  dfpctrl->alf = 0;
  //*/
  mini.client = this;
  mini.setup(dfpctrl, DFPVM);
  mini.x[0] = x.rvalue(2);
  mini.x[1] = x.rvalue(3);
  mini.x[2] = theta;

 
  mini.fit();
 
  x.rvalue(0) = cos(mini.x[2]);
  x.rvalue(1) = sin(mini.x[2]);
  x.rvalue(2) = mini.x[0];
  x.rvalue(3) = mini.x[1];
  double e = error(mini.x[2], mini.x[0], mini.x[1]);
  delete dfpctrl;
  cout<<e<<endl;
  return e;
}


double Matcher::match(cmatrix &x){
  svdMatch(x);
  return dfpMatch(x);
}

void Matcher::match2(cmatrix &x, double thresh, void**seq){
  match(x);
  removeOutlier(x, thresh, seq);
  match(x);
}

void Matcher::removeOutlier(cmatrix &x, double thresh, void**seq){
  cmatrix bufa;
  bufa.gemm( x, aut.view(0, 4, 0, m) );
  bufa -= b.view(0,1,0,m);

  for(int i=0; i<m; i+=2){
    if( sqrt(sqr(bufa.rvalue(i)) + sqr(bufa.rvalue(i+1))) > thresh ){
      m -= 2;

      if(seq) tswap(seq[i/2], seq[m/2]);

      tswap( aut.value(0,i) , aut.value(0, m) );
      tswap( aut.value(0,i+1) , aut.value(0, m+1) );

      tswap( aut.value(1,i) , aut.value(1, m) );
      tswap( aut.value(1,i+1) , aut.value(1, m+1) );

      tswap( b.value(0,i) , b.value(0, m) );
      tswap( b.value(0,i+1) , b.value(0, m+1) );

      tswap( bufa.value(0,i) , bufa.value(0, m) );
      tswap( bufa.value(0,i+1) , bufa.value(0, m+1) );

      i -= 2;
    }
  }
}
void Matcher::match2(cmatrix &x, double thresh, int *seq){
  match(x);
  removeOutlier(x, thresh, seq);
  match(x);
}

void Matcher::removeOutlier(cmatrix &x, double thresh, int *seq){
  cmatrix bufa;
  bufa.gemm( x, aut.view(0, 4, 0, m) );
  bufa -= b.view(0,1,0,m);

  for(int i=0; i<m; i+=2){
    if( sqrt(sqr(bufa.rvalue(i)) + sqr(bufa.rvalue(i+1))) > thresh ){
      m -= 2;

      if(seq) tswap(seq[i/2], seq[m/2]);

      tswap( aut.value(0,i) , aut.value(0, m) );
      tswap( aut.value(0,i+1) , aut.value(0, m+1) );

      tswap( aut.value(1,i) , aut.value(1, m) );
      tswap( aut.value(1,i+1) , aut.value(1, m+1) );

      tswap( b.value(0,i) , b.value(0, m) );
      tswap( b.value(0,i+1) , b.value(0, m+1) );

      tswap( bufa.value(0,i) , bufa.value(0, m) );
      tswap( bufa.value(0,i+1) , bufa.value(0, m+1) );

      i -= 2;
    }
  }
}
