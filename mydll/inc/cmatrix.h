#ifndef __CMATRIX_HH
#define __CMATRIX_HH

//#include "mydll.h"



struct IdxArr
{
  int *ia; //int array for storing the indices
  int ia_size; //number of indices stored in IdxArr
  int memsize; //allocated size of ia

  IdxArr(int msz=32);
  IdxArr(const IdxArr &in); //copy constructor
  ~IdxArr();

  inline int operator[](int i) const; //get index #i (0 to n-1)
  inline IdxArr& operator=(const IdxArr&);

  int size() const{ return ia_size; };
  void clear() { ia_size = 0; };

  void addIdx(int x); //add an index to the end(no repeated entries)
  bool exist(int) const;
  void operator+=(int);
  IdxArr& append(const IdxArr &in);
  void hack(int nn){ ia_size = nn; }

private:
  void _increaseMemSize(); //helper function used by addIdx
  void _copy(const IdxArr &in);

};

ostream& operator<<(ostream& out, const IdxArr& idxarr);

//--------------------------------------
// overload operators 

inline int IdxArr::operator[](int i) const {
  assert(i < ia_size);
  return ia[i];
}

inline IdxArr& IdxArr::operator=(const IdxArr& in){
  _copy(in);
  return *this;
}






enum CMAT_METHOD { CMAT_COPY, CMAT_VIEW};
typedef double (*dFunc0i)();
typedef double (*dFunc2i)(int, int);
typedef double (*dFunc1d)(double);
typedef double (*dFunc1d2i)(double,int,int);


class cmatrix {
  double *d;
  int ld;
  int inc;
  int nr;
  int nc;
  int memsize;

  //friend cmatrix& gemm(double alpha, const cmatrix& a , const cmatrix& b , double beta, cmatrix& c);
  friend ostream& operator<<( ostream& out, const cmatrix&);
 public:
  void memdebug() const;
  bool isOwner() const {return memsize>0; }
  int rows() const { return nr; }
  int cols() const { return nc; }
  double * dptr() { return d;}
  cmatrix(){memsize=-1; nr=nc=0; d = NULL;}
  // copy constructor
  cmatrix(const cmatrix& in, enum CMAT_METHOD method );//= CMAT_COPY

  // it maybe a user or owner depending __nr, __nc
  cmatrix(int, int);

  ~cmatrix();

  // for owner type, simple realloc
  void memDirectRealloc(int);
  // mem (de)alloc by given parameter, compact form out
  void destructiveLazyResize(int, int);
  // mem (re/de/m)alloced with given parameter, compact form out
  void destructiveResize(int, int);

  void cautiousResize(int __nr, int __nc);

  void release(){ memsize = -1; }
  // compact owner or assert dimension match
  cmatrix& operator=(const cmatrix& in);

  double asum();
  double sum();
  double amax();
  double norm22();

  cmatrix& resize1(int, int);

  // view related
  void view(int, int, int, int, cmatrix&) const;
  cmatrix view(int, int, int, int) const;
  cmatrix trans() const;
  cmatrix row(int) const;
  cmatrix col(int) const;
  cmatrix rowT(int) const;
  cmatrix colT(int) const;

  cmatrix& incr(const cmatrix& in, double a=1);
  cmatrix& operator+=(const cmatrix& in){ incr(in, 1); return *this;}
  cmatrix& operator-=(const cmatrix& in){ incr(in,-1); return *this;}
  cmatrix& dotTimes(const cmatrix&);
  // the following functions should not care owner type
  cmatrix& makeIdentity();
  cmatrix& makeRand(dFunc0i);
  cmatrix& makeCustom(dFunc2i);
  cmatrix& op(dFunc1d);
  cmatrix& op(dFunc1d2i);
  cmatrix& operator*=(double);
  cmatrix& operator=(double);

  //
  cmatrix& inverse(const cmatrix& in);
  cmatrix& inverse(const cmatrix& in, double& det);
  cmatrix& inverse(){ return inverse(*this); }
  cmatrix& combine(double alpha, const cmatrix& a, double beta, const cmatrix& b);
  cmatrix& add(const cmatrix& a, const cmatrix& b);
  cmatrix& sub(const cmatrix& a, const cmatrix& b);
  // alpha * A * B + beta * C
  cmatrix& gemm(const cmatrix& a, const cmatrix& b, double beta=0, double alpha=1);
  // cmatrix& gemm(double beta, double alpha, const cmatrix& a, const cmatrix& b){ return gemm(a,b,beta,alpha); }
  
  cmatrix& block(const IdxArr& ria, const IdxArr& cia,       cmatrix& blk          );
  void setblock( const IdxArr& ria, const IdxArr& cia, const cmatrix& blk          ) const;
  void incrblock(const IdxArr& ria, const IdxArr& cia, const cmatrix& blk, double a=1) const;
  //  cmatrix& block(int r1, const IdxArr& cia, const cmatrix& in);
  //cmatrix& incrblock( double a,const IdxArr& ria, const IdxArr& cia, const cmatrix& in);
  void symmetrize();
  cmatrix& jacobi_nr(double *eval, cmatrix& evec);
  double det_nr();
  double logdet_nr();

  void symmEig_nr(cmatrix& eigval, cmatrix& eigvec);
  // A = U W VT, *this is AT
  void tsvd_nr(cmatrix &w, cmatrix &vt);
  void svbksb(const cmatrix &w, const cmatrix &vt, const cmatrix &b, cmatrix &x);

  double& value(int r, int c) const {assert(r>=0 && r<nr && c>=0 && c<nc);  return *(d + r*ld + c*inc); }
  double& dvalue(int i) const { return *(d + (inc+ld)*i); }
  double& rvalue(int c) const { assert(nr==1 && nc>0 && c<nc && c>=0);return *(d + c*inc); }
  double& cvalue(int r) const { return *(d + r*ld); }

  void open(FILE *fp);
  void store(FILE *fp);

};

#define bufCheckMemSize0 cautiousResize
#define memResize1 cautiousResize

#endif
