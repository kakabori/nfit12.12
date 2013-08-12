#include "mydll.h"
#include <assert.h>
#include <cstring>

#define CMAT_MEMDEBUG

#ifdef CMAT_MEMDEBUG
int CmatrixNumAlloc=0;
#endif


void cmatrix::memdebug() const{
#ifdef CMAT_MEMDEBUG
  cout<<"numAlloc: "<<CmatrixNumAlloc<<endl;
#endif
}

cmatrix::cmatrix(int __nr, int __nc):
  nr(__nr), 
  nc(__nc)
{
  memsize = nr * nc;
  if( memsize > 0 ){
    d = (double *) malloc ( sizeof(double) * memsize );
#ifdef CMAT_MEMDEBUG
    CmatrixNumAlloc++;
#endif
  }
  ld = nc;
  inc = 1;
}

cmatrix::cmatrix(const cmatrix& in, enum CMAT_METHOD method){
  if( method == CMAT_COPY ) {
    nr = in.nr;
    ld = nc = in.nc;
    inc = 1;
    assert( ( memsize = nr * nc ) > 0 );
    d = (double *) malloc ( sizeof(double) * memsize );
#ifdef CMAT_MEMDEBUG
    CmatrixNumAlloc++;
#endif
    yl_blas_dmcopy( in.d, in.ld, in.inc, d, ld, inc, nr, nc);
  } else {
    memcpy(this, &in, sizeof(cmatrix) );
    memsize = -1;
  }
}
  

// should we take care of referenc count here
cmatrix::~cmatrix(){
  //cout<<"destructor "<<this<<' '<<d<<endl;
  if( isOwner() ){
    free(d);
#ifdef CMAT_MEMDEBUG
    CmatrixNumAlloc--;
#endif
  }
}


///////////////////////////////////////////////////////////
// size related
///////////////////////////////////////////////////////////

// need to debate about mem reallocation on
// the issue of reference count


// to be modified for nr*nc consistancy check
void cmatrix::memDirectRealloc(int __memsize){
  assert( memsize > 0 );
  memsize = __memsize;
  d = (double*) realloc(d , memsize);
}

//post: mem alloced if needed
void cmatrix::destructiveLazyResize(int __nr, int __nc){
  bool wasOwner = isOwner();
  int newsize = (nr = __nr) * (ld = nc = __nc);
  if( newsize <= 0 ){
    if( wasOwner ){
      free(d);
#ifdef CMAT_MEMDEBUG
      CmatrixNumAlloc--;
#endif
    }
    memsize = -1;
  } else {
    if( newsize > memsize ){
      if( wasOwner ){
	free(d);
#ifdef CMAT_MEMDEBUG
	CmatrixNumAlloc--;
#endif
      }
      memsize = newsize;
      d = (double*) malloc( sizeof(double) * memsize ); 
#ifdef CMAT_MEMDEBUG
      CmatrixNumAlloc++;
#endif
    }
  }
  inc = 1;
}

void cmatrix::destructiveResize(int __nr, int __nc){
  nr = __nr;
  nc = __nc;
  if( nr * nc <= 0 ){
    if( isOwner() ){
      free(d);
#ifdef CMAT_MEMDEBUG
      CmatrixNumAlloc--;
#endif
      memsize = -1;
    }
  } else {
    if( isOwner() ){
      d = (double*) realloc( d, sizeof(double) * (memsize = nr * nc) );
    } else {
      d = (double*) malloc( sizeof(double) * (memsize = nr * nc) ); 
#ifdef CMAT_MEMDEBUG
      CmatrixNumAlloc++;
#endif
    }
  }
  ld = nc;
  inc = 1;
}

void cmatrix::cautiousResize(int __nr, int __nc){
  if( __nr != nr || __nc != nc ){
    assert(isOwner() || nr*nc <=0 );
    destructiveLazyResize(__nr, __nc);
  }
}


cmatrix& cmatrix::resize1(int __nr, int __nc){
  assert( isOwner() );
  if(__nr <= nr ){
    if( __nc <= nc){
      nr = __nr;
      nc = __nc;
    } else {
      //      double *td = (double*)malloc(__nr*__nc*sizeof(double));
    }
  }
  return *this;
}



/////////////////////////////////////////////////////////////
// view related
/////////////////////////////////////////////////////////////

void cmatrix::view(int r1, int r2, int c1, int c2, cmatrix &viewer) const{
  assert( ! viewer.isOwner() );
  assert(r2 <= nr);
  assert(c2 <= nc);
  viewer.d = d + c1 * inc + r1 * ld;
  viewer.inc = inc;
  viewer.ld = ld;
  viewer.nr = r2 - r1;
  viewer.nc = c2 - c1;
}

cmatrix cmatrix::view(int r1, int r2, int c1, int c2) const{
  cmatrix viewer;
  viewer.d = d + c1 * inc + r1 * ld;
  viewer.inc = inc;
  viewer.ld = ld;
  viewer.nr = r2 - r1;
  viewer.nc = c2 - c1;
  return viewer;
}

cmatrix cmatrix::trans() const{
  cmatrix viewer;
  viewer.d = d;
  viewer.ld = inc;
  viewer.inc = ld;
  viewer.nr = nc;
  viewer.nc = nr;
  return viewer;
}

cmatrix cmatrix::row(int i) const{
  cmatrix viewer(*this, CMAT_VIEW);
  viewer.d = d+i*ld;
  viewer.nr = 1;
  return viewer;
}

cmatrix cmatrix::col(int i) const{
  cmatrix viewer(*this, CMAT_VIEW);
  viewer.d = d+i*inc;
  viewer.nc = 1;
  return viewer;
}

cmatrix cmatrix::rowT(int i) const{
  cmatrix viewer;
  viewer.d = d+i*ld;
  viewer.nc = 1;
  viewer.nr = nc;
  viewer.ld = inc;
  viewer.inc = ld;
  return viewer;
}

cmatrix cmatrix::colT(int i) const{
  cmatrix viewer;
  viewer.d = d+i*inc;
  viewer.nr = 1;
  viewer.nc = nr;
  viewer.ld = inc;
  viewer.inc = ld;
  return viewer;
}


/////////////////////////////////////////////////////////////////
// property related
/////////////////////////////////////////////////////////////////
double cmatrix::norm22(){
  double sum=0;
  for(int i=0; i<nr; i++){
    sum += yl_blas_dnrm22(nc, d+i*ld, inc);
  }
  return sum;
}

double cmatrix::sum(){
  double sum=0;
  for(int i=0; i<nr; i++){
    sum += yl_blas_dsum(nc, d+i*ld, inc);
  }
  return sum;
}

double cmatrix::asum(){
  double sum=0;
  for(int i=0; i<nr; i++){
    sum += yl_blas_dasum(nc, d+i*ld, inc);
  }
  return sum;
}
double cmatrix::amax(){
  double max=value(0,0);
  for(int i=0; i<nr; i++)
    yl_blas_idamax2(nc, d+i*ld, inc, &max);
  return max;
}

cmatrix& cmatrix::inverse(const cmatrix& inv, double& det){
  assert(inv.nr==inv.nc);
  cmatrix acpy( inv, CMAT_COPY );
  cautiousResize(inv.nr, inv.nr);
  yl_cblas_nr_minverse(acpy.d, acpy.ld, acpy.inc,
		       d, ld, inc,
		       nr);
  det = acpy.value(0,0);
  for(int i=1; i<acpy.nr; i++)
    det *= acpy.dvalue(i);
  return *this;
}


cmatrix& cmatrix::inverse(const cmatrix& inv){
  assert(inv.nr==inv.nc);
  cmatrix acpy( inv,  CMAT_COPY );
  cautiousResize(inv.nr, inv.nr);
  //cout<<acpy<<endl;
  yl_cblas_nr_minverse(acpy.d, acpy.ld, acpy.inc,
		       d, ld, inc,
		       nr);
  return *this;
}

cmatrix& cmatrix::jacobi_nr(double *eval, cmatrix& evec){
  assert(nr == nc);
  int nrot;
  cmatrix acpy( *this, CMAT_COPY);
  evec.cautiousResize(nr, nr);
  yl_cblas_nr_jacobi(nr,
		     acpy.d, acpy.ld, acpy.inc,
		     evec.d, evec.ld, evec.inc,
		     eval, &nrot);
  return evec;
}


double cmatrix::det_nr(){
  assert(nr==nc);
  cmatrix acpy( *this,  CMAT_COPY );
  return yl_cblas_nr_mdeterminant(nr, acpy.d, acpy.ld, acpy.inc);
}


double cmatrix::logdet_nr(){
  assert(nr==nc);
  cmatrix acpy( *this, CMAT_COPY );
  return yl_cblas_nr_mlogdet(nr, acpy.d, acpy.ld, acpy.inc);
}

void cmatrix::symmEig_nr(cmatrix& eigval, cmatrix& eigvec){
  assert(nr == nc);
  eigval.cautiousResize(1 , nr);
  eigvec.cautiousResize(nr, nr);
  eigvec = *this;
  double dd[nr];	// auxiliary vectors
  double ee[nr];	// auxiliary vectors
  
  yl_blas_nr_tred2(nr, eigvec.d, eigvec.ld, eigvec.inc, dd, ee);
  yl_blas_nr_tqli( nr, eigvec.d, eigvec.ld, eigvec.inc, dd, ee);
  yl_blas_dcopy(nr, dd, 1, eigval.d, eigval.inc);
}


/////////////////////////////////////////////
// A = U W VT, input: *this is AT, 
// output: *this is UT, vt is VT
/////////////////////////////////////////////
void cmatrix::tsvd_nr(cmatrix &w, cmatrix &vt){
  w.cautiousResize( 1, nr );
  vt.cautiousResize( nr, nr );
  yl_cblas_nr_svdcmp(  nc, nr,
		        d, ld,    inc,
		      w.d, w.inc,
		     vt.d, vt.ld,  vt.inc);
}

void cmatrix::svbksb(const cmatrix &w, const cmatrix &vt, const cmatrix &b, cmatrix &x){
  assert(w.nc == nr && vt.nr == nr && b.nc == nc);
  x.cautiousResize(1, nr);
  yl_cblas_nr_svbksb( nc, nr,
		      d, ld, inc,
		      w.d, w.inc,
		      vt.d, vt.ld, vt.inc,
		      b.d, b.inc,
		      x.d, x.inc);
}

/////////////////////////////////////////////////////////////////
// operation related
/////////////////////////////////////////////////////////////////


cmatrix& cmatrix::operator=(double a){
  yl_blas_dmset(a , d, ld, inc, nr, nc*inc);
  return *this;
}

cmatrix& cmatrix::makeIdentity(){
  *this=(0.0);
  yl_blas_dset( yfmin(nr, nc), 1 , d, ld+inc);
  return *this;
}

cmatrix& cmatrix::makeRand(dFunc0i f){
  for(int i=0; i<nr; i++)
    for(int j=0; j<nc; j++)
      value(i,j) = f();
  return *this;
}

cmatrix& cmatrix::makeCustom(dFunc2i f){
  for(int i=0; i<nr; i++)
    for(int j=0; j<nc; j++)
      value(i,j) = f(i,j);
  return *this;
}

cmatrix& cmatrix::op(dFunc1d2i f){
  for(int i=0; i<nr; i++)
    for(int j=0; j<nc; j++)
      value(i,j) = f(value(i,j),i,j);
  return *this;
}

cmatrix& cmatrix::op(dFunc1d f){
  for(int i=0; i<nr; i++)
    for(int j=0; j<nc; j++)
      value(i,j) = f(value(i,j));
  return *this;
}


cmatrix& cmatrix::operator=(const cmatrix& in){
  cautiousResize(in.nr, in.nc);
  yl_blas_dmcopy( in.d, in.ld, in.inc,
		  d, ld, inc,
		  nr, nc);
  return *this;
}

cmatrix& cmatrix::operator*=(double a){
  yl_blas_dmscal(a, d, ld, inc, nr, nc*inc);
  return *this;
}

cmatrix& cmatrix::incr(const cmatrix& in, double a){
  assert(nr==in.nr && nc==in.nc);
  yl_blas_dmaxpy( in.d, in.ld, in.inc,
		  d, ld, inc,
		  nr, nc, a);
  return *this;

}

cmatrix& cmatrix::dotTimes(const cmatrix& in){
  assert(nr==in.nr && nc==in.nc);
  for(int i=0; i<nr; i++)
    for(int j=0; j<nc; j++)
      value(i,j) *= in.value(i,j);
  return *this;
}


void cmatrix::symmetrize(){
  assert(nr == nc);
  for(int i=0; i<nr; i++){
    for(int j=i+1; j<nr; j++){
      //assert( fabs(value(i,j)- value(j,i))<1e-7);
      value(i,j) = value(j,i) = 0.5*(value(i,j)+value(j,i));
    }
  }
}


/////////////////////////////////////////////////////////////////
// buffer operation related
/////////////////////////////////////////////////////////////////


cmatrix& cmatrix::combine(double alpha,  const cmatrix& a,  double beta, const cmatrix& b){
  assert(a.nr==b.nr && a.nc==b.nc);
  cautiousResize(a.nr, a.nc);
  yl_blas_dmaxpbytc(a.d, a.ld, a.inc,
		  b.d, b.ld, b.inc,
		  d, ld, inc,
		  nr, nc, alpha, beta);
  return *this;
}



cmatrix& cmatrix::add( const cmatrix& a, const cmatrix& b){
  assert(a.nr==b.nr && a.nc==b.nc);
  cautiousResize(a.nr, a.nc);
  yl_blas_dmxpytc(a.d, a.ld, a.inc,
		  b.d, b.ld, b.inc,
		  d, ld, inc,
		  nr, nc);
  return *this;
}


cmatrix& cmatrix::sub(const cmatrix& a, const cmatrix& b){
  assert(a.nr==b.nr && a.nc==b.nc);
  cautiousResize(a.nr, a.nc);
  yl_blas_dmxmytc(a.d, a.ld, a.inc,
		  b.d, b.ld, b.inc,
		  d, ld, inc,
		  nr, nc);
  return *this;
}


// alpha * A * B + beta * C
cmatrix& cmatrix::gemm(const cmatrix& a, const cmatrix& b, double beta, double alpha){
  assert(a.nc == b.nr);
  cautiousResize(a.nr, b.nc);
  yl_cblas_dgemm(a.d, a.ld, a.inc,
		 b.d, b.ld, b.inc,
		 d, ld, inc,
		 nr, a.nc, nc, alpha, beta);
  return *this;
}


/////////////////////////////////////////////////////////////////
// block operation related
/////////////////////////////////////////////////////////////////


cmatrix& cmatrix::block(const IdxArr& ria, const IdxArr& cia, cmatrix& blk){
  int _nr = ria.size();
  int _nc = cia.size();

  if( _nr > 0 && _nc > 0 ){
    blk.cautiousResize(_nr, _nc);
    for(int i = 0; i < _nr; i++){
      double *td = d + ria[i] * ld;
      double *blkd = blk.d + i * blk.ld;
      for(int j = 0; j < _nc; j++)
	blkd[ j * blk.inc ] = td[ cia[j] * inc ];
    }
  } else if( _nr <= 0 ){
    _nr = nr;
    blk.cautiousResize(_nr, _nc);
    for(int i = 0; i < _nr; i++){
      double *td = d + i * ld;
      double *blkd = blk.d + i * blk.ld;
      for(int j = 0; j < _nc; j++)
	blkd[ j * blk.inc ] = td[ cia[j] * inc ];
    }
  } else {
    _nc = nc;
    blk.cautiousResize(_nr, _nc);
    for(int i = 0; i < _nr; i++)
      yl_blas_dcopy(_nc, d + ria[i] * ld, inc, blk.d + i * blk.ld, blk.inc);
  }
  return blk;
}

void cmatrix::setblock(const IdxArr& ria, const IdxArr& cia, const cmatrix& blk) const{
  int _nr = ria.size();
  int _nc = cia.size();

  if( _nr > 0 && _nc > 0 ){
    for(int i = 0; i < _nr; i++){
      double *td = d + ria[i] * ld;
      double *blkd = blk.d + i * blk.ld;
      for(int j = 0; j < _nc; j++)
	td[ cia[j] * inc ] = blkd[ j * blk.inc ];
    }
  } else if( _nr <= 0 ){
    _nr = nr;
    for(int i = 0; i < _nr; i++){
      double *td = d + i * ld;
      double *blkd = blk.d + i * blk.ld;
      for(int j = 0; j < _nc; j++)
	td[ cia[j] * inc ] = blkd[ j * blk.inc ];
    }
  } else {
    _nc = nc;
    for(int i = 0; i < _nr; i++)
      yl_blas_dcopy(_nc,  blk.d + i * blk.ld, blk.inc, d + ria[i] * ld, inc);
  }
  assert(_nr == blk.nr && _nc == blk.nc);
}

void cmatrix::incrblock(const IdxArr& ria, const IdxArr& cia, const cmatrix& blk, double a) const{
  int _nr = ria.size();
  int _nc = cia.size();

  if( _nr > 0 && _nc > 0 ){
    for(int i = 0; i < _nr; i++){
      double *td = d + ria[i] * ld;
      double *blkd = blk.d + i * blk.ld;
      for(int j = 0; j < _nc; j++)
	td[ cia[j] * inc ] += a * blkd[ j * blk.inc ];
    }
  } else if( _nr <= 0 ){
    _nr = nr;
    for(int i = 0; i < _nr; i++){
      double *td = d + i * ld;
      double *blkd = blk.d + i * blk.ld;
      for(int j = 0; j < _nc; j++)
	td[ cia[j] * inc ] += a * blkd[ j * blk.inc ];
    }
  } else {
    _nc = nc;
    for(int i = 0; i < _nr; i++)
      yl_blas_daxpy(_nc, a,
		    blk.d + i * blk.ld,
		    blk.inc,
		    d + ria[i] * ld,
		    inc);
  }
  assert(_nr == blk.nr && _nc == blk.nc);
}

/////////////////////////////////////////////////////////////////
// io related
/////////////////////////////////////////////////////////////////



ostream& operator<<( ostream& out, const cmatrix& cm){
  yl_blas_mprint(out, cm.nr, cm.nc, cm.d, cm.ld, cm.inc, 1);
  out<<(cm.isOwner() ? " owner " : " user ");
  out<<" nr: "<<cm.nr;
  out<<" nc: "<<cm.nc;
  out<<" ld: "<<cm.ld;
  out<<" inc: "<<cm.inc;
  out<<" d: "<<cm.d<<endl;
  return out;
}

void cmatrix::open(FILE *fp){
  int _nr, _nc;
  fread( &_nr, sizeof(int), 1, fp);
  fread( &_nc, sizeof(int), 1, fp);
  cautiousResize(_nr,_nc);
  for(int i=0; i<nr; i++){
    for(int j=0; j<nc; j++){
      fread( &value(i,j) , sizeof(double), 1, fp);
    }
  }
} 

void cmatrix::store(FILE *fp){
  fwrite( &nr, sizeof(int), 1, fp);
  fwrite( &nc, sizeof(int), 1, fp);
  for(int i=0; i<nr; i++){
    for(int j=0; j<nc; j++){
      fwrite( &value(i,j) , sizeof(double), 1, fp);
    }
  }
}


//--------------------------------
// constructor,  destructor(inline)

IdxArr::IdxArr(int msz): ia(NULL), ia_size(0), memsize(msz){
  if( memsize > 0 ){
    ia = (int *) malloc( memsize * sizeof(int) );
    assert( ia );
  }
}

IdxArr::IdxArr(const IdxArr &in): ia(NULL), ia_size(0), memsize(0){
  _copy( in );
}

IdxArr::~IdxArr() {
  // cout<<ia<<endl;
  free(ia);
  // cout<<"IdxArrd"<<endl;
}

//----------------------------------
//private helpers

void IdxArr::_increaseMemSize(){
  memsize = memsize>0 ? (2*memsize) : 8 ;
  ia = (int*) realloc( ia, memsize * sizeof(int) );
  assert( ia );
}

void IdxArr::_copy(const IdxArr &in) {
  ia_size = in.ia_size;
  if( memsize < ia_size ){
    memsize = in.ia_size;
    free( ia );
    ia = (int*) malloc( memsize * sizeof(int) );
    assert( ia );
  }
  memcpy(ia, in.ia, ia_size * sizeof(int));
}

//--------------------------------------
//element manipulations

//addIdx checks to see that the index array doesn't have the index x already
void IdxArr::addIdx(int x) {
  for(int i = 0; i < ia_size; i++)
    if(ia[i] == x) return;
  if(ia_size >= memsize) _increaseMemSize();
  ia[ia_size] = x;
  ia_size++;
}


// assume no overlap
IdxArr& IdxArr::append(const IdxArr& in) {
  for(int i=0; i<in.ia_size; i++)    addIdx(in.ia[i]);
  return *this;
}


//---------------------------------
bool IdxArr::exist(int idx) const{
  for(int i=0; i<ia_size; i++)
    if(ia[i] == idx) return true;
  return false;
}

void IdxArr::operator+=(int c){
  for(int i=0; i<ia_size; i++)
    ia[i] += c;
}

ostream& operator<<(ostream& out, const IdxArr& idxarr){
  out<<"IdxArr "<<idxarr.ia<<'{';
  for(int i=0; i<idxarr.ia_size; i++)
    out<<' '<<idxarr.ia[i];
  out<<'}';
  return out;
}




