#ifndef __UNIFORM_HH
#define __UNIFORM_HH


// ***********************************************************************
// ***********************************************************************
// ***********************************************************************
// 
// uniform - class for uniform distribution
// 
// ***********************************************************************
// ***********************************************************************
// ***********************************************************************


class uniform{

public:

  // constructor(s)
  uniform(int __dimension);

  // destructor(s)
  ~uniform();

  // copy operator
  //uniform(const uniform& x); 

  // direct access to the elements of the uniform distribution
  double &lower(int i);
  double &upper(int i);

  // sample
  void generateSample(cmatrix &data, int numSamples = 1);


private:
  // memory allocation
  void allocate();
  void deallocate();
  // internal data structures
  int dimension;
  dduple *limits;
};



#endif
