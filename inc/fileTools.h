#ifndef GUARD_FILETOOLS_H
#define GUARD_FILETOOLS_H

#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <assert.h>
#include <stdexcept>


using std::vector; using std::string; using std::ifstream;
using std::getline; using std::istringstream; using std::ofstream;
using std::cout; using std::endl; using std::cerr;

/******************************************************************************
Read an ASCII formatted matrix into the input 1D vector, vec. 
Returns the total number of columns. The total number of rows
can be calculated via vec.size() divided by the function's 
returned value.

vec[i * width + j] returns the matrix element in ith row and 
jth column
******************************************************************************/
template <class T>
int readMatrixFromFile(const char *filename, vector<T>& vec)
{
  T tmp;
  string line;
  int count = 0, width = 0;
  
  ifstream infile;
  infile.open(filename);
  if (!infile) {
    cerr << filename << ": not found" << endl;
    //throw domain_error("ERROR: could not open the file");
  }
  
  // Go through the first line and determine the width
  getline(infile, line);
  istringstream iss(line);
  while (iss >> tmp) {
    vec.push_back(tmp);
    width++;
  }
  
  // Go through the rest of the file
  while (getline(infile, line)) {
    istringstream iss2(line);
    count = 0;
    while (iss2 >> tmp) {
      vec.push_back(tmp);
      count++;
    }
    // Complain if the width is not constant
    if (count != width) {
      cout << count << " " << width << endl;
      cerr << "what the heck!?" << endl;
      //throw domain_error("the input file was not formatted correctly");
    }
  }
  
  return width;
}



#endif
