#include <cstdlib>
#include <cstdio>
#include <assert.h>
#include <cmath>
#include <values.h>
#include <stack>
#include <signal.h>
#include <unistd.h>

#include <typeinfo>
#include <fstream>
#include <iostream>
#include <string>


#include "mydll.h"


using namespace std;

double normalSample(double mean, double variance){
  assert(variance > 0);
  return mean+variance*gasdev();
}
