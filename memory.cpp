#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <complex>
#include <omp.h>

using std::cout;
using std::endl;
using std::vector;
using std::complex;
using std::ofstream;

#include "include.h"

using namespace global;

namespace memory{

   std::vector< DArray<4> > R4;

   /**
    * initialize the memory
    */
   void init(){

      R4.resize(Lx-1);

   }

}
