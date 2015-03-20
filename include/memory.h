#ifndef MEMORY_H
#define MEMORY_H

#include <iostream>
#include <fstream>
#include <vector>

#include <btas/common/blas_cxx_interface.h>
#include <btas/common/TVector.h>
#include <btas/DENSE/TArray.h>

using std::ostream;
using std::vector;

/**
 * @author Brecht Verstichel
 * @data 19-03-2015\n\n
 * namesapce used for memory allocation, because it's slow!
 */
namespace memory {

      void init();

      //!public variable for environment construction
      extern std::vector< DArray<4> > R4;

};

#endif
