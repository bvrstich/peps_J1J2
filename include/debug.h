#ifndef DEBUG_H
#define DEBUG_H

#include <iostream>
#include <iomanip>

#include <btas/common/blas_cxx_interface.h>
#include <btas/common/TVector.h>
#include <btas/DENSE/TArray.h>

using namespace btas;
using namespace propagate;

namespace debug {

   //cost function for testing the ALS solution of the tensor update
   template<size_t M>
      double cost_function(const PROP_DIR &,int,int,PEPS<double> &,const DArray<6> &,const DArray<6> &,

            const DArray<M> &,const DArray<M> &,const DArray<M+2> &,const DArray<M+2> &,const DArray<M+2> &,const DArray<M+2> &);

   void diagonalize(DArray<8> &,DArray<1> &);

}

#endif
