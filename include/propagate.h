#ifndef PROPAGATE_H
#define PROPAGATE_H

#include <iostream>
#include <iomanip>

#include <btas/common/blas_cxx_interface.h>
#include <btas/common/TVector.h>
#include <btas/DENSE/TArray.h>

using namespace btas;

namespace propagate {

   void step(PEPS<double> &,int);

   void update_vertical(int,int,PEPS<double> &,const DArray<5> &,const DArray<5> &,int);

   void construct_N_eff_vertical(int,int,PEPS<double> &,const DArray<5> &,const DArray<5> &,DArray<8> &,bool);

   void construct_rhs_vertical(int,int,PEPS<double> &,const DArray<5> &,const DArray<5> &,DArray<5> &,bool);

   void solve(DArray<4> &,DArray<3> &);

}

#endif
