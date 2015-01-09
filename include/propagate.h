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

   void construct_reduced_tensor(char,char,const DArray<5> &,DArray<4> &,DArray<3> &);

   void calc_vertical_N_eff(char option,int col,const DArray<5> &L,const DArray<4> &QL,const DArray<5> &R,const DArray<4> &QR,DArray<4> &N_eff_n);

   void calc_horizontal_N_eff(char option,int col,const DArray<5> &L,const DArray<4> &QL,const DArray<5> &R,const DArray<4> &QR,DArray<4> &N_eff_n);

   void canonicalize_n(DArray<4> &,DArray<3> &a_L,DArray<4> &QL,DArray<3> &a_R,DArray<4> &QR);

   void invert(DArray<2> &);

   void update_n(const DArray<4> &,DArray<3> &,DArray<3> &,int);

   double cost_function_n(const DArray<4> &,const DArray<4> &,const DArray<3> &,const DArray<3> &);

   void solve(DArray<4> &,DArray<3> &);

}

#endif
