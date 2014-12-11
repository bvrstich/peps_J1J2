#ifndef PROPAGATE_H
#define PROPAGATE_H

#include <iostream>
#include <iomanip>

#include <btas/common/blas_cxx_interface.h>
#include <btas/common/TVector.h>
#include <btas/DENSE/TArray.h>

using namespace btas;

namespace propagate {

   void step(bool,PEPS<double> &,int);

   void construct_reduced_tensor(char,char,const DArray<5> &,DArray<4> &,DArray<3> &);

   void invert(DArray<2> &);

   void prop_local(PEPS<double> &);

   void solve(DArray<4> &N_eff,DArray<3> &b);

   void calc_N_eff(char,int,const DArray<3> &,const DArray<4> &,const DArray<3> &, const DArray<4> &,DArray<4> &);

   void calc_N_eff(char,int,int,const DArray<4> &,const DArray<4> &,const DArray<4> &,const DArray<4> &,DArray<4> &);

   void canonicalize(bool,DArray<4> &,DArray<3> &a_L,DArray<4> &QL,DArray<3> &a_R,DArray<4> &QR);

   void update(bool,const DArray<4> &N_eff,DArray<3> &,DArray<3> &,int);

   void check_N_eff(const DArray<4> &N_eff,const DArray<3> &,const DArray<3> &);

   double cost_function(const DArray<4> &N_eff,const DArray<4> &b,const DArray<3> &a_L,const DArray<3> &a_R);

}

#endif
