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

   void update_horizontal(int,int,PEPS<double> &,const DArray<5> &,const DArray<5> &,int);

   void construct_lin_sys_vertical(int,int,PEPS<double> &,const DArray<6> &,const DArray<6> &,DArray<8> &,DArray<5> &,const DArray<7> &,const DArray<7> &,bool);

   void construct_lin_sys_horizontal(int,int,PEPS<double> &,const DArray<6> &,const DArray<6> &,DArray<8> &,DArray<5> &,const DArray<7> &,const DArray<7> &,bool);

   double cost_function_vertical(int,int,PEPS<double> &,const DArray<6> &,const DArray<6> &,const DArray<7> &,const DArray<7> &);

   double cost_function_horizontal(int,int,PEPS<double> &,const DArray<6> &,const DArray<6> &,const DArray<7> &,const DArray<7> &);

   void solve(DArray<8> &,DArray<5> &);

}

#endif
