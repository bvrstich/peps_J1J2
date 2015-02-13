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

   void solve(DArray<8> &,DArray<5> &);

   //updates
   template<size_t M>
      void update_vertical(int,int,PEPS<double> &,const DArray<M> &,const DArray<M> &,int);

   template<size_t M>
      void update_horizontal(int,int,PEPS<double> &,const DArray<M> &,const DArray<M> &,int);

   void update_diagonal_lurd(int,int,PEPS<double> &,const DArray<5> &,const DArray<5> &,int);

   //sweeping sections
   template<size_t M>
      void sweep_vertical(int,int,PEPS<double> &,const DArray<6> &,const DArray<6> &,const DArray<M> &,

            const DArray<M> &,const DArray<M+2> &,const DArray<M+2> &, int);

   template<size_t M>
      void sweep_horizontal(int,int,PEPS<double> &,const DArray<6> &,const DArray<6> &,const DArray<M> &,

            const DArray<M> &,const DArray<M+2> &,const DArray<M+2> &, int);

   //construct intermediate objects for N_eff construction
   template<size_t M>
      void construct_intermediate_vertical(int,int,const PEPS<double> &,const DArray<M> &,const DArray<M> &,DArray<M+2> &,DArray<M+2> &);

   template<size_t M>
      void construct_intermediate_horizontal(int,int,const PEPS<double> &,const DArray<M> &,const DArray<M> &,DArray<M+2> &,DArray<M+2> &);

   //linear systems construct
   template<size_t M>
      void construct_lin_sys_vertical(int,int,PEPS<double> &,const DArray<6> &,const DArray<6> &, DArray<8> &,DArray<5> &,

            const DArray<M> &,const DArray<M> &,const DArray<M+2> &,const DArray<M+2> &,bool);

   template<size_t M>
      void construct_lin_sys_horizontal(int,int,PEPS<double> &,const DArray<6> &,const DArray<6> &, DArray<8> &,DArray<5> &,

            const DArray<M> &,const DArray<M> &,const DArray<M+2> &,const DArray<M+2> &,bool);

   void construct_lin_sys_diagonal_lurd(int,int,PEPS<double> &,const DArray<6> &,const DArray<6> &,DArray<8> &,DArray<5> &,

         const DArray<5> &,const DArray<5> &,const DArray<7> &,const DArray<7> &,bool);

   //cost functions
   template<size_t M>
      double cost_function_vertical(int,int,PEPS<double> &,const DArray<6> &,const DArray<6> &,

            const DArray<M> &,const DArray<M> &,const DArray<M+2> &,const DArray<M+2> &);

   template<size_t M>
      double cost_function_horizontal(int,int,PEPS<double> &,const DArray<6> &,const DArray<6> &,
            
            const DArray<M> &, const DArray<M> &,const DArray<M+2> &,const DArray<M+2> &);

   double cost_function_diagonal_lurd(int,int,PEPS<double> &,const DArray<6> &,const DArray<6> &,const DArray<5> &,

         const DArray<5> &,const DArray<7> &,const DArray<7> &);

   //initialization by SVD
   void initialize_vertical(const DArray<6> &,const DArray<6> &,DArray<5> &,DArray<5> &);

   void initialize_horizontal(const DArray<6> &,const DArray<6> &,DArray<5> &,DArray<5> &);

   //restore after update
   void equilibrate_vertical(DArray<5> &,DArray<5> &);

   void equilibrate_horizontal(DArray<5> &,DArray<5> &);

}

#endif
