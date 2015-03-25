#ifndef PROPAGATE_H
#define PROPAGATE_H

#include <iostream>
#include <iomanip>

#include <btas/common/blas_cxx_interface.h>
#include <btas/common/TVector.h>
#include <btas/DENSE/TArray.h>

using namespace btas;

namespace propagate {

   //!define a variable for the direction of the propagation
   enum PROP_DIR {

      //!vertical gate on two nearest neigbour peps (row,col) --> (row+1,col)
      VERTICAL=0,

      //!horizontal gate on two nearest neigbour peps (row,col) --> (row,col+1)
      HORIZONTAL=1,
      
      //!diagonal gate connect two next-nearest neigbour peps, from top left to bottom right (row+1,col) --> (row,col+1)
      DIAGONAL_LURD=2,

      //!diagonal gate connect two next-nearest neigbour peps, from bottom left to top right (row,col) --> (row+1,col+1)
      DIAGONAL_LDRU=3
   
   };

   void step(PEPS<double> &,int);

   void solve(DArray<8> &,DArray<5> &);

   //updates
   template<size_t M>
      void update(const PROP_DIR &,int,int,PEPS<double> &,const DArray<M> &,const DArray<M> &,int);

   //quasi-canonicalization of the environment
   template<size_t M>
      void canonicalize(const PROP_DIR &,int,int,PEPS<double> &,const DArray<M> &,const DArray<M> &,const DArray<M+2> &,const DArray<M+2> &);

   //sweeping sections
   template<size_t M>
      void sweep(const PROP_DIR &,int,int,PEPS<double> &,const DArray<6> &,const DArray<6> &,
            
            const DArray<M> &,const DArray<M> &,const DArray<M+2> &,const DArray<M+2> &,const DArray<M+2> &,const DArray<M+2> &,int);

   //construct intermediate objects for N_eff construction
   template<size_t M>
      void construct_intermediate(const PROP_DIR &,int,int,const PEPS<double> &,const DArray<5> &,
            
            const DArray<M> &,const DArray<M> &,DArray<M+2> &,DArray<M+2> &,DArray<M+2> &,DArray<M+2> &);

   //linear systems construct
   template<size_t M>
      void construct_lin_sys(const PROP_DIR &,int,int,PEPS<double> &,const DArray<6> &,const DArray<6> &, DArray<8> &,DArray<5> &,

            const DArray<M> &,const DArray<M> &,const DArray<M+2> &,const DArray<M+2> &,const DArray<M+2> &,const DArray<M+2> &,bool);
   
   //calculate the effective environment
   template<size_t M>
      void calc_N_eff(const PROP_DIR &,int,int,PEPS<double> &, DArray<8> &, const DArray<M> &,const DArray<M> &, 

            const DArray<M+2> &,const DArray<M+2> &,bool);
 
   //calculate right hand of the linear system of equations
   template<size_t M>
      void calc_rhs(const PROP_DIR &,int,int,PEPS<double> &,const DArray<6> &,const DArray<6> &,DArray<5> &,

            const DArray<M> &,const DArray<M> &,const DArray<M+2> &,const DArray<M+2> &,const DArray<M+2> &,const DArray<M+2> &,bool);

   //initialization by SVD
   void initialize(const PROP_DIR &,int,int,const DArray<6> &,const DArray<6> &,PEPS<double> &);

   //restore after update
   void equilibrate(const PROP_DIR &,int,int,PEPS<double> &);

   void diagonalize(DArray<8> &,DArray<1> &);

   void get_X(const DArray<8> &,const DArray<1> &,DArray<5> &);

}

#endif
