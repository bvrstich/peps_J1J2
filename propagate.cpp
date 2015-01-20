#include <iostream>
#include <iomanip>
#include <fstream>
#include <complex>
#include <chrono>

using std::cout;
using std::endl;
using std::ostream;
using std::complex;
using std::ofstream;

#include "include.h"

using namespace btas;
using namespace global;

namespace propagate {

   /**
    * propagate the peps one imaginary time step
    * @param peps the PEPS to be propagated
    * @param n_sweeps the number of sweeps performed for the solution of the linear problem
    */
   void step(PEPS<double> &peps,int n_sweeps){

      enum {i,j,k,l,m,n,o};

      // -------------------------------------------//
      // --- !!! (1) the bottom two rows (1) !!! ---// 
      // -------------------------------------------//

      //containers for the renormalized operators
      vector< DArray<5> > R(Lx - 1);
      DArray<5> L;

      //construct the full top environment:
      env.calc('T',peps);

      //and the bottom row environment
      env.gb(0).fill('b',peps);

      //initialize the right operators for the bottom row
      contractions::init_ro('b',peps,R);

      // for(int col = 0;col < Lx - 1;++col){
      int col = 0;

      //--- (1) update the vertical pair on column 'col' ---
      update_vertical(0,col,peps,L,R[0],n_sweeps); 
/*
      // --- (2) update the horizontal pair on column 'col'-'col+1' ---

      contractions::update_L('b',col,peps,L);

      //}

      //update the bottom row for the new peps
      env.gb(0).fill('b',peps);

      // ---------------------------------------------------//
      // --- !!! (2) the middle rows (1 -> Ly-2) (2) !!! ---// 
      // ---------------------------------------------------//

      //renormalized operators for the middle sites
      vector< DArray<4> > RO(Lx - 1);
      DArray<4> LO;

      for(int row = 1;row < Ly-1;++row){

      //first create right renormalized operator
      contractions::init_ro(false,'H',row,peps,RO);

      for(int col = 0;col < Lx - 1;++col){

      //first construct the reduced tensors of the first pair to propagate
      construct_reduced_tensor('H','L',peps(row,col),QL,a_L);
      construct_reduced_tensor('H','R',peps(row,col+1),QR,a_R);

      //calculate the effective environment N_eff
      calc_N_eff('H',row,col,LO,QL,RO[col + 1],QR,N_eff);

      //make environment close to unitary before the update
      canonicalize(full,N_eff,a_L,QL,a_R,QR);

      //now do the update! Apply the gates!
      update(full,N_eff,a_L,a_R,n_sweeps);

      //and expand back to the full tensors
      Contract(1.0,QL,shape(i,j,k,o),a_L,shape(o,m,n),0.0,peps(row,col),shape(i,j,m,k,n));
      Contract(1.0,a_R,shape(i,j,k),QR,shape(k,o,m,n),0.0,peps(row,col+1),shape(i,o,j,m,n));

      //first construct a double layer object for the newly updated bottom 
      contractions::update_L('H',row,col,peps,LO);

      }

      //finally update the 'bottom' environment for the row
      env.add_layer('b',row,peps);

      }

      // ------------------------------------------//
      // --- !!! (3) the top row (Ly-1) (3) !!! ---// 
      // ------------------------------------------//

      //make the right operators
      contractions::init_ro(false,'t',peps,R);

      for(int col = 0;col < Lx - 1;++col){

      //first construct the reduced tensors of the first pair to propagate
      construct_reduced_tensor('H','L',peps(Ly-1,col),QL,a_L);
      construct_reduced_tensor('H','R',peps(Ly-1,col + 1),QR,a_R);

      //calculate the effective environment N_eff
      calc_N_eff('t',col,L,QL,R[col + 1],QR,N_eff);

      //make environment close to unitary before the update
      canonicalize(full,N_eff,a_L,QL,a_R,QR);

      //now do the update! Apply the gates!
      update(full,N_eff,a_L,a_R,n_sweeps);

      //and expand back to the full tensors
      Contract(1.0,QL,shape(i,j,k,o),a_L,shape(o,m,n),0.0,peps(Ly-1,col),shape(i,j,m,k,n));
      Contract(1.0,a_R,shape(i,j,k),QR,shape(k,o,m,n),0.0,peps(Ly-1,col+1),shape(i,o,j,m,n));

      contractions::update_L('t',col,peps,L);

   }

   //get the norm matrix
   contractions::update_L('t',Lx-1,peps,L);

   //scale the peps
   peps.scal(1.0/sqrt(L(0,0,0)));

   */
   }

 
   /** 
    * wrapper function solve linear system: N_eff * x = b
    * @param N_eff input matrix
    * @param rhs right hand side input and x output
    */
   void solve(DArray<8> &N_eff,DArray<5> &rhs){

      //int n = N_eff.shape(0) * N_eff.shape(1);

      //lapack::potrf(CblasRowMajor,'U',n, N_eff.data(), n);

      //lapack::potrs(CblasRowMajor,'U',n,d, N_eff.data(), n,b.data(),d);

   }

   /**
    * update the tensors in a sweeping fashion
    * @param row , the row index of the bottom site
    * @param col column index of the vertical column
    * @param peps, full PEPS object before update
    * @param L Left environment contraction
    * @param R Right environment contraction
    * @param n_iter nr of sweeps in the ALS algorithm
    */
   void update_vertical(int row,int col,PEPS<double> &peps,const DArray<5> &L,const DArray<5> &R,int n_iter){

      DArray<5> rhs;
      DArray<8> N_eff;

      DArray<7> tmp7;
      Gemm(CblasNoTrans,CblasNoTrans,1.0,env.gt(0)[0],R,0.0,tmp7);

      cout << tmp7 << endl;

      //start sweeping
      int iter = 0;

      while(iter < n_iter){

         // --(1)-- top site

         //construct effective environment, matrix for linear system for top site
         construct_N_eff_vertical(row,col,peps,L,R,N_eff,true);

         //construct right hand side of linear system for top site
         construct_rhs_vertical(row,col,peps,L,R,rhs,true);

         //solve the system

         // --(2)-- bottom site

         //construct effective environment for bottom site
         construct_N_eff_vertical(row,col,peps,L,R,N_eff,false);

         //construct right hand side of linear system for bottom site
         construct_rhs_vertical(row,col,peps,L,R,rhs,false);

         //solve the system

         //repeat until converged
         ++iter;

      }

   }

   /**
    * construct the single-site effective environment for the vertical gate, for top or bottom site
    * @param row , the row index of the bottom site
    * @param col column index of the vertical column
    * @param peps, full PEPS object before update
    * @param L Left environment contraction
    * @param R Right environment contraction
    * @param N_eff output object, contains N_eff on output
    * @param top boolean flag for top or bottom site of vertical gate
    */
   void construct_N_eff_vertical(int row,int col,PEPS<double> &peps,const DArray<5> &L,const DArray<5> &R,DArray<8> &N_eff,bool top){

      if(top){//top site of vertical, so site (row,col+1) environment

      }
      else{//bottom site (row,col)

      }

   }

   /**
    * construct the single-site right hand side of the linear system N_eff x = rhs, for top or bottom site
    * @param row , the row index of the bottom site
    * @param col column index of the vertical column
    * @param peps, full PEPS object before update
    * @param L Left environment contraction
    * @param R Right environment contraction
    * @param rhs output object, contains rhs on output
    * @param top boolean flag for top or bottom site of vertical gate
    */
   void construct_rhs_vertical(int row,int col,PEPS<double> &peps,const DArray<5> &L,const DArray<5> &R,DArray<5> &rhs,bool top){

      if(top){//top site of vertical, so site (row,col+1) environment

      }
      else{//bottom site (row,col)

      }

   }

}
