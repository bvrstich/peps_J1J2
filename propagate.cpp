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
    * wrapper function solve general linear system: N_eff * x = b
    * @param N_eff input matrix
    * @param rhs right hand side input and x output
    */
   void solve(DArray<8> &N_eff,DArray<5> &rhs){

      int n = N_eff.shape(0) * N_eff.shape(1) * N_eff.shape(2) * N_eff.shape(3);

      int *ipiv = new int [n];

      lapack::getrf(CblasRowMajor,n,n, N_eff.data(), n,ipiv);

      lapack::getrs(CblasRowMajor,'N',n,d, N_eff.data(),n,ipiv, rhs.data(),d);

      delete [] ipiv;

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

      enum {i,j,k,l,m,n,o};

      if(row == 0){

         DArray<5> rhs;
         DArray<8> N_eff;

         //first make left and right intermediary objects
         DArray<7> LI7;
         DArray<7> RI7;

         DArray<6> lop;
         DArray<6> rop;

         if(col == 0){

            Gemm(CblasNoTrans,CblasNoTrans,1.0,env.gt(0)[0],R,0.0,RI7);

            //act with operators on left and right peps
            Contract(1.0,peps(row,col),shape(i,j,k,l,m),global::trot.gLO_n(),shape(k,o,n),0.0,lop,shape(i,j,n,o,l,m));
            Contract(1.0,peps(row+1,col),shape(i,j,k,l,m),global::trot.gRO_n(),shape(k,o,n),0.0,rop,shape(i,j,n,o,l,m));

            cout << cost_function_vertical(row,col,peps,lop,rop,L,R,LI7,RI7) << endl;

            //start sweeping
            int iter = 0;

            while(iter < n_iter){

               // --(1)-- top site

               //construct effective environment and right hand side for linear system of top site
               construct_lin_sys_vertical(row,col,peps,lop,rop,L,R,N_eff,rhs,LI7,RI7,true);

               //solve the system
               solve(N_eff,rhs);

               //update upper peps
               Permute(rhs,shape(1,2,0,3,4),peps(row+1,col));

               // --(2)-- bottom site

               //construct effective environment and right hand side for linear system of bottom site
               construct_lin_sys_vertical(row,col,peps,lop,rop,L,R,N_eff,rhs,LI7,RI7,false);

               //solve the system
               solve(N_eff,rhs);

               //update lower peps
               Permute(rhs,shape(1,2,0,3,4),peps(row,col));

               //repeat until converged
               ++iter;

            }

         }
         else{//col != 0

         }

      }
      else{//whatever other options, row != 0

      }

   }

   /**
    * construct the single-site effective environment and right hand side needed for the linear system of the vertical gate, for top or bottom site
    * @param row , the row index of the bottom site
    * @param col column index of the vertical column
    * @param peps, full PEPS object before update
    * @param L Left environment contraction
    * @param R Right environment contraction
    * @param N_eff output object, contains N_eff on output
    * @param rhs output object, contains N_eff on output
    * @param top boolean flag for top or bottom site of vertical gate
    */
   void construct_lin_sys_vertical(int row,int col,PEPS<double> &peps,const DArray<6> &lop,const DArray<6> &rop,const DArray<5> &L,

         const DArray<5> &R,DArray<8> &N_eff,DArray<5> &rhs,const DArray<7> &LI7,const DArray<7> &RI7,bool top){

      if(row == 0){

         if(col == 0){

            if(top){//top site of vertical, so site (row,col+1) environment

               //paste bottom peps to right intermediate
               DArray<10> tmp10;
               Gemm(CblasNoTrans,CblasTrans,1.0,RI7,peps(row,col),0.0,tmp10);

               //construct right hand side, attach operator to tmp10:
               DArray<7> tmp7 = tmp10.reshape_clear( shape(D,D,D,D,D,D,d) );

               //another bottom peps to this one
               DArray<8> tmp8;
               Contract(1.0,tmp7,shape(6,4),peps(row,col),shape(2,4),0.0,tmp8);

               Permute(tmp8,shape(5,0,6,2,7,1,4,3),N_eff);

               //right hand side: add left operator to tmp7
               DArray<9> tmp9;
               Contract(1.0,tmp7,shape(6,4),lop,shape(2,5),0.0,tmp9);

               //remove the dimension-one legs
               tmp7 = tmp9.reshape_clear( shape(D,D,D,D,D,D,global::trot.gLO_n().shape(1)) );

               DArray<5> tmp5;
               Contract(1.0,tmp7,shape(0,6,5,2),rop,shape(1,3,4,5),0.0,tmp5);

               rhs.clear();
               Permute(tmp5,shape(4,3,0,2,1),rhs);

            }
            else{//bottom site (row,col)
               
               //paste top peps to right intermediate
               DArray<6> tmp6;
               Contract(1.0,RI7,shape(0,2,4),peps(row+1,col),shape(0,1,4),0.0,tmp6);

               //add another top peps for N_eff
               DArray<5> tmp5;
               Contract(1.0,tmp6,shape(0,4,1),peps(row+1,col),shape(1,2,4),0.0,tmp5);

               DArray<5> tmp5bis;
               Permute(tmp5,shape(3,4,0,2,1),tmp5bis);

               N_eff = tmp5bis.reshape_clear( shape(1,D,1,D,1,D,1,D) );

               //right hand side
               DArray<6> tmp6bis;
               Contract(1.0,tmp6,shape(0,4,1),rop,shape(1,2,5),0.0,tmp6bis);

               DArray<4> tmp4;
               Contract(1.0,tmp6bis,shape(3,5,4,0),lop,shape(0,1,3,5),0.0,tmp4);

               DArray<4> tmp4bis;
               Permute(tmp4,shape(2,1,0,3),tmp4bis);

               rhs = tmp4bis.reshape_clear( shape(d,1,D,1,D) );

            }

         }
         else{//later

         }

      }
      else{//later

      }

   }

   /**
    * evaluate the cost function of the linear system for two vertically connected PEPS: -2 <\Psi|\Psi'> + <\Psi'|\Psi'> where \Psi full PEPS with operator and \Psi is old PEPS
    * @param row , the row index of the bottom site
    * @param col column index of the vertical column
    * @param peps, full PEPS object before update
    * @param L Left environment contraction
    * @param R Right environment contraction
    */
   double cost_function_vertical(int row,int col,PEPS<double> &peps,const DArray<6> &lop,const DArray<6> &rop,const DArray<5> &L,const DArray<5> &R,const DArray<7> &LI7,const DArray<7> &RI7){

      if(row == 0){

         if(col == 0){

            // --- (1) calculate overlap of approximated state:

            //paste bottom peps to right intermediate
            DArray<10> tmp10;
            Gemm(CblasNoTrans,CblasTrans,1.0,RI7,peps(row,col),0.0,tmp10);

            DArray<7> tmp7 = tmp10.reshape_clear( shape(D,D,D,D,D,D,d) );

            //another bottom peps to this one
            DArray<8> tmp8;
            Contract(1.0,tmp7,shape(6,4),peps(row,col),shape(2,4),0.0,tmp8);

            DArray<6> tmp6 = tmp8.reshape_clear( shape(D,D,D,D,D,D) );

            //add upper peps
            DArray<5> tmp5;
            Contract(1.0,tmp6,shape(0,4,2),peps(row+1,col),shape(1,3,4),0.0,tmp5);

            DArray<5> tmp5bis;
            Permute(tmp5,shape(3,0,4,2,1),tmp5bis);

            double val = Dot(tmp5bis,peps(row+1,col));
            cout << val << endl;

            cout << val << endl;

            // --- (2) calculate 'b' part of overlap
            
            //right hand side: add left operator to tmp7
            DArray<9> tmp9;
            Contract(1.0,tmp7,shape(6,4),lop,shape(2,5),0.0,tmp9);

            //remove the dimension-one legs
            tmp7 = tmp9.reshape_clear( shape(D,D,D,D,D,D,global::trot.gLO_n().shape(1)) );

            tmp5.clear();
            Contract(1.0,tmp7,shape(0,6,5,2),rop,shape(1,3,4,5),0.0,tmp5);

            tmp5bis.clear();
            Permute(tmp5,shape(3,0,4,2,1),tmp5bis);

            val -= 2.0 * Dot(tmp5bis,peps(row+1,col));

            return val;

         }
         else{//col != 0

         }

      }
      else{//row != 0

      }

      return 0.0;

   }

}
