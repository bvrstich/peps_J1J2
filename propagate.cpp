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

      for(int col = 0;col < Lx - 1;++col){

         // --- (1) update the vertical pair on column 'col' ---
         update_vertical(0,col,peps,L,R[col],n_sweeps); 

         // --- (2) update the horizontal pair on column 'col'-'col+1' ---
         update_horizontal(0,col,peps,L,R[col+1],n_sweeps); 

         // --- (3) update diagonal LU-RD
         // todo

         // --- (4) update diagonal LD-RU
         // todo

         contractions::update_L('b',col,peps,L);

      }

      //one last vertical update
      update_vertical(0,Lx-1,peps,L,R[Lx-2],n_sweeps); 

      //update the bottom row for the new peps
      env.gb(0).fill('b',peps);

      // ---------------------------------------------------//
      // --- !!! (2) the middle rows (1 -> Ly-2) (2) !!! ---// 
      // ---------------------------------------------------//

      //renormalized operators for the middle sites
      vector< DArray<6> > RO(Lx - 1);
      DArray<6> LO;

      //for(int row = 1;row < Ly-1;++row){
      int row = 1;

         //first create right renormalized operator
         contractions::init_ro(row,peps,RO);


         for(int col = 0;col < 2;++col){

            // --- (1) update vertical pair on column 'col', with lowest site on row 'row'
            update_vertical(row,col,peps,LO,RO[col],n_sweeps); 


            // --- (2) update the horizontal pair on column 'col'-'col+1' ---
            update_horizontal(row,col,peps,LO,RO[col+1],n_sweeps); 

            // --- (3) update diagonal LU-RD
            // todo

            // --- (4) update diagonal LD-RU
            // todo


            //first construct a double layer object for the newly updated bottom 
            contractions::update_L(row,col,peps,LO);

         }
         
         /*

         //finally update the 'bottom' environment for the row
         env.add_layer('b',row,peps);

         //    }

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
    * update the tensors in a sweeping fashion, for bottom or top rows, i.e. with R and L environments of order 5
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

            Gemm(CblasNoTrans,CblasNoTrans,1.0,env.gt(0)[col],R,0.0,RI7);

            //act with operators on left and right peps
            Contract(1.0,peps(row,col),shape(i,j,k,l,m),global::trot.gLO_n(),shape(k,o,n),0.0,lop,shape(i,j,n,o,l,m));
            Contract(1.0,peps(row+1,col),shape(i,j,k,l,m),global::trot.gRO_n(),shape(k,o,n),0.0,rop,shape(i,j,n,o,l,m));

            //start sweeping
            int iter = 0;

            while(iter < n_iter){

               // --(1)-- top site

               //construct effective environment and right hand side for linear system of top site
               construct_lin_sys_vertical(row,col,peps,lop,rop,N_eff,rhs,L,R,LI7,RI7,true);

               //solve the system
               solve(N_eff,rhs);

               //update upper peps
               Permute(rhs,shape(0,1,4,2,3),peps(row+1,col));

               // --(2)-- bottom site

               //construct effective environment and right hand side for linear system of bottom site
               construct_lin_sys_vertical(row,col,peps,lop,rop,N_eff,rhs,L,R,LI7,RI7,false);

               //solve the system
               solve(N_eff,rhs);

               //update lower peps
               Permute(rhs,shape(0,1,4,2,3),peps(row,col));

               //repeat until converged
               ++iter;

            }

         }
         else if(col < Lx - 1){//col != 0

            Gemm(CblasTrans,CblasNoTrans,1.0,L,env.gt(0)[col],0.0,LI7);

            //act with operators on left and right peps
            Contract(1.0,peps(row,col),shape(i,j,k,l,m),global::trot.gLO_n(),shape(k,o,n),0.0,lop,shape(i,j,n,o,l,m));
            Contract(1.0,peps(row+1,col),shape(i,j,k,l,m),global::trot.gRO_n(),shape(k,o,n),0.0,rop,shape(i,j,n,o,l,m));

            //start sweeping
            int iter = 0;

            while(iter < n_iter){

               // --(1)-- top site

               //construct effective environment and right hand side for linear system of top site
               construct_lin_sys_vertical(row,col,peps,lop,rop,N_eff,rhs,L,R,LI7,RI7,true);

               //solve the system
               solve(N_eff,rhs);

               //update upper peps
               Permute(rhs,shape(0,1,4,2,3),peps(row+1,col));

               // --(2)-- bottom site

               //construct effective environment and right hand side for linear system of bottom site
               construct_lin_sys_vertical(row,col,peps,lop,rop,N_eff,rhs,L,R,LI7,RI7,false);

               //solve the system
               solve(N_eff,rhs);

               //update lower peps
               Permute(rhs,shape(0,1,4,2,3),peps(row,col));

               //repeat until converged
               ++iter;

            }

         }
         else{//col == Lx-1

            Gemm(CblasTrans,CblasNoTrans,1.0,L,env.gt(0)[col],0.0,LI7);

            Contract(1.0,peps(row,col),shape(i,j,k,l,m),global::trot.gLO_n(),shape(k,o,n),0.0,lop,shape(i,j,n,o,l,m));
            Contract(1.0,peps(row+1,col),shape(i,j,k,l,m),global::trot.gRO_n(),shape(k,o,n),0.0,rop,shape(i,j,n,o,l,m));

            //start sweeping
            int iter = 0;

            while(iter < n_iter){

               // --(1)-- top site

               //construct effective environment and right hand side for linear system of top site
               construct_lin_sys_vertical(row,col,peps,lop,rop,N_eff,rhs,L,R,LI7,RI7,true);

               //solve the system
               solve(N_eff,rhs);

               //update upper peps
               Permute(rhs,shape(0,1,4,2,3),peps(row+1,col));

               // --(2)-- bottom site

               //construct effective environment and right hand side for linear system of bottom site
               construct_lin_sys_vertical(row,col,peps,lop,rop,N_eff,rhs,L,R,LI7,RI7,false);

               //solve the system
               solve(N_eff,rhs);

               //update lower peps
               Permute(rhs,shape(0,1,4,2,3),peps(row,col));

               //repeat until converged
               ++iter;

            }

         }

      }
      else{//row == Lx - 2

      }

   }

   /**
    * update the tensors in a sweeping fashion, for middle rows, i.e. with R and L environments of order 6
    * @param row , the row index of the bottom site
    * @param col column index of the vertical column
    * @param peps, full PEPS object before update
    * @param LO Left environment contraction
    * @param RO Right environment contraction
    * @param n_iter nr of sweeps in the ALS algorithm
    */
   void update_vertical(int row,int col,PEPS<double> &peps,const DArray<6> &LO,const DArray<6> &RO,int n_iter){

      enum {i,j,k,l,m,n,o};

      DArray<5> rhs;
      DArray<8> N_eff;

      //first make left and right intermediary objects
      DArray<8> LI8;
      DArray<8> RI8;

      DArray<6> lop;
      DArray<6> rop;

      if(col == 0){

         DArray<8> tmp8;
         Gemm(CblasNoTrans,CblasNoTrans,1.0,env.gt(row)[col],RO,0.0,tmp8);

         //inefficient but it's on the side, so it doesn't matter
         Contract(1.0,tmp8,shape(0,7),env.gb(row-1)[col],shape(0,3),0.0,RI8);

         //act with operators on left and right peps
         Contract(1.0,peps(row,col),shape(i,j,k,l,m),global::trot.gLO_n(),shape(k,o,n),0.0,lop,shape(i,j,n,o,l,m));
         Contract(1.0,peps(row+1,col),shape(i,j,k,l,m),global::trot.gRO_n(),shape(k,o,n),0.0,rop,shape(i,j,n,o,l,m));

         //start sweeping
         int iter = 0;

         while(iter < n_iter){

            // --(1)-- top site

            //construct effective environment and right hand side for linear system of top site
            construct_lin_sys_vertical(row,col,peps,lop,rop,N_eff,rhs,LO,RO,LI8,RI8,true);

            //solve the system
            solve(N_eff,rhs);

            //update upper peps
            Permute(rhs,shape(0,1,4,2,3),peps(row+1,col));

            // --(2)-- bottom site

            //construct effective environment and right hand side for linear system of bottom site
            construct_lin_sys_vertical(row,col,peps,lop,rop,N_eff,rhs,LO,RO,LI8,RI8,false);

            //solve the system
            solve(N_eff,rhs);

            //update lower peps
            Permute(rhs,shape(0,1,4,2,3),peps(row,col));

            //repeat until converged
            ++iter;

         }

      }
      else{//col != 0

         //RI8
         Gemm(CblasNoTrans,CblasTrans,1.0,env.gb(row-1)[col],RO,0.0,RI8);

         //LI8
         Gemm(CblasTrans,CblasNoTrans,1.0,LO,env.gt(row)[col],0.0,LI8);

         //act with operators on left and right peps
         Contract(1.0,peps(row,col),shape(i,j,k,l,m),global::trot.gLO_n(),shape(k,o,n),0.0,lop,shape(i,j,n,o,l,m));
         Contract(1.0,peps(row+1,col),shape(i,j,k,l,m),global::trot.gRO_n(),shape(k,o,n),0.0,rop,shape(i,j,n,o,l,m));

         //start sweeping
         int iter = 0;

         while(iter < n_iter){

            // --(1)-- top site
            cout << iter << "\t" << cost_function_vertical(row,col,peps,lop,rop,LO,RO,LI8,RI8) << endl;

            //construct effective environment and right hand side for linear system of top site
            construct_lin_sys_vertical(row,col,peps,lop,rop,N_eff,rhs,LO,RO,LI8,RI8,true);

            //solve the system
            solve(N_eff,rhs);

            //update upper peps
            Permute(rhs,shape(0,1,4,2,3),peps(row+1,col));

            // --(2)-- bottom site

            //construct effective environment and right hand side for linear system of bottom site
            construct_lin_sys_vertical(row,col,peps,lop,rop,N_eff,rhs,LO,RO,LI8,RI8,false);

            //solve the system
            solve(N_eff,rhs);

            //update lower peps
            Permute(rhs,shape(0,1,4,2,3),peps(row,col));

            //repeat until converged
            ++iter;

         }

      }

   }

   /**
    * update the tensors in a sweeping fashion, horizontal pairs on the bottom or top rows (with R and L of order 5)
    * @param row , the row index of the horizontal row
    * @param col column index of left site
    * @param peps, full PEPS object before update
    * @param L Left environment contraction
    * @param R Right environment contraction
    * @param n_iter nr of sweeps in the ALS algorithm
    */
   void update_horizontal(int row,int col,PEPS<double> &peps,const DArray<5> &L,const DArray<5> &R,int n_iter){

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

            //create left and right intermediary operators: right
            Gemm(CblasNoTrans,CblasNoTrans,1.0,env.gt(0)[col+1],R,0.0,RI7);

            DArray<8> tmp8;
            Contract(1.0,RI7,shape(1,3),peps(row+1,col+1),shape(1,4),0.0,tmp8);

            DArray<7> tmp7;
            Contract(1.0,tmp8,shape(1,6,2),peps(row+1,col+1),shape(1,2,4),0.0,tmp7);

            Permute(tmp7,shape(0,3,5,4,6,1,2),RI7);

            //left
            DArray<5> tmp5;
            Gemm(CblasTrans,CblasNoTrans,1.0,env.gt(0)[col],peps(row+1,col),0.0,tmp5);

            DArray<6> tmp6;
            Contract(1.0,tmp5,shape(0,2),peps(row+1,col),shape(1,2),0.0,tmp6);

            DArray<6> tmp6bis;
            Permute(tmp6,shape(0,2,5,1,4,3),tmp6bis);

            LI7 = tmp6bis.reshape_clear( shape(env.gt(0)[col].shape(3),D,D,D,D,1,1) );

            //act with operators on left and right peps
            Contract(1.0,peps(row,col),shape(i,j,k,l,m),global::trot.gLO_n(),shape(k,o,n),0.0,lop,shape(i,j,n,o,l,m));
            Contract(1.0,peps(row,col+1),shape(i,j,k,l,m),global::trot.gRO_n(),shape(k,o,n),0.0,rop,shape(i,j,n,o,l,m));

            //start sweeping
            int iter = 0;

            while(iter < n_iter){

               // --(1)-- left site

               //construct effective environment and right hand side for linear system of left site
               construct_lin_sys_horizontal(row,col,peps,lop,rop,N_eff,rhs,L,R,LI7,RI7,true);

               //solve the system
               solve(N_eff,rhs);

               //update upper peps
               Permute(rhs,shape(0,1,4,2,3),peps(row,col));

               // --(2)-- right site

               //construct effective environment and right hand side for linear system of right site
               construct_lin_sys_horizontal(row,col,peps,lop,rop,N_eff,rhs,L,R,LI7,RI7,false);

               //solve the system
               solve(N_eff,rhs);

               //update lower peps
               Permute(rhs,shape(0,1,4,2,3),peps(row,col+1));

               //repeat until converged
               ++iter;

            }

         }
         else if(col < Lx - 2){//col != 0

            //create left and right intermediary operators: right
            Gemm(CblasNoTrans,CblasNoTrans,1.0,env.gt(0)[col+1],R,0.0,RI7);

            DArray<8> tmp8;
            Contract(1.0,RI7,shape(1,3),peps(row+1,col+1),shape(1,4),0.0,tmp8);

            DArray<7> tmp7;
            Contract(1.0,tmp8,shape(1,6,2),peps(row+1,col+1),shape(1,2,4),0.0,tmp7);

            RI7.clear();
            Permute(tmp7,shape(0,3,5,4,6,1,2),RI7);

            //left
            Gemm(CblasTrans,CblasNoTrans,1.0,L,env.gt(0)[col],0.0,LI7);

            tmp8.clear();
            Contract(1.0,LI7,shape(0,4),peps(row+1,col),shape(0,1),0.0,tmp8);

            tmp7.clear();
            Contract(1.0,tmp8,shape(0,3,5),peps(row+1,col),shape(0,1,2),0.0,tmp7);

            LI7.clear();
            Permute(tmp7,shape(2,4,6,3,5,0,1),LI7);

            //act with operators on left and right peps
            Contract(1.0,peps(row,col),shape(i,j,k,l,m),global::trot.gLO_n(),shape(k,o,n),0.0,lop,shape(i,j,n,o,l,m));
            Contract(1.0,peps(row,col+1),shape(i,j,k,l,m),global::trot.gRO_n(),shape(k,o,n),0.0,rop,shape(i,j,n,o,l,m));

            //start sweeping
            int iter = 0;

            while(iter < n_iter){

               // --(1)-- left site

               //construct effective environment and right hand side for linear system of left site
               construct_lin_sys_horizontal(row,col,peps,lop,rop,N_eff,rhs,L,R,LI7,RI7,true);

               //solve the system
               solve(N_eff,rhs);

               //update upper peps
               Permute(rhs,shape(0,1,4,2,3),peps(row,col));

               // --(2)-- right site

               //construct effective environment and right hand side for linear system of right site
               construct_lin_sys_horizontal(row,col,peps,lop,rop,N_eff,rhs,L,R,LI7,RI7,false);

               //solve the system
               solve(N_eff,rhs);

               //update lower peps
               Permute(rhs,shape(0,1,4,2,3),peps(row,col+1));

               //repeat until converged
               ++iter;

            }

         }
         else{//col == Lx - 2

            //create left and right intermediary operators: right
            DArray<5> tmp5;
            Contract(1.0,env.gt(0)[Lx - 1],shape(2,3),peps(1,Lx - 1),shape(1,4),0.0,tmp5);

            DArray<6> tmp6;
            Contract(1.0,tmp5,shape(1,3),peps(1,Lx - 1),shape(1,2),0.0,tmp6);

            DArray<6> tmp6bis;
            Permute(tmp6,shape(0,3,1,4,2,5),tmp6bis);

            RI7 = tmp6bis.reshape_clear( shape(env.gt(0)[Lx - 1].shape(0),D,D,D,D,1,1) );

            //left
            Gemm(CblasTrans,CblasNoTrans,1.0,L,env.gt(0)[col],0.0,LI7);

            DArray<8> tmp8;
            Contract(1.0,LI7,shape(0,4),peps(row+1,col),shape(0,1),0.0,tmp8);

            DArray<7> tmp7;
            Contract(1.0,tmp8,shape(0,3,5),peps(row+1,col),shape(0,1,2),0.0,tmp7);

            LI7.clear();
            Permute(tmp7,shape(2,4,6,3,5,0,1),LI7);

            //act with operators on left and right peps
            Contract(1.0,peps(row,col),shape(i,j,k,l,m),global::trot.gLO_n(),shape(k,o,n),0.0,lop,shape(i,j,n,o,l,m));
            Contract(1.0,peps(row,col+1),shape(i,j,k,l,m),global::trot.gRO_n(),shape(k,o,n),0.0,rop,shape(i,j,n,o,l,m));

            //start sweeping
            int iter = 0;

            while(iter < n_iter){

               // --(1)-- left site

               //construct effective environment and right hand side for linear system of left site
               construct_lin_sys_horizontal(row,col,peps,lop,rop,N_eff,rhs,L,R,LI7,RI7,true);

               //solve the system
               solve(N_eff,rhs);

               //update upper peps
               Permute(rhs,shape(0,1,4,2,3),peps(row,col));

               // --(2)-- right site

               //construct effective environment and right hand side for linear system of right site
               construct_lin_sys_horizontal(row,col,peps,lop,rop,N_eff,rhs,L,R,LI7,RI7,false);

               //solve the system
               solve(N_eff,rhs);

               //update lower peps
               Permute(rhs,shape(0,1,4,2,3),peps(row,col+1));

               //repeat until converged
               ++iter;

            }

         }

      }
      else{//whatever other options, row != 0

      }

   }


   /**
    * update the tensors in a sweeping fashion, horizontal pairs on the middle rows (with R and L of order 6)
    * @param row , the row index of the horizontal row
    * @param col column index of left site
    * @param peps, full PEPS object before update
    * @param LO Left environment contraction
    * @param RO Right environment contraction
    * @param n_iter nr of sweeps in the ALS algorithm
    */
   void update_horizontal(int row,int col,PEPS<double> &peps,const DArray<6> &LO,const DArray<6> &RO,int n_iter){

      enum {i,j,k,l,m,n,o};

      DArray<5> rhs;
      DArray<8> N_eff;

      //first make left and right intermediary objects
      DArray<8> LI8;
      DArray<8> RI8;

      DArray<6> lop;
      DArray<6> rop;

      if(col == 0){

         //create left and right intermediary operators: right
         Gemm(CblasNoTrans,CblasNoTrans,1.0,env.gt(row)[col+1],RO,0.0,RI8);

         DArray<9> tmp9;
         Contract(1.0,RI8,shape(1,3),peps(row+1,col+1),shape(1,4),0.0,tmp9);

         DArray<8> tmp8;
         Contract(1.0,tmp9,shape(1,7,2),peps(row+1,col+1),shape(1,2,4),0.0,tmp8);

         RI8.clear();
         Permute(tmp8,shape(0,4,6,5,7,1,2,3),RI8);

         //left
         DArray<5> tmp5;
         Gemm(CblasTrans,CblasNoTrans,1.0,env.gt(row)[col],peps(row+1,col),0.0,tmp5);

         DArray<6> tmp6;
         Contract(1.0,tmp5,shape(0,2),peps(row+1,col),shape(1,2),0.0,tmp6);

         DArray<6> tmp6bis;
         Permute(tmp6,shape(0,2,5,1,4,3),tmp6bis);

         LI8 = tmp6bis.reshape_clear( shape(env.gt(0)[col].shape(3),D,D,D,D,1,1,1) );

         //act with operators on left and right peps
         Contract(1.0,peps(row,col),shape(i,j,k,l,m),global::trot.gLO_n(),shape(k,o,n),0.0,lop,shape(i,j,n,o,l,m));
         Contract(1.0,peps(row,col+1),shape(i,j,k,l,m),global::trot.gRO_n(),shape(k,o,n),0.0,rop,shape(i,j,n,o,l,m));

         //start sweeping
         int iter = 0;

         while(iter < n_iter){

            // --(1)-- left site

            //construct effective environment and right hand side for linear system of left site
            construct_lin_sys_horizontal(row,col,peps,lop,rop,N_eff,rhs,LO,RO,LI8,RI8,true);

            //solve the system
            solve(N_eff,rhs);

            //update upper peps
            Permute(rhs,shape(0,1,4,2,3),peps(row,col));

            // --(2)-- right site

            //construct effective environment and right hand side for linear system of right site
            construct_lin_sys_horizontal(row,col,peps,lop,rop,N_eff,rhs,LO,RO,LI8,RI8,false);

            //solve the system
            solve(N_eff,rhs);

            //update lower peps
            Permute(rhs,shape(0,1,4,2,3),peps(row,col+1));

            //repeat until converged
            ++iter;

         }
      }

   }

   /**
    * construct the single-site effective environment and right hand side needed for the linear system of the vertical gate, for top or bottom site on top or bottom rows
    * @param row the row index of the bottom site
    * @param col column index of the vertical column
    * @param peps, full PEPS object before update
    * @param N_eff output object, contains N_eff on output
    * @param rhs output object, contains N_eff on output
    * @param L Left environment contraction
    * @param R Right environment contraction
    * @param LI7 left intermediate object
    * @param RI7 right intermediate object
    * @param top boolean flag for top or bottom site of vertical gate
    */
   void construct_lin_sys_vertical(int row,int col,PEPS<double> &peps,const DArray<6> &lop,const DArray<6> &rop,DArray<8> &N_eff,DArray<5> &rhs,

         const DArray<5> &L, const DArray<5> &R, const DArray<7> &LI7,const DArray<7> &RI7,bool top){

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
               Permute(tmp5,shape(3,0,2,1,4),rhs);

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
               Permute(tmp4,shape(1,3,0,2),tmp4bis);

               rhs = tmp4bis.reshape_clear( shape(1,D,1,D,d) );

            }

         }
         else if(col < Lx - 1){//col != 0

            if(top){//top site

               // (1) calculate N_eff

               //paste bottom peps to right
               DArray<8> tmp8;
               Gemm(CblasNoTrans,CblasTrans,1.0,peps(row,col),R,0.0,tmp8);

               //and another!
               DArray<7> tmp7;
               Contract(1.0,peps(row,col),shape(2,3,4),tmp8,shape(2,3,7),0.0,tmp7);

               //add to LI7
               DArray<8> tmp8bis;
               Contract(1.0,LI7,shape(2,3,6),tmp7,shape(0,2,4),0.0,tmp8bis);

               N_eff.clear();
               Permute(tmp8bis,shape(0,2,4,6,1,3,5,7),N_eff);

               // (2) right hand side

               //attach left operator to tmp8
               tmp8bis.clear();
               Contract(1.0,lop,shape(2,4,5),tmp8,shape(2,3,7),0.0,tmp8bis);

               //and right operator
               tmp8.clear();
               Contract(1.0,rop,shape(3,4,5),tmp8bis,shape(2,1,6),0.0,tmp8);

               DArray<5> tmp5;
               Contract(1.0,LI7,shape(6,2,3,0,4),tmp8,shape(6,3,4,0,1),0.0,tmp5);

               rhs.clear();
               Permute(tmp5,shape(0,1,3,4,2),rhs);

            }
            else{//bottom site

               //(1) first N_eff

               //paste top peps to left
               DArray<8> tmp8;
               Contract(1.0,LI7,shape(1,5),peps(row+1,col),shape(0,1),0.0,tmp8);

               //and another: watch out, order is reversed!
               DArray<7> tmp7;
               Contract(1.0,tmp8,shape(0,3,5),peps(row+1,col),shape(0,1,2),0.0,tmp7);

               //now add right side to it
               DArray<6> tmp6;
               Contract(1.0,tmp7,shape(2,6,4),R,shape(0,1,2),0.0,tmp6);

               DArray<6> tmp6bis;
               Permute(tmp6,shape(0,3,4,1,2,5),tmp6bis);

               N_eff = tmp6bis.reshape_clear(shape(D,D,1,D,D,D,1,D));

               // (2) right hand side

               //add right operator to tmp8
               DArray<8> tmp8bis;
               Contract(1.0,tmp8,shape(0,3,5),rop,shape(0,1,2),0.0,tmp8bis);

               //add left operator
               tmp8.clear();
               Contract(1.0,tmp8bis,shape(0,6,5),lop,shape(0,1,3),0.0,tmp8);

               //finally contract with right side
               DArray<5> tmp5;
               Contract(1.0,tmp8,shape(1,4,3,7),R,shape(0,1,2,3),0.0,tmp5);

               rhs.clear();
               Permute(tmp5,shape(0,1,3,4,2),rhs);

            }

         }
         else{ //col == Lx - 1

            if(top){//top site

               // (1) calculate N_eff

               //paste bottom peps to left
               DArray<8> tmp8;
               Contract(1.0,LI7,shape(3,6),peps(row,col),shape(0,4),0.0,tmp8);

               //and another
               DArray<7> tmp7;
               Contract(1.0,tmp8,shape(2,6,7),peps(row,col),shape(0,2,3),0.0,tmp7);

               DArray<7> tmp7bis;
               Permute(tmp7,shape(0,2,5,1,3,4,6),tmp7bis);

               N_eff = tmp7bis.reshape_clear(shape(D,D,D,1,D,D,D,1));

               // (2) right hand side

               //attach left operator to tmp8
               DArray<8> tmp8bis;
               Contract(1.0,tmp8,shape(2,6,7),lop,shape(0,2,4),0.0,tmp8bis);

               //and right operator
               DArray<4> tmp4;
               Contract(1.0,tmp8bis,shape(0,2,6,5,7),rop,shape(0,1,3,4,5),0.0,tmp4);

               rhs = tmp4.reshape_clear( shape(D,D,D,1,d) );

            }
            else{//bottom site

               //(1) first N_eff

               //paste top peps to left
               DArray<8> tmp8;
               Contract(1.0,LI7,shape(1,5),peps(row+1,col),shape(0,1),0.0,tmp8);

               //and another: watch out, order is reversed!
               DArray<7> tmp7;
               Contract(1.0,tmp8,shape(0,3,5),peps(row+1,col),shape(0,1,2),0.0,tmp7);

               DArray<7> tmp7bis;
               Permute(tmp7,shape(0,5,1,3,2,4,6),tmp7bis);

               N_eff = tmp7bis.reshape_clear( shape(D,D,1,1,D,D,1,1) );

               // (2) right hand side

               //add right operator to tmp8
               DArray<8> tmp8bis;
               Contract(1.0,tmp8,shape(0,3,5),rop,shape(0,1,2),0.0,tmp8bis);

               //add left operator
               tmp8.clear();
               Contract(1.0,tmp8bis,shape(0,6,5),lop,shape(0,1,3),0.0,tmp8);

               //finally contract with right side
               rhs = tmp8.reshape_clear( shape(D,D,1,1,d) );

            }

         }

      }
      else{//row != 0

      }

   }

   /**
    * construct the single-site effective environment and right hand side needed for the linear system of the vertical gate, for top or bottom site on middle rows
    * @param row the row index of the bottom site
    * @param col column index of the vertical column
    * @param peps, full PEPS object before update
    * @param N_eff output object, contains N_eff on output
    * @param rhs output object, contains N_eff on output
    * @param LO Left environment contraction
    * @param RO Right environment contraction
    * @param LI8 left intermediate object
    * @param RI8 right intermediate object
    * @param top boolean flag for top or bottom site of vertical gate
    */
   void construct_lin_sys_vertical(int row,int col,PEPS<double> &peps,const DArray<6> &lop,const DArray<6> &rop,DArray<8> &N_eff,DArray<5> &rhs,

         const DArray<6> &LO, const DArray<6> &RO, const DArray<8> &LI8,const DArray<8> &RI8,bool top){

      if(col == 0){

         if(top){//top site environment

            // (1) calculate N_eff

            //add bottom peps  to intermediate
            DArray<9> tmp9;
            Contract(1.0,peps(row,col),shape(3,4),RI8,shape(7,5),0.0,tmp9);

            //and another
            DArray<8> tmp8;
            Contract(1.0,peps(row,col),shape(2,3,4),tmp9,shape(2,8,7),0.0,tmp8);

            N_eff.clear();
            Permute(tmp8,shape(0,4,1,6,2,5,3,7),N_eff);

            // (2) right hand side

            //add left operator to intermediate
            DArray<9> tmp9bis;
            Contract(1.0,lop,shape(2,4,5),tmp9,shape(2,8,7),0.0,tmp9bis);

            //and right operator
            DArray<5> tmp5;
            Contract(1.0,rop,shape(0,1,3,4,5),tmp9bis,shape(0,5,2,1,7),0.0,tmp5);

            rhs.clear();
            Permute(tmp5,shape(1,3,2,4,0),rhs);

         }
         else{//bottom site

            // (1) calculate N_eff

            //add top to intermediate
            DArray<9> tmp9;
            Contract(1.0,peps(row+1,col),shape(1,4),RI8,shape(1,3),0.0,tmp9);

            //and another
            DArray<8> tmp8;
            Contract(1.0,peps(row+1,col),shape(1,2,4),tmp9,shape(3,1,4),0.0,tmp8);

            N_eff.clear();
            Permute(tmp8,shape(0,1,6,4,2,3,7,5),N_eff);

            // (2) right hand side

            //add right operator
            DArray<9> tmp9bis;
            Contract(1.0,rop,shape(1,2,5),tmp9,shape(3,1,4),0.0,tmp9bis);

            //and right operator
            DArray<5> tmp5;
            Contract(1.0,lop,shape(0,1,3,4,5),tmp9bis,shape(0,2,1,7,5),0.0,tmp5);

            rhs.clear();
            Permute(tmp5,shape(1,2,4,3,0),rhs);

         }

      }
      else{//col != 0
         
         if(top){//top site environment

            // (1) calculate N_eff

            //add bottom peps  to intermediate right
            DArray<9> tmp9;
            Contract(1.0,peps(row,col),shape(3,4),RI8,shape(2,7),0.0,tmp9);

            //and another
            DArray<8> tmp8;
            Contract(1.0,peps(row,col),shape(2,3,4),tmp9,shape(2,4,8),0.0,tmp8);

            //now contract with left side
            DArray<8> tmp8bis;
            Contract(1.0,LI8,shape(7,2,3,4),tmp8,shape(5,0,2,4),0.0,tmp8bis);
            
            N_eff.clear();
            Permute(tmp8bis,shape(0,2,4,6,1,3,5,7),N_eff);

            // (2) right hand side

            //add left operator to intermediate right
            DArray<9> tmp9bis;
            Contract(1.0,lop,shape(2,4,5),tmp9,shape(2,4,8),0.0,tmp9bis);

            //and right operator
            tmp9.clear();
            Contract(1.0,rop,shape(3,4,5),tmp9bis,shape(2,1,7),0.0,tmp9);

            //contract with left hand side
            DArray<5> tmp5;
            Contract(1.0,LI8,shape(7,5,0,2,3,4),tmp9,shape(7,1,0,3,4,6),0.0,tmp5);

            rhs.clear();
            Permute(tmp5,shape(0,1,3,4,2),rhs);

         }
         else{//bottom site

            // (1) calculate N_eff

            //add upper peps to LI8
            DArray<9> tmp9;
            Contract(1.0,LI8,shape(1,6),peps(row+1,col),shape(0,1),0.0,tmp9);

            //and another
            DArray<8> tmp8;
            Contract(1.0,tmp9,shape(0,4,6),peps(row+1,col),shape(0,1,2),0.0,tmp8);

            //contract with right intermediate
            DArray<8> tmp8bis;
            Contract(1.0,tmp8,shape(3,7,5,2),RI8,shape(3,4,5,0),0.0,tmp8bis);

            N_eff.clear();
            Permute(tmp8bis,shape(0,3,4,6,1,2,5,7),N_eff);

            // (2) right hand side

            //add right operator to intermediate
            DArray<9> tmp9bis;
            Contract(1.0,tmp9,shape(0,4,6),rop,shape(0,1,2),0.0,tmp9bis);

            //next add left operator
            tmp9.clear();
            Contract(1.0,tmp9bis,shape(0,7,6),lop,shape(0,1,3),0.0,tmp9);

            DArray<5> tmp5;
            Contract(1.0,tmp9,shape(2,5,4,8,7,1),RI8,shape(3,4,5,6,1,0),0.0,tmp5);

            rhs.clear();
            Permute(tmp5,shape(0,1,3,4,2),rhs);

         }

      }

   }


   /**
    * construct the single-site effective environment and right hand side needed for the linear system of the horinzontal gate, for left or right site on top ro bottom row ( with R,L order 5)
    * @param row , the row index of the bottom site
    * @param col column index of the vertical column
    * @param peps, full PEPS object before update
    * @param N_eff output object, contains N_eff on output
    * @param rhs output object, contains N_eff on output
    * @param L Left environment contraction
    * @param R Right environment contraction
    * @param LI7 Left intermediate object
    * @param RI7 Right intermediate object
    * @param left boolean flag for left (true) or right (false) site of horizontal gate
    */
   void construct_lin_sys_horizontal(int row,int col,PEPS<double> &peps,const DArray<6> &lop,const DArray<6> &rop,DArray<8> &N_eff,DArray<5> &rhs,

         const DArray<5> &L,const DArray<5> &R,const DArray<7> &LI7,const DArray<7> &RI7,bool left){

      if(row == 0){

         if(left){//left site of horizontal gate, so site (row,col) environment

            //(1) construct N_eff

            //add right peps to intermediate
            DArray<8> tmp8;
            Contract(1.0,RI7,shape(4,6),peps(row,col+1),shape(1,4),0.0,tmp8);

            //add second peps
            DArray<5> tmp5;
            Contract(1.0,tmp8,shape(3,6,7,4),peps(row,col+1),shape(1,2,3,4),0.0,tmp5);

            //contract with left hand side
            DArray<6> tmp6;
            Gemm(CblasTrans,CblasNoTrans,1.0,LI7,tmp5,0.0,tmp6);

            DArray<6> tmp6bis;
            Permute(tmp6,shape(2,0,5,3,1,4),tmp6bis);

            int DL = peps(row,col).shape(0);
            int DU = peps(row,col).shape(1);
            int DD = peps(row,col).shape(3);
            int DR = peps(row,col).shape(4);

            N_eff = tmp6bis.reshape_clear( shape(DL,DU,DD,DR,DL,DU,DD,DR) );

            // (2) construct right hand side

            //add right operator to tmp8
            tmp6.clear();
            Contract(1.0,tmp8,shape(3,6,7,4),rop,shape(1,2,4,5),0.0,tmp6);

            //attach LI7 to right side
            DArray<7> tmp7;
            Gemm(CblasTrans,CblasNoTrans,1.0,LI7,tmp6,0.0,tmp7);

            //now paste left operator in
            tmp5.clear();
            Contract(1.0,tmp7,shape(2,0,6,5),lop,shape(0,1,3,5),0.0,tmp5);

            rhs.clear();
            Permute(tmp5,shape(1,0,4,2,3),rhs);

         }
         else{//right site of horizontal gate, so site (row+1,col) environment

            //(1) constsruct N_eff

            //add left peps to LI7
            DArray<8> tmp8;
            Contract(1.0,LI7,shape(6,4),peps(row,col),shape(0,1),0.0,tmp8);

            //and another peps
            DArray<5> tmp5;
            Contract(1.0,tmp8,shape(4,3,5,6),peps(row,col),shape(0,1,2,3),0.0,tmp5);

            //contract with right hand side
            DArray<6> tmp6;
            Gemm(CblasTrans,CblasNoTrans,1.0,tmp5,RI7,0.0,tmp6);

            DArray<6> tmp6bis;
            Permute(tmp6,shape(1,2,4,0,3,5),tmp6bis);

            int DL = peps(row,col+1).shape(0);
            int DU = peps(row,col+1).shape(1);
            int DD = peps(row,col+1).shape(3);
            int DR = peps(row,col+1).shape(4);

            N_eff = tmp6bis.reshape_clear( shape(DL,DU,DD,DR,DL,DU,DD,DR) );

            // (2) construct right hand side

            //and another peps
            tmp6.clear();
            Contract(1.0,tmp8,shape(4,3,5,6),lop,shape(0,1,2,4),0.0,tmp6);

            //contract with RI7
            DArray<7> tmp7;
            Gemm(CblasTrans,CblasNoTrans,1.0,tmp6,RI7,0.0,tmp7);

            //now paste right operator in
            tmp5.clear();
            Contract(1.0,tmp7,shape(2,3,1,5),rop,shape(0,1,3,5),0.0,tmp5);

            rhs.clear();
            Permute(tmp5,shape(0,1,4,2,3),rhs);

         }

      }
      else{//row != 0

      }

   }

   /**
    * construct the single-site effective environment and right hand side needed for the linear system of the horinzontal gate, for left or right site on middle rows (with R,L order 6)
    * @param row , the row index of the bottom site
    * @param col column index of the vertical column
    * @param peps, full PEPS object before update
    * @param N_eff output object, contains N_eff on output
    * @param rhs output object, contains N_eff on output
    * @param LO Left environment contraction
    * @param RO Right environment contraction
    * @param LI8 left intermediate object
    * @param RI8 right intermediate object
    * @param left boolean flag for left (true) or right (false) site of horizontal gate
    */
   void construct_lin_sys_horizontal(int row,int col,PEPS<double> &peps,const DArray<6> &lop,const DArray<6> &rop,DArray<8> &N_eff,DArray<5> &rhs,

         const DArray<6> &LO,const DArray<6> &RO,const DArray<8> &LI8,const DArray<8> &RI8,bool left){

      if(left){//left site of horizontal gate, so site (row,col) environment

         //(1) construct N_eff

         //add right peps to intermediate
         DArray<9> tmp9;
         Contract(1.0,RI8,shape(4,6),peps(row,col+1),shape(1,4),0.0,tmp9);

         //add second peps
         DArray<8> tmp8;
         Contract(1.0,tmp9,shape(3,7,4),peps(row,col+1),shape(1,2,4),0.0,tmp8);

         //add bottom environment
         DArray<6> tmp6;
         Contract(1.0,tmp8,shape(7,5,3),env.gb(row-1)[col+1],shape(1,2,3),0.0,tmp6);

         //and next bottom environment
         tmp8.clear();
         Gemm(CblasNoTrans,CblasTrans,1.0,tmp6,env.gb(row-1)[col],0.0,tmp8);

         //contract with left hand side
         DArray<8> tmp8bis;
         Contract(1.0,LI8,shape(0,1,2,7),tmp8,shape(0,1,2,5),0.0,tmp8bis);

         N_eff.clear();
         Permute(tmp8bis,shape(2,0,6,5,3,1,7,4),N_eff);

         // (2) construct right hand side

         //add right operator
         DArray<9> tmp9bis;
         Contract(1.0,tmp9,shape(3,7,4),rop,shape(1,2,5),0.0,tmp9bis);

         //add bottom environment
         DArray<7> tmp7;
         Contract(1.0,tmp9bis,shape(8,5,3),env.gb(row-1)[col+1],shape(1,2,3),0.0,tmp7);

         //and next bottom environment
         tmp9.clear();
         Gemm(CblasNoTrans,CblasTrans,1.0,tmp7,env.gb(row-1)[col],0.0,tmp9);

         //now add left operator
         tmp9bis.clear();
         Contract(1.0,lop,shape(3,4,5),tmp9,shape(5,7,4),0.0,tmp9bis);

         //attach LI8 to right side
         DArray<5> tmp5;
         Contract(1.0,LI8,shape(0,1,2,3,5,7),tmp9bis,shape(3,4,5,1,0,7),0.0,tmp5);

         rhs.clear();
         Permute(tmp5,shape(1,0,4,3,2),rhs);

      }
      else{//right site of horizontal gate, so site (row+1,col) environment

         //(1) constsruct N_eff

         //add left peps to LI8
         DArray<9> tmp9;
         Contract(1.0,LI8,shape(6,4),peps(row,col),shape(0,1),0.0,tmp9);

         //and another peps
         DArray<8> tmp8;
         Contract(1.0,tmp9,shape(4,3,6),peps(row,col),shape(0,1,2),0.0,tmp8);

         //now contract with bottom environment
         DArray<6> tmp6;
         Contract(1.0,tmp8,shape(3,6,4),env.gb(row-1)[col],shape(0,1,2),0.0,tmp6);

         //and next bottom environment
         tmp8.clear();
         Gemm(CblasNoTrans,CblasNoTrans,1.0,tmp6,env.gb(row-1)[col+1],0.0,tmp8);

         //finally contract with RI8
         DArray<8> tmp8bis;
         Contract(1.0,tmp8,shape(0,1,2,7),RI8,shape(0,1,2,7),0.0,tmp8bis);

         N_eff.clear();
         Permute(tmp8bis,shape(1,4,2,6,0,5,3,7),N_eff);

         // (2) construct right hand side

         //add left operator
         DArray<9> tmp9bis;
         Contract(1.0,tmp9,shape(4,3,6),lop,shape(0,1,2),0.0,tmp9bis);

         //now contract with bottom environment
         DArray<7> tmp7;
         Contract(1.0,tmp9bis,shape(3,7,4),env.gb(row-1)[col],shape(0,1,2),0.0,tmp7);

         //and next bottom environment
         tmp9.clear();
         Gemm(CblasNoTrans,CblasNoTrans,1.0,tmp7,env.gb(row-1)[col+1],0.0,tmp9);

         //add right operator
         tmp9bis.clear();
         Contract(1.0,tmp9,shape(5,4,6),rop,shape(0,3,4),0.0,tmp9bis);

         //finally contract with RI8
         DArray<5> tmp5;
         Contract(1.0,tmp9bis,shape(0,1,2,6,8,5),RI8,shape(0,1,2,3,5,7),0.0,tmp5);

         rhs.clear();
         Permute(tmp5,shape(0,3,1,4,2),rhs);

      }

   }

   /**
    * evaluate the cost function of the linear system for two vertically connected PEPS: -2 <\Psi|\Psi'> + <\Psi'|\Psi'> where \Psi full PEPS with operator and \Psi is old PEPS
    * for top or bottom rows, i.e. with L and R of order 5 and intermediates of order 7
    * @param row , the row index of the bottom site
    * @param col column index of the vertical column
    * @param peps, full PEPS object before update
    * @param L Left environment contraction
    * @param R Right environment contraction
    * @param LI7 left intermediate object
    * @param RI7 left intermediate object
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
            Contract(1.0,tmp6,shape(0,5,2),peps(row+1,col),shape(1,3,4),0.0,tmp5);

            DArray<5> tmp5bis;
            Permute(tmp5,shape(3,0,4,2,1),tmp5bis);

            double val = Dot(tmp5bis,peps(row+1,col));

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
         else if(col < Lx - 1){//col != 0

            // --- (1) calculate overlap of approximated state:

            //add upper peps to LI7
            DArray<8> tmp8;
            Contract(1.0,LI7,shape(1,5),peps(row+1,col),shape(0,1),0.0,tmp8);

            //and another
            DArray<7> tmp7;
            Contract(1.0,tmp8,shape(0,3,5),peps(row+1,col),shape(0,1,2),0.0,tmp7);

            DArray<8> tmp8bis;
            Contract(1.0,tmp7,shape(0,5),peps(row,col),shape(0,1),0.0,tmp8bis);

            DArray<5> tmp5;
            Contract(1.0,tmp8bis,shape(0,2,5,6),peps(row,col),shape(0,1,2,3),0.0,tmp5);

            DArray<5> tmp5bis;
            Permute(tmp5,shape(0,2,1,3,4),tmp5bis);

            double val = Dot(tmp5bis,R);

            // --- (2) calculate 'b' part of overlap

            //add right operator to tmp8
            tmp8bis.clear();
            Contract(1.0,tmp8,shape(0,3,5),rop,shape(0,1,2),0.0,tmp8bis);

            //then add left operator
            tmp8.clear();
            Contract(1.0,tmp8bis,shape(0,6,5),lop,shape(0,1,3),0.0,tmp8);

            //finally add lop
            tmp5.clear();
            Contract(1.0,tmp8,shape(0,2,5,6),peps(row,col),shape(0,1,2,3),0.0,tmp5);

            tmp5bis.clear();
            Permute(tmp5,shape(0,2,1,3,4),tmp5bis);

            val -= 2.0 * Dot(tmp5bis,R);

            return val;

         }
         else{//col == Lx - 1

            DArray<8> tmp8;
            Contract(1.0,LI7,shape(1,5),peps(row+1,col),shape(0,1),0.0,tmp8);

            //and another
            DArray<7> tmp7;
            Contract(1.0,tmp8,shape(0,3,5),peps(row+1,col),shape(0,1,2),0.0,tmp7);

            DArray<8> tmp8bis;
            Contract(1.0,tmp7,shape(0,5),peps(row,col),shape(0,1),0.0,tmp8bis);

            DArray<5> tmp5 = tmp8bis.reshape_clear( shape(D,D,d,1,1) );

            double val =  Dot(tmp5,peps(row,col));

            // --- (2) calculate 'b' part of overlap

            //add right operator to tmp8
            tmp8bis.clear();
            Contract(1.0,tmp8,shape(0,3,5),rop,shape(0,1,2),0.0,tmp8bis);

            //then add left operator
            tmp8.clear();
            Contract(1.0,tmp8bis,shape(0,6,5),lop,shape(0,1,3),0.0,tmp8);

            //finally add lop
            tmp5 = tmp8.reshape_clear( shape(D,D,d,1,1) );

            val -= 2.0 * Dot(tmp5,peps(row,col));

            return val;

         }

      }
      else{//row != 0

      }

      return 0.0;

   }

   /**
    * evaluate the cost function of the linear system for two vertically connected PEPS: -2 <\Psi|\Psi'> + <\Psi'|\Psi'> where \Psi full PEPS with operator and \Psi is old PEPS
    * for middle rows, i.e. with R and L order 6 and intermediates of order 8
    * @param row , the row index of the bottom site
    * @param col column index of the vertical column
    * @param peps, full PEPS object before update
    * @param LO Left environment contraction
    * @param RO Right environment contraction
    * @param LI8 left intermediate object
    * @param RI8 right intermediate object
    */
   double cost_function_vertical(int row,int col,PEPS<double> &peps,const DArray<6> &lop,const DArray<6> &rop,const DArray<6> &LO,const DArray<6> &RO,const DArray<8> &LI8,const DArray<8> &RI8){

      // (1) calculate N_eff

      if(col == 0){

         //add bottom peps  to intermediate
         DArray<9> tmp9;
         Contract(1.0,peps(row,col),shape(3,4),RI8,shape(7,5),0.0,tmp9);

         //and another
         DArray<8> tmp8;
         Contract(1.0,peps(row,col),shape(2,3,4),tmp9,shape(2,8,7),0.0,tmp8);

         //and top one
         DArray<5> tmp5;
         Contract(1.0,peps(row+1,col),shape(0,1,3,4),tmp8,shape(0,4,1,6),0.0,tmp5);

         DArray<5> tmp5bis;
         Permute(tmp5,shape(1,3,0,2,4),tmp5bis);

         double val = Dot(tmp5bis,peps(row+1,col));

         // (2) right hand side

         //add left operator to intermediate
         DArray<9> tmp9bis;
         Contract(1.0,lop,shape(2,4,5),tmp9,shape(2,8,7),0.0,tmp9bis);

         //and right operator
         tmp5.clear();
         Contract(1.0,rop,shape(0,1,3,4,5),tmp9bis,shape(0,5,2,1,7),0.0,tmp5);

         tmp5bis.clear(); 
         Permute(tmp5,shape(1,3,0,2,4),tmp5bis);

         val -= 2.0 * Dot(tmp5bis,peps(row+1,col));

         return val;

      }
      else{//col != 0

         //add bottom peps  to intermediate right
         DArray<9> tmp9;
         Contract(1.0,peps(row,col),shape(3,4),RI8,shape(2,7),0.0,tmp9);

         //and another
         DArray<8> tmp8;
         Contract(1.0,peps(row,col),shape(2,3,4),tmp9,shape(2,4,8),0.0,tmp8);

         DArray<8> rn;
         Permute(tmp8,shape(5,6,7,1,3,0,2,4),rn);

         //add upper peps to LI8
         DArray<9> tmp9tris;
         Contract(1.0,LI8,shape(1,6),peps(row+1,col),shape(0,1),0.0,tmp9tris);

         //and another
         tmp8.clear();
         Contract(1.0,tmp9tris,shape(0,4,6),peps(row+1,col),shape(0,1,2),0.0,tmp8);

         DArray<8> ln;
         Permute(tmp8,shape(3,7,5,6,4,0,1,2),ln);

         double val = Dot(ln,rn);

         // (2) right hand side

         //add left operator to intermediate right
         DArray<9> tmp9bis;
         Contract(1.0,lop,shape(2,4,5),tmp9,shape(2,4,8),0.0,tmp9bis);

         //and right operator
         tmp9.clear();
         Contract(1.0,rop,shape(3,4,5),tmp9bis,shape(2,1,7),0.0,tmp9);

         //contract with left hand side
         DArray<5> tmp5;
         Contract(1.0,LI8,shape(7,5,0,2,3,4),tmp9,shape(7,1,0,3,4,6),0.0,tmp5);

         val -= 2.0 * Dot(tmp5,peps(row+1,col));

         return val;

      }

   }

   /**
    * evaluate the cost function of the linear system for two horizontally connected PEPS: -2 <\Psi|\Psi'> + <\Psi'|\Psi'> where \Psi full PEPS with operator and \Psi is old PEPS
    * for top or bottom rows, i.e. with L and R of order 5 and intermediates of order 7
    * @param row the row index of the bottom site
    * @param col column index of the vertical column
    * @param peps, full PEPS object before update
    * @param L Left environment contraction
    * @param R Right environment contraction
    * @param LI7 left intermediate object
    * @param RI7 right intermediate object
    */
   double cost_function_horizontal(int row,int col,PEPS<double> &peps,const DArray<6> &lop,const DArray<6> &rop,const DArray<5> &L,const DArray<5> &R,const DArray<7> &LI7,const DArray<7> &RI7){

      if(row == 0){

         // --- (1) calculate overlap of approximated state:

         //add peps to right side
         DArray<8> tmp8;
         Contract(1.0,RI7,shape(4,6),peps(row,col+1),shape(1,4),0.0,tmp8);

         //add second peps
         DArray<5> tmp5;
         Contract(1.0,tmp8,shape(3,6,7,4),peps(row,col+1),shape(1,2,3,4),0.0,tmp5);

         //add peps to left side
         DArray<8> tmp8bis;
         Contract(1.0,LI7,shape(6,4),peps(row,col),shape(0,1),0.0,tmp8bis);

         DArray<5> tmp5bis;
         Contract(1.0,tmp8bis,shape(4,3,5,6),peps(row,col),shape(0,1,2,3),0.0,tmp5bis);

         double val = Dot(tmp5bis,tmp5);

         // --- (2) calculate 'b' part of overlap

         //add right operator to tmp8
         DArray<6> tmp6;
         Contract(1.0,tmp8,shape(3,6,7,4),rop,shape(1,2,4,5),0.0,tmp6);

         //attach LI7 to right side
         DArray<7> tmp7;
         Gemm(CblasTrans,CblasNoTrans,1.0,LI7,tmp6,0.0,tmp7);

         //now paste left operator in
         tmp5.clear();
         Contract(1.0,tmp7,shape(2,0,6,5),lop,shape(0,1,3,5),0.0,tmp5);

         //and contract with peps(row,col)
         tmp5bis.clear();
         Permute(tmp5,shape(1,0,3,4,2),tmp5bis);

         val -= 2.0 * Dot(tmp5bis,peps(row,col));

         return val;


      }
      else{//row != 0

      }

      return 0.0;

   }

   /**
    * evaluate the cost function of the linear system for two horizontally connected PEPS: -2 <\Psi|\Psi'> + <\Psi'|\Psi'> where \Psi full PEPS with operator and \Psi is old PEPS
    * for middle rows, with R and L of order 6 and intermediates of order 8
    * @param row the row index of the bottom site
    * @param col column index of the vertical column
    * @param peps, full PEPS object before update
    * @param LO Left environment contraction
    * @param RO Right environment contraction
    * @param LI8 left interediate object
    * @param RI8 right interediate object
    */
   double cost_function_horizontal(int row,int col,PEPS<double> &peps,const DArray<6> &lop,const DArray<6> &rop,const DArray<6> &LO,const DArray<6> &RO,const DArray<8> &LI8,const DArray<8> &RI8){

      //(1) overlap term

      //add right peps to intermediate
      DArray<9> tmp9;
      Contract(1.0,RI8,shape(4,6),peps(row,col+1),shape(1,4),0.0,tmp9);

      //add second peps
      DArray<8> tmp8;
      Contract(1.0,tmp9,shape(3,7,4),peps(row,col+1),shape(1,2,4),0.0,tmp8);

      //add bottom environment
      DArray<6> tmp6;
      Contract(1.0,tmp8,shape(7,5,3),env.gb(row-1)[col+1],shape(1,2,3),0.0,tmp6);

      //add left peps to LI8
      DArray<9> tmp9bis;
      Contract(1.0,LI8,shape(6,4),peps(row,col),shape(0,1),0.0,tmp9bis);

      //and another peps
      tmp8.clear();
      Contract(1.0,tmp9bis,shape(4,3,6),peps(row,col),shape(0,1,2),0.0,tmp8);

      //now contract with bottom environment
      DArray<6> tmp6bis;
      Contract(1.0,tmp8,shape(3,6,4),env.gb(row-1)[col],shape(0,1,2),0.0,tmp6bis);

      double val = Dot(tmp6,tmp6bis);

      // (2) construct right hand side

      //add right operator
      DArray<9> tmp9op;
      Contract(1.0,tmp9,shape(3,7,4),rop,shape(1,2,5),0.0,tmp9op);

      //add bottom environment
      DArray<7> tmp7;
      Contract(1.0,tmp9op,shape(8,5,3),env.gb(row-1)[col+1],shape(1,2,3),0.0,tmp7);

      tmp9op.clear();
      Contract(1.0,tmp9bis,shape(4,3,6),lop,shape(0,1,2),0.0,tmp9op);

      //now contract with bottom environment
      DArray<7> tmp7bis;
      Contract(1.0,tmp9op,shape(3,7,4),env.gb(row-1)[col],shape(0,1,2),0.0,tmp7bis);

      DArray<7> perm7;
      Permute(tmp7bis,shape(0,1,2,3,5,4,6),perm7) ;

      val -= 2.0 * Dot(perm7,tmp7);

      return val;

   }

}
