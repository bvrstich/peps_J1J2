#include <iostream>
#include <iomanip>
#include <fstream>
#include <complex>
#include <chrono>
#include <omp.h>

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
    * update the tensors in a sweeping fashion, for bottom or top rows, i.e. with R and L environments of order 5
    * @param dir vertical, horizontal,diagonal lurd or diagonal ldru update
    * @param row , the row index of the bottom site
    * @param col column index of the vertical column
    * @param peps, full PEPS object before update
    * @param L Left environment contraction
    * @param R Right environment contraction
    * @param n_iter nr of sweeps in the ALS algorithm
    */
   template<size_t M>
      void update(const PROP_DIR &dir,int row,int col,PEPS<double> &peps,const DArray<M> &L,const DArray<M> &R,int n_iter){

         enum {i,j,k,l,m,n,o};

         //containers for left and right intermediary objects
         DArray<M+2> LI;
         DArray<M+2> RI;

         DArray<M+2> b_L;
         DArray<M+2> b_R;

         //construct for left and right operators acting the correct peps elements
         DArray<6> lop;
         DArray<6> rop;

         //for diagonal update we need a 'middle' site connecting left and right
         DArray<5> mop;

         if(dir == VERTICAL){// (row,col) --> (row+1,col)

            Contract(1.0,peps(row,col),shape(i,j,k,l,m),global::trot.gLO_n(),shape(k,o,n),0.0,lop,shape(i,j,n,o,l,m));
            Contract(1.0,peps(row+1,col),shape(i,j,k,l,m),global::trot.gRO_n(),shape(k,o,n),0.0,rop,shape(i,j,n,o,l,m));

         }
         else if(dir == HORIZONTAL){// (row,col) --> (row,col+1)

            Contract(1.0,peps(row,col),shape(i,j,k,l,m),global::trot.gLO_n(),shape(k,o,n),0.0,lop,shape(i,j,n,o,l,m));
            Contract(1.0,peps(row,col+1),shape(i,j,k,l,m),global::trot.gRO_n(),shape(k,o,n),0.0,rop,shape(i,j,n,o,l,m));

         }
         else if(dir == DIAGONAL_LURD){//(row+1,col) --> (row,col+1)

            //middle peps is left bottom 
            mop = peps(row,col);

            Contract(1.0,peps(row+1,col),shape(i,j,k,l,m),global::trot.gLO_nn(),shape(k,o,n),0.0,lop,shape(i,j,n,o,l,m));
            Contract(1.0,peps(row,col+1),shape(i,j,k,l,m),global::trot.gRO_nn(),shape(k,o,n),0.0,rop,shape(i,j,n,o,l,m));

         }
         else{//(row,col) --> (row+1,col+1)

            //middle peps is bottom right
            mop = peps(row,col+1);

            Contract(1.0,peps(row,col),shape(i,j,k,l,m),global::trot.gLO_nn(),shape(k,o,n),0.0,lop,shape(i,j,n,o,l,m));
            Contract(1.0,peps(row+1,col+1),shape(i,j,k,l,m),global::trot.gRO_nn(),shape(k,o,n),0.0,rop,shape(i,j,n,o,l,m));
         }

         // --- (a) --- initial guess:
         initialize(dir,row,col,lop,rop,peps); 

         // --- (b) --- create intermediary object using during N_eff construction, doesn't change during update
         construct_intermediate(dir,row,col,peps,mop,L,R,LI,RI,b_L,b_R);
       
         // --- (c) --- sweeping update
         sweep(dir,row,col,peps,lop,rop,L,R,LI,RI,b_L,b_R,n_iter);

         // --- (d) --- set top and bottom back on equal footing
         //equilibrate(dir,row,col,peps);

      }

   /**
    * Sweep back and forward between the two peps to be updated, solving the linear system for compression until convergence is reached
    * @param dir vertical, horizontal,diagonal lurd or diagonal ldru update
    * @param row row index of the bottom left peps
    * @param col col index of the bottom left peps
    * @param lop bottom peps acted on with left trotter operator
    * @param rop upper peps acted on with right trotter operator
    * @param L left contracted environment 
    * @param R right contracted environment 
    * @param LI intermediate object created to simplify N_eff construction
    * @param RI intermediate object created to simplify N_eff construction
    * @param n_sweeps number of sweeps to execute
    */
   template<size_t M>
      void sweep(const PROP_DIR &dir,int row,int col,PEPS<double> &peps,const DArray<6> &lop,const DArray<6> &rop,

            const DArray<M> &L,const DArray<M> &R,const DArray<M+2> &LI,const DArray<M+2> &RI,
            
            const DArray<M+2> &b_L,const DArray<M+2> &b_R, int n_sweeps){

         //indices of sites between which to jump back and forth
         int l_row(row),l_col(col),r_row(row),r_col(col);

         if(dir == VERTICAL)//(row,col) --> (row+1,col)
            r_row++;
         else if(dir == HORIZONTAL)//(row,col) --> (row,col+1)
            r_col++;
         else if(dir == DIAGONAL_LURD){//(row+1,col) --> (row,col+1)

            l_row++;
            r_col++;

         }
         else{//(row,col) --> (row+1,col+1)

            r_row++;
            r_col++;

         }

         //storage
         DArray<8> N_eff;
         DArray<5> rhs;

         int iter = 0;

         //while(iter < n_sweeps){

            // cout << iter << "\t" << debug::cost_function(dir,row,col,peps,lop,rop,L,R,LI,RI,b_L,b_R) << endl;

            // --(1)-- 'left' site

            //construct effective environment and right hand side for linear system of top site
            construct_lin_sys(dir,row,col,peps,lop,rop,N_eff,rhs,L,R,LI,RI,b_L,b_R,true);

            DArray<1> eig;
            debug::diagonalize(N_eff,eig);

            cout << endl;
            cout << "left\t" << eig(0) << "\t" << eig(eig.size() - 1) << "\t|\t" << eig(eig.size()-1)/eig(0) << endl;
/*
            //solve the system
            solve(N_eff,rhs);

            //update 'left' peps
            Permute(rhs,shape(0,1,4,2,3),peps(l_row,l_col));
    */
            // --(2)-- 'right' site

            //construct effective environment and right hand side for linear system of bottom site
            construct_lin_sys(dir,row,col,peps,lop,rop,N_eff,rhs,L,R,LI,RI,b_L,b_R,false);

            eig.clear();
            debug::diagonalize(N_eff,eig);

            cout << "right\t" << eig(0) << "\t" << eig(eig.size() - 1) << "\t|\t" << eig(eig.size()-1)/eig(0) << endl;
            cout << endl;
/*
            //solve the system
            solve(N_eff,rhs);

            //update 'right' peps
            Permute(rhs,shape(0,1,4,2,3),peps(r_row,r_col));
  */
            //repeat until converged
            ++iter;
  
         //}

      }

   /**
    * propagate the peps one imaginary time step
    * @param peps the PEPS to be propagated
    * @param n_sweeps the number of sweeps performed for the solution of the linear problem
    */
   void step(PEPS<double> &peps,int n_sweeps){

      enum {i,j,k,l,m,n,o};

      // ****************************************** //
      // ****************************************** //
      //    (A)     FIRST THE EVEN ROWS             //
      // ****************************************** //
      // ****************************************** //

      //calculate top and bottom environment
      env.calc('A',peps);

#pragma omp parallel for schedule(dynamic,1)
      for(int row = 0;row < Ly;row+=2){

         if(row == 0){

            //containers for the renormalized operators
            vector< DArray<5> > R(Lx - 1);
            DArray<5> L;

            //initialize the right operators for the bottom row
            contractions::init_ro('b',peps,R);

            //for(int col = 0;col < Lx - 1;++col){
           int col = 0; 

               // --- (1) update the vertical pair on column 'col' ---
               update(VERTICAL,0,col,peps,L,R[col],n_sweeps); 

               // --- (2) update the horizontal pair on column 'col'-'col+1' ---
               //update(HORIZONTAL,0,col,peps,L,R[col+1],n_sweeps); 

               // --- (3) update diagonal LU-RD
               //update(DIAGONAL_LURD,0,col,peps,L,R[col+1],n_sweeps); 

               // --- (4) update diagonal LD-RU
               //update(DIAGONAL_LDRU,0,col,peps,L,R[col+1],n_sweeps); 

               contractions::update_L('b',col,peps,L);

            //}

            //one last vertical update
//            update(VERTICAL,0,Lx-1,peps,L,R[Lx-2],n_sweeps); 

         }
         else if(row < Lx - 2){

            //renormalized operators for the middle sites
            vector< DArray<6> > RO(Lx - 1);
            DArray<6> LO;

            //first create right renormalized operator
            contractions::init_ro(row,peps,RO);

            for(int col = 0;col < Lx - 1;++col){

               // --- (1) update vertical pair on column 'col', with lowest site on row 'row'
               update(VERTICAL,row,col,peps,LO,RO[col],n_sweeps); 

               // --- (2) update the horizontal pair on column 'col'-'col+1' ---
               update(HORIZONTAL,row,col,peps,LO,RO[col+1],n_sweeps); 

               // --- (3) update diagonal LU-RD
               //update(DIAGONAL_LURD,row,col,peps,LO,RO[col+1],n_sweeps); 

               // --- (4) update diagonal LD-RU
               //update(DIAGONAL_LDRU,row,col,peps,LO,RO[col+1],n_sweeps); 

               //first construct a double layer object for the newly updated bottom 
               contractions::update_L(row,col,peps,LO);

            }

            //finally, last vertical gate
            update(VERTICAL,row,Lx-1,peps,LO,RO[Lx-2],n_sweeps); 

         }
         else{//row == Lx - 2

            //containers for the renormalized operators
            vector< DArray<5> > R(Lx - 1);
            DArray<5> L;

            //make the right operators
            contractions::init_ro('t',peps,R);

            for(int col = 0;col < Lx - 1;++col){

               // --- (1) update vertical pair on column 'col' on upper two rows
               //update(VERTICAL,Ly-2,col,peps,L,R[col],n_sweeps); 

               // --- (2a) update the horizontal pair on row 'row' and colums 'col'-'col+1' ---
               //update(HORIZONTAL,Ly-2,col,peps,L,R[col+1],n_sweeps); 

               // --- (2b) update the horizontal pair on row 'row+1' and colums 'col'-'col+1' ---
               //update(HORIZONTAL,Ly-1,col,peps,L,R[col+1],n_sweeps); 

               // --- (3) update diagonal LU-RD
               //update(DIAGONAL_LURD,Ly-2,col,peps,L,R[col+1],n_sweeps); 

               // --- (4) update diagonal LD-RU
               //update(DIAGONAL_LDRU,Ly-2,col,peps,L,R[col+1],n_sweeps); 

               contractions::update_L('t',col,peps,L);

            }

            //finally the very last vertical gate
            //update(VERTICAL,Ly-2,Lx-1,peps,L,R[Lx-2],n_sweeps); 

         }

      }

      // ****************************************** //
      // ****************************************** //
      //    (B)     THEN THE ODD ROWS               //
      // ****************************************** //
      // ****************************************** //
/*
      //update top and bottom environment
      global::env.calc('A',peps);

#pragma omp parallel for schedule(static,1)
      for(int row = 1;row < Ly-2;row+=2){

         //renormalized operators for the middle sites
         vector< DArray<6> > RO(Lx - 1);
         DArray<6> LO;

         //first create right renormalized operator
         contractions::init_ro(row,peps,RO);

         for(int col = 0;col < Lx - 1;++col){

            // --- (1) update vertical pair on column 'col', with lowest site on row 'row'
            update(VERTICAL,row,col,peps,LO,RO[col],n_sweeps); 

            // --- (2) update the horizontal pair on column 'col'-'col+1' ---
            update(HORIZONTAL,row,col,peps,LO,RO[col+1],n_sweeps); 

            // --- (3) update diagonal LU-RD
            //update(DIAGONAL_LURD,row,col,peps,LO,RO[col+1],n_sweeps); 

            // --- (4) update diagonal LD-RU
            //update(DIAGONAL_LDRU,row,col,peps,LO,RO[col+1],n_sweeps); 

            //first construct a double layer object for the newly updated bottom 
            contractions::update_L(row,col,peps,LO);

         }

         //finally, last vertical gate
         update(VERTICAL,row,Lx-1,peps,LO,RO[Lx-2],n_sweeps); 

      }
*/
   }

   /** 
    * wrapper function solve symmetric linear system: N_eff * x = b, symmetrize the N_eff first
    * @param N_eff input matrix
    * @param rhs right hand side input and x output
    */
   void solve(DArray<8> &N_eff,DArray<5> &rhs){

      int matdim = N_eff.shape(0) * N_eff.shape(1) * N_eff.shape(2) * N_eff.shape(3);

      //symmetrize
      for(int i = 0;i < matdim;++i)
         for(int j = i + 1;j < matdim;++j){

            N_eff.data()[i*matdim + j] = 0.5 * (N_eff.data()[i*matdim + j]  + N_eff.data()[j*matdim + i]);
            N_eff.data()[j*matdim + i] = N_eff.data()[i*matdim + j];

         }

      int *ipiv = new int [matdim];

      lapack::sytrf(CblasRowMajor,'U',matdim, N_eff.data(), matdim,ipiv);

      lapack::sytrs(CblasRowMajor,'U',matdim,d, N_eff.data(),matdim,ipiv, rhs.data(),d);

      delete [] ipiv;

   }

   /**
    * construct the single-site effective environment and right hand side needed for 
    * the linear system of any gate direction specified by 'dir'
    * @param dir vertical, horizontal,diagonal lurd or diagonal ldru update
    * @param row the row index of the bottom site
    * @param col column index of the vertical column
    * @param peps, full PEPS object before update
    * @param N_eff output object, contains N_eff on output
    * @param rhs output object, contains N_eff on output
    * @param L Left environment contraction
    * @param R Right environment contraction
    * @param LI7 left intermediate object
    * @param RI7 right intermediate object
    * @param b_L left intermediate object
    * @param b_R right intermediate object
    * @param left boolean flag for peps with left operator or right operator acting on it
    */
   template<>
      void construct_lin_sys(const PROP_DIR &dir,int row,int col,PEPS<double> &peps,const DArray<6> &lop,const DArray<6> &rop,
            
            DArray<8> &N_eff,DArray<5> &rhs, const DArray<5> &L, const DArray<5> &R, const DArray<7> &LI7,const DArray<7> &RI7,
            
            const DArray<7> &b_L,const DArray<7> &b_R,bool left){

         if(dir == VERTICAL){

            if(row == 0){

               if(col == 0){

                  if(left){//bottom site of vertical, so site (row,col) environment

                     calc_N_eff();

                     calc_rhs();
                  }
                  else{//top site (row,col)

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

               }
               else if(col < Lx - 1){//col != 0

                  if(left){//bottom site

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
                  else{//bottom site

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

               }
               else{ //col == Lx - 1

                  if(left){//top site

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
                  else{//bottom site

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

               }

            }
            else{//row = Lx - 2

               if(col == 0){

                  if(left){//top site of vertical, so site (row,col+1) environment

                     //(1) construct N_eff

                     //paste top peps to right intermediate
                     DArray<10> tmp10;
                     Contract(1.0,peps(row+1,col),shape(4),RI7,shape(4),0.0,tmp10);

                     DArray<7> tmp7;
                     Contract(1.0,peps(row+1,col),shape(0,1,2,4),tmp10,shape(0,1,2,7),0.0,tmp7);

                     DArray<7> tmp7bis;
                     Permute(tmp7,shape(2,0,3,5,1,4,6),tmp7bis);

                     N_eff = tmp7bis.reshape_clear(shape(1,D,D,D,1,D,D,D));

                     //(2) right hand side

                     //paste right operator to tmp10
                     DArray<8> tmp8;
                     Contract(1.0,rop,shape(0,1,2,5),tmp10,shape(0,1,2,7),0.0,tmp8);

                     //add left operator to tmp8
                     DArray<4> tmp4;
                     Contract(1.0,lop,shape(0,1,3,4,5),tmp8,shape(3,1,0,4,6),0.0,tmp4);

                     DArray<4> tmp4bis;
                     Permute(tmp4,shape(1,2,3,0),tmp4bis);

                     rhs = tmp4bis.reshape_clear( shape(1,D,D,D,d) );
                  }
                  else{//bottom site (row,col)

                     // (1) construct N_eff

                     //paste bottom peps to right intermediate
                     DArray<8> tmp8;
                     Contract(1.0,peps(row,col),shape(3,4),RI7,shape(2,6),0.0,tmp8);

                     //and another bottom peps to tmp8
                     DArray<7> tmp7;
                     Contract(1.0,peps(row,col),shape(2,3,4),tmp8,shape(2,4,7),0.0,tmp7);

                     DArray<7> tmp7bis;
                     Permute(tmp7,shape(0,4,1,5,2,3,6),tmp7bis);

                     N_eff = tmp7bis.reshape_clear(shape(1,D,1,D,1,D,1,D));

                     //(2) right hand side:

                     //add left operator to tmp8
                     DArray<6> tmp6;
                     Contract(1.0,lop,shape(0,2,4,5),tmp8,shape(0,2,4,7),0.0,tmp6);

                     //add right operator
                     DArray<6> tmp6bis;
                     Contract(1.0,rop,shape(3,4,5),tmp6,shape(1,0,4),0.0,tmp6bis);

                     tmp6.clear();
                     Permute(tmp6bis,shape(3,5,2,0,1,4),tmp6);

                     rhs = tmp6.reshape_clear( shape(1,1,D,D,d) );

                  }

               }
               else if(col < Lx - 1){//middle columns

                  if(left){//top site of vertical, so site (row,col+1) environment

                     //(1) construct N_eff

                     //first add top peps to RI7
                     DArray<10> tmp10;
                     Contract(1.0,peps(row+1,col),shape(4),RI7,shape(4),0.0,tmp10);

                     //then add another top peps
                     DArray<9> tmp9;
                     Contract(1.0,peps(row+1,col),shape(1,2,4),tmp10,shape(1,2,7),0.0,tmp9);

                     //now contract with left side
                     DArray<8> tmp8;
                     Contract(1.0,L,shape(0,1,4),tmp9,shape(0,2,4),0.0,tmp8);

                     N_eff.clear();
                     Permute(tmp8,shape(0,2,4,6,1,3,5,7),N_eff);

                     //(2) right hand side

                     //add right operator to tmp10
                     DArray<10> tmp10bis;
                     Contract(1.0,rop,shape(1,2,5),tmp10,shape(1,2,7),0.0,tmp10bis);

                     //and left operator
                     tmp8.clear();
                     Contract(1.0,tmp10bis,shape(2,1,6,8),lop,shape(1,3,4,5),0.0,tmp8);

                     //now contract with left side
                     rhs.clear();
                     Contract(1.0,L,shape(0,1,2,4),tmp8,shape(0,1,6,3),0.0,rhs);

                  }
                  else{//bottom site (row,col)

                     // (1) construct N_eff

                     //paste bottom peps to right intermediate
                     DArray<8> tmp8;
                     Contract(1.0,peps(row,col),shape(3,4),RI7,shape(2,6),0.0,tmp8);

                     //and another bottom peps to tmp8
                     DArray<7> tmp7;
                     Contract(1.0,peps(row,col),shape(2,3,4),tmp8,shape(2,4,7),0.0,tmp7);

                     //contract with left side
                     DArray<6> tmp6;
                     Contract(1.0,L,shape(2,3,4),tmp7,shape(0,2,4),0.0,tmp6);

                     DArray<6> tmp6bis;
                     Permute(tmp6,shape(0,2,4,1,3,5),tmp6bis);

                     N_eff = tmp6bis.reshape_clear(shape(D,1,D,D,D,1,D,D));

                     //(2) right hand side:

                     //add left operator to tmp8
                     DArray<8> tmp8bis;
                     Contract(1.0,lop,shape(2,4,5),tmp8,shape(2,4,7),0.0,tmp8bis);

                     //add right operator
                     tmp8.clear();
                     Contract(1.0,rop,shape(3,4,5),tmp8bis,shape(2,1,6),0.0,tmp8);

                     //contract with left side
                     DArray<5> tmp5;
                     Contract(1.0,L,shape(0,2,3,4),tmp8,shape(0,3,4,6),0.0,tmp5);

                     Permute(tmp5,shape(0,1,3,4,2),rhs);

                  }

               }
               else{//col == Lx - 1

                  if(left){//top site of vertical, so site (row,col+1) environment

                     //(1) construct N_eff

                     //first add top peps to LI7
                     DArray<8> tmp8;
                     Contract(1.0,LI7,shape(1,6),peps(row+1,col),shape(0,4),0.0,tmp8);

                     //then add another top peps
                     DArray<7> tmp7;
                     Contract(1.0,tmp8,shape(0,5,6),peps(row+1,col),shape(0,1,2),0.0,tmp7);

                     DArray<7> tmp7bis;
                     Permute(tmp7,shape(0,5,2,6,1,4,3),tmp7bis);

                     N_eff = tmp7bis.reshape_clear( shape(D,D,D,1,D,D,D,1) );

                     //(2) right hand side

                     //add right operator to tmp8
                     DArray<8> tmp8bis;
                     Contract(1.0,tmp8,shape(0,5,6),rop,shape(0,1,2),0.0,tmp8bis);

                     //and left operator
                     DArray<4> tmp4;
                     Contract(1.0,tmp8bis,shape(0,6,5,2,7),lop,shape(0,1,3,4,5),0.0,tmp4);

                     DArray<4> tmp4bis;
                     Permute(tmp4,shape(0,2,1,3),tmp4bis);

                     rhs = tmp4bis.reshape_clear( shape(D,D,D,1,d) );

                  }
                  else{//bottom site (row,col)

                     // (1) construct N_eff

                     //paste bottom peps to left intermediate
                     DArray<6> tmp6;
                     Contract(1.0,peps(row,col),shape(0,3,4),LI7,shape(3,5,6),0.0,tmp6);

                     //and another
                     DArray<5> tmp5;
                     Contract(1.0,peps(row,col),shape(0,2,3),tmp6,shape(4,1,5),0.0,tmp5);

                     DArray<5> tmp5bis;
                     Permute(tmp5,shape(3,0,4,2,1),tmp5bis);

                     N_eff = tmp5bis.reshape_clear( shape(D,1,D,1,D,1,D,1) );

                     //(2) right hand side:

                     //add left operator to tmp6
                     DArray<6> tmp6bis;
                     Contract(1.0,lop,shape(0,2,4),tmp6,shape(4,1,5),0.0,tmp6bis);

                     //add right operator
                     DArray<4> tmp4;
                     Contract(1.0,rop,shape(0,3,4,5),tmp6bis,shape(4,1,0,2),0.0,tmp4);

                     DArray<4> tmp4bis;
                     Permute(tmp4,shape(3,0,2,1),tmp4bis);

                     rhs = tmp4bis.reshape_clear( shape(D,1,D,1,d) );

                  }

               }

            }

         }//close VERTICAL
         else if(dir == HORIZONTAL){

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

                  //(1) construct N_eff

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
            else if(row == Ly - 2){//bottom horizontal peps of topmost update

               if(left){//left site of horizontal gate, so site (row,col) environment

                  //(1) construct N_eff

                  //add right peps to intermediate
                  DArray<8> tmp8;
                  Contract(1.0,RI7,shape(3,5),peps(row,col+1),shape(1,4),0.0,tmp8);

                  //add second peps
                  DArray<7> tmp7;
                  Contract(1.0,tmp8,shape(2,6,3),peps(row,col+1),shape(1,2,4),0.0,tmp7);

                  //add bottom environment
                  DArray<5> tmp5;
                  Contract(1.0,env.gb(Ly-3)[col+1],shape(1,2,3),tmp7,shape(6,4,2),0.0,tmp5);

                  //add next bottom
                  tmp7.clear();
                  Contract(1.0,env.gb(Ly-3)[col],shape(3),tmp5,shape(0),0.0,tmp7);

                  //contract with left hand side
                  DArray<8> tmp8bis;
                  Contract(1.0,LI7,shape(0,1,6),tmp7,shape(3,4,0),0.0,tmp8bis);

                  N_eff.clear();
                  Permute(tmp8bis,shape(2,0,4,7,3,1,5,6),N_eff);

                  // (2) construct right hand side

                  //add right operator to tmp8
                  tmp8bis.clear();
                  Contract(1.0,tmp8,shape(2,6,3),rop,shape(1,2,5),0.0,tmp8bis);

                  //add bottom environment
                  DArray<6> tmp6;
                  Contract(1.0,env.gb(Ly-3)[col+1],shape(1,2,3),tmp8bis,shape(7,4,2),0.0,tmp6);

                  //add next bottom
                  tmp8.clear();
                  Contract(1.0,env.gb(Ly-3)[col],shape(3),tmp6,shape(0),0.0,tmp8);

                  //now paste left operator in
                  tmp8bis.clear();
                  Contract(1.0,lop,shape(3,4,5),tmp8,shape(7,1,6),0.0,tmp8bis);

                  //contract with left hand side
                  tmp5.clear();
                  Contract(1.0,LI7,shape(0,1,2,4,6),tmp8bis,shape(5,6,1,0,3),0.0,tmp5);

                  Permute(tmp5,shape(1,0,3,4,2),rhs);

               }
               else{//right site of horizontal gate, so site (row+1,col) environment

                  //(1) constsruct N_eff

                  //add left peps to LI7
                  DArray<8> tmp8;
                  Contract(1.0,LI7,shape(5,3),peps(row,col),shape(0,1),0.0,tmp8);

                  //and another peps
                  DArray<7> tmp7;
                  Contract(1.0,tmp8,shape(3,2,5),peps(row,col),shape(0,1,2),0.0,tmp7);

                  //add bottom environment
                  DArray<5> tmp5;
                  Contract(1.0,tmp7,shape(2,5,3),env.gb(Ly-3)[col],shape(0,1,2),0.0,tmp5);

                  //add next bottom environment
                  tmp7.clear();
                  Contract(1.0,tmp5,shape(4),env.gb(Ly-3)[col+1],shape(0),0.0,tmp7);

                  //now contract with right
                  DArray<8> tmp8bis;
                  Contract(1.0,tmp7,shape(0,1,6),RI7,shape(0,1,6),0.0,tmp8bis);

                  N_eff.clear();
                  Permute(tmp8bis,shape(1,4,2,6,0,5,3,7),N_eff);

                  // (2) construct right hand side

                  //add left operator
                  tmp8bis.clear();
                  Contract(1.0,tmp8,shape(3,2,5),lop,shape(0,1,2),0.0,tmp8bis);

                  //add bottom environment
                  DArray<6> tmp6;
                  Contract(1.0,tmp8bis,shape(2,6,3),env.gb(Ly-3)[col],shape(0,1,2),0.0,tmp6);

                  //add next bottom environment
                  tmp8.clear();
                  Contract(1.0,tmp6,shape(5),env.gb(Ly-3)[col+1],shape(0),0.0,tmp8);

                  //now contract with right operator
                  tmp8bis.clear();
                  Contract(1.0,tmp8,shape(4,3,5),rop,shape(0,3,4),0.0,tmp8bis);

                  tmp5.clear();
                  Contract(1.0,tmp8bis,shape(0,1,5,7,4),RI7,shape(0,1,2,4,6),0.0,tmp5);

                  rhs.clear();
                  Permute(tmp5,shape(0,3,1,4,2),rhs);

               }

            }
            else{//row == Ly - 1

               if(left){//left site of horizontal gate, so site (row,col) environment

                  //(1) construct N_eff

                  //add right peps to intermediate
                  DArray<8> tmp8;
                  Contract(1.0,peps(row,col+1),shape(3,4),RI7,shape(3,1),0.0,tmp8);

                  //add second peps
                  DArray<5> tmp5;
                  Contract(1.0,peps(row,col+1),shape(1,2,3,4),tmp8,shape(1,2,4,3),0.0,tmp5);

                  //contract with left hand side
                  DArray<6> tmp6;
                  Gemm(CblasNoTrans,CblasTrans,1.0,LI7,tmp5,0.0,tmp6);

                  DArray<6> tmp6bis;
                  Permute(tmp6,shape(0,2,4,1,3,5),tmp6bis);

                  int DL = peps(row,col).shape(0);

                  N_eff = tmp6bis.reshape_clear( shape(DL,1,D,D,DL,1,D,D) );

                  // (2) construct right hand side

                  //add right operator to tmp8
                  tmp6.clear();
                  Contract(1.0,rop,shape(1,2,4,5),tmp8,shape(1,2,4,3),0.0,tmp6);

                  //now paste left operator in
                  tmp8.clear();
                  Contract(1.0,lop,shape(5,3),tmp6,shape(0,1),0.0,tmp8);

                  //contract with left hand side
                  tmp5.clear();
                  Contract(1.0,LI7,shape(0,2,4,5,6),tmp8,shape(0,3,5,6,7),0.0,tmp5);

                  rhs.clear();
                  Permute(tmp5,shape(0,2,1,4,3),rhs);

               }
               else{//right site of horizontal gate, so site (row+1,col) environment

                  //(1) construct N_eff

                  //add left to intermediate
                  DArray<8> tmp8;
                  Contract(1.0,LI7,shape(1,3),peps(row,col),shape(0,3),0.0,tmp8);

                  //and another
                  DArray<5> tmp5;
                  Contract(1.0,tmp8,shape(0,5,6,1),peps(row,col),shape(0,1,2,3),0.0,tmp5);

                  //contract with right side
                  DArray<6> tmp6;
                  Gemm(CblasNoTrans,CblasNoTrans,1.0,RI7,tmp5,0.0,tmp6);

                  DArray<6> tmp6bis;
                  Permute(tmp6,shape(5,2,0,4,3,1),tmp6bis);

                  int DR = peps(row,col+1).shape(4);

                  N_eff = tmp6bis.reshape_clear( shape(D,1,D,DR,D,1,D,DR) );

                  // (2) construct right hand side

                  //add left operator
                  tmp6.clear();
                  Contract(1.0,tmp8,shape(0,5,6,1),lop,shape(0,1,2,4),0.0,tmp6);

                  //and right
                  tmp8.clear();
                  Contract(1.0,tmp6,shape(5,4),rop,shape(0,3),0.0,tmp8);

                  //contract with RI7 hand side
                  tmp5.clear();
                  Contract(1.0,tmp8,shape(7,6,0,1,2),RI7,shape(0,2,4,5,6),0.0,tmp5);

                  rhs.clear();
                  Permute(tmp5,shape(0,1,4,3,2),rhs);

               }

            }

         }//close HORIZONTAL
         else if(dir == DIAGONAL_LURD){

            if(row == 0){

               if(left){//left-up site of diagonal gate, so site (row+1,col) environment

                  // (1) construct N_eff
                  DArray<8> tmp8;
                  Contract(1.0,RI7,shape(4,6),peps(row,col+1),shape(1,4),0.0,tmp8);

                  //und again
                  DArray<5> tmp5;
                  Contract(1.0,tmp8,shape(3,6,7,4),peps(row,col+1),shape(1,2,3,4),0.0,tmp5);

                  //add top environment to intermediate
                  DArray<7> tmp7;
                  Contract(1.0,env.gt(row)[col],shape(3),tmp5,shape(0),0.0,tmp7);

                  //add LI7 to it
                  DArray<8> tmp8bis;
                  Contract(1.0,LI7,shape(0,5,6),tmp7,shape(0,6,5),0.0,tmp8bis);

                  N_eff.clear();
                  Permute(tmp8bis,shape(0,4,2,6,1,5,3,7),N_eff);

                  // (2) construct right hand side

                  //add right operator to tmp8
                  DArray<6> tmp6;
                  Contract(1.0,tmp8,shape(3,6,7,4),rop,shape(1,2,4,5),0.0,tmp6);

                  //add top environment to intermediate
                  tmp8.clear();
                  Contract(1.0,env.gt(row)[col],shape(3),tmp6,shape(0),0.0,tmp8);

                  //add left operator to b_L
                  DArray<9> tmp9;
                  Contract(1.0,b_L,shape(1,3),lop,shape(0,4),0.0,tmp9);

                  //contract both sides to form right hand side of equation
                  tmp5.clear();
                  Contract(1.0,tmp9,shape(0,5,8,7,3,4),tmp8,shape(0,1,3,7,6,5),0.0,tmp5);

                  rhs.clear();
                  Permute(tmp5,shape(0,3,1,4,2),rhs);

               }
               else{//right-down site of diagonal gate, so site (row,col + 1) environment

                  // (1) construct N_eff

                  //add upper-left peps to intermediate left
                  DArray<8> tmp8;
                  Contract(1.0,LI7,shape(2,4),peps(row+1,col),shape(0,3),0.0,tmp8);

                  //und again
                  DArray<7> tmp7;
                  Contract(1.0,tmp8,shape(1,6,2),peps(row+1,col),shape(0,2,3),0.0,tmp7);

                  //add top environment to intermediate
                  DArray<5> tmp5;
                  Contract(1.0,tmp7,shape(0,5,3),env.gt(row)[col],shape(0,1,2),0.0,tmp5);

                  //add left and right together
                  DArray<6> tmp6;
                  Contract(1.0,tmp5,shape(4,3,2),RI7,shape(0,1,2),0.0,tmp6);

                  DArray<6> tmp6bis;
                  Permute(tmp6,shape(0,2,4,1,3,5),tmp6bis);

                  int DL = peps(row,col+1).shape(0);
                  int DU = peps(row,col+1).shape(1);
                  int DD = peps(row,col+1).shape(3);
                  int DR = peps(row,col+1).shape(4);

                  N_eff = tmp6bis.reshape_clear( shape(DL,DU,DD,DR,DL,DU,DD,DR) );

                  // (2) construct right hand side

                  //add upper-left peps to intermediate b_L
                  Contract(1.0,b_L,shape(2,4),peps(row+1,col),shape(0,3),0.0,tmp8);

                  //add left operator to tmp8
                  DArray<8> tmp8bis;
                  Contract(1.0,tmp8,shape(1,6,2),lop,shape(0,2,4),0.0,tmp8bis);

                  //add top environment to intermediate
                  tmp6.clear();
                  Contract(1.0,tmp8bis,shape(0,5,3),env.gt(row)[col],shape(0,1,2),0.0,tmp6);

                  //add right operator to RI7
                  DArray<9> tmp9;
                  Contract(1.0,RI7,shape(3,5),rop,shape(1,5),0.0,tmp9);

                  //now contract left and right
                  tmp5.clear();
                  Contract(1.0,tmp6,shape(5,4,3,2,0),tmp9,shape(0,1,7,2,5),0.0,tmp5);

                  rhs.clear();
                  Permute(tmp5,shape(0,1,4,2,3),rhs);

               }

            }
            else{//top two rows

               if(left){//left-up site of diagonal gate, so site (row+1,col) environment

                  // (1) construct N_eff

                  //add right peps to intermediate
                  DArray<8> tmp8;
                  Contract(1.0,RI7,shape(3,5),peps(row,col+1),shape(1,4),0.0,tmp8);

                  //add second peps
                  DArray<7> tmp7;
                  Contract(1.0,tmp8,shape(2,6,3),peps(row,col+1),shape(1,2,4),0.0,tmp7);

                  //add bottom environment
                  DArray<5> tmp5;
                  Contract(1.0,env.gb(Ly-3)[col+1],shape(1,2,3),tmp7,shape(6,4,2),0.0,tmp5);

                  //now contract left and right
                  DArray<6> tmp6;
                  Contract(1.0,LI7,shape(4,5,6),tmp5,shape(4,3,0),0.0,tmp6);

                  DArray<6> tmp6bis;
                  Permute(tmp6,shape(0,2,4,1,3,5),tmp6bis);

                  int DL = peps(row+1,col).shape(0);
                  int DR = peps(row+1,col).shape(4);

                  N_eff = tmp6bis.reshape_clear( shape(DL,1,D,DR,DL,1,D,DR) );

                  // (2) construct right hand side

                  //add right operator to tmp8
                  DArray<8> tmp8bis;
                  Contract(1.0,tmp8,shape(2,6,3),rop,shape(1,2,5),0.0,tmp8bis);

                  //add bottom environment
                  tmp6.clear();
                  Contract(1.0,env.gb(row-1)[col+1],shape(1,2,3),tmp8bis,shape(7,4,2),0.0,tmp6);

                  //now add left operator
                  tmp8.clear();
                  Contract(1.0,lop,shape(3,5),tmp6,shape(5,1),0.0,tmp8);

                  //contract with b_L
                  tmp5.clear(); 
                  Contract(1.0,b_L,shape(0,2,4,5,6),tmp8,shape(0,3,7,6,4),0.0,tmp5);

                  rhs.clear();
                  Permute(tmp5,shape(0,2,1,4,3),rhs);

               }
               else{//right-down site of diagonal gate, so site (row,col + 1) environment

                  // (1) construct N_eff

                  //add upper-left peps to intermediate left
                  DArray<8> tmp8;
                  Contract(1.0,LI7,shape(1,3),peps(row+1,col),shape(0,3),0.0,tmp8);

                  //und again
                  DArray<5> tmp5;
                  Contract(1.0,tmp8,shape(0,5,6,1),peps(row+1,col),shape(0,1,2,3),0.0,tmp5);

                  //add next bottom environment to intermediate
                  DArray<7> tmp7;
                  Contract(1.0,tmp5,shape(2),env.gb(row-1)[col+1],shape(0),0.0,tmp7);

                  //add left and right together
                  tmp8.clear();
                  Contract(1.0,tmp7,shape(3,2,6),RI7,shape(0,1,6),0.0,tmp8);

                  N_eff.clear();
                  Permute(tmp8,shape(0,4,2,6,1,5,3,7),N_eff);

                  // (2) construct right hand side

                  //add upper-left peps to b_L
                  tmp8.clear();
                  Contract(1.0,b_L,shape(1,3),peps(row+1,col),shape(0,3),0.0,tmp8);

                  //next add left operator to tmp8
                  DArray<6> tmp6;
                  Contract(1.0,tmp8,shape(0,5,6,1),lop,shape(0,1,2,4),0.0,tmp6);

                  //add next bottom environment to intermediate
                  tmp8.clear();
                  Contract(1.0,tmp6,shape(2),env.gb(row-1)[col+1],shape(0),0.0,tmp8);

                  //and add right operator to tmp8
                  DArray<8> tmp8bis;
                  Contract(1.0,tmp8,shape(0,3,5),rop,shape(0,3,4),0.0,tmp8bis);

                  tmp5.clear();
                  Contract(1.0,tmp8bis,shape(2,1,5,7,4),RI7,shape(0,1,2,4,6),0.0,tmp5);

                  rhs.clear();
                  Permute(tmp5,shape(0,3,1,4,2),rhs);

               }

            }

         }
         else{//DIAGONAL LDRU

            if(row == 0){

               if(left){//left site of diagonal gate, so site (row,col) environment

                  // (1) construct N_eff
                  DArray<8> tmp8;
                  Contract(1.0,peps(row+1,col+1),shape(3,4),RI7,shape(4,2),0.0,tmp8);

                  //und again
                  DArray<7> tmp7;
                  Contract(1.0,peps(row+1,col+1),shape(2,3,4),tmp8,shape(2,5,4),0.0,tmp7);

                  //and add top environment to intermediate
                  DArray<5> tmp5;
                  Contract(1.0,env.gt(row)[col+1],shape(1,2,3),tmp7,shape(1,3,4),0.0,tmp5);

                  //add LI7 to it
                  DArray<6> tmp6;
                  Gemm(CblasTrans,CblasNoTrans,1.0,LI7,tmp5,0.0,tmp6);

                  DArray<6> tmp6bis;
                  Permute(tmp6,shape(2,0,4,3,1,5),tmp6bis);

                  int DL = peps(row,col).shape(0);
                  int DU = peps(row,col).shape(1);
                  int DD = peps(row,col).shape(3);
                  int DR = peps(row,col).shape(4);

                  N_eff = tmp6bis.reshape_clear( shape(DL,DU,DR,DD,DL,DU,DD,DR) );

                  // (2) construct right hand side

                  //add top right peps to b_R
                  Contract(1.0,peps(row+1,col+1),shape(3,4),b_R,shape(4,2),0.0,tmp8);

                  //add right operator to tmp8
                  DArray<8> tmp8bis;
                  Contract(1.0,rop,shape(2,4,5),tmp8,shape(2,5,4),0.0,tmp8bis);

                  //add top environment to intermediate
                  tmp6.clear();
                  Contract(1.0,env.gt(row)[col+1],shape(1,2,3),tmp8bis,shape(1,4,5),0.0,tmp6);

                  //add left operator to LI7
                  DArray<9> tmp9;
                  Contract(1.0,LI7,shape(5,3),lop,shape(0,1),0.0,tmp9);

                  //contract both sides to form right hand side of equation
                  tmp5.clear();
                  Contract(1.0,tmp9,shape(0,1,2,6,8),tmp6,shape(0,1,3,2,4),0.0,tmp5);

                  rhs.clear();
                  Permute(tmp5,shape(1,0,3,4,2),rhs);

               }
               else{//right site of diagonal gate, so site (row+1,col+1) environment

                  // (1) construct N_eff
                  DArray<8> tmp8;
                  Contract(1.0,LI7,shape(6,4),peps(row,col),shape(0,1),0.0,tmp8);

                  //and again
                  DArray<5> tmp5;
                  Contract(1.0,tmp8,shape(4,3,5,6),peps(row,col),shape(0,1,2,3),0.0,tmp5);

                  //now add top environemt to tmp5
                  DArray<7> tmp7;
                  Gemm(CblasTrans,CblasNoTrans,1.0,tmp5,env.gt(row)[col+1],0.0,tmp7);

                  //add RI7 to it
                  DArray<8> tmp8bis;
                  Contract(1.0,tmp7,shape(6,3,2),RI7,shape(0,5,6),0.0,tmp8bis);

                  N_eff.clear();
                  Permute(tmp8bis,shape(0,2,6,4,1,3,7,5),N_eff);

                  // (2) construct right hand side

                  //add left operator to tmp8
                  DArray<6> tmp6;
                  Contract(1.0,tmp8,shape(4,3,5,6),lop,shape(0,1,2,4),0.0,tmp6);

                  tmp8.clear();
                  Gemm(CblasTrans,CblasNoTrans,1.0,tmp6,env.gt(row)[col+1],0.0,tmp8);

                  //add right operator to b_R
                  DArray<9> tmp9;
                  Contract(1.0,rop,shape(4,5),b_R,shape(3,1),0.0,tmp9);

                  //contract both sides to form right hand side of equation
                  tmp5.clear();
                  Contract(1.0,tmp8,shape(7,5,0,4,3,2),tmp9,shape(4,1,0,7,3,8),0.0,tmp5);

                  rhs.clear();
                  Permute(tmp5,shape(0,1,4,3,2),rhs);

               }

            }
            else{//row = Ly - 2

               if(left){//left site of diagonal gate, so site (row,col) environment

                  // (1) construct N_eff
                  DArray<8> tmp8;
                  Contract(1.0,peps(row+1,col+1),shape(3,4),RI7,shape(3,1),0.0,tmp8);

                  //and again
                  DArray<5> tmp5;
                  Contract(1.0,peps(row+1,col+1),shape(1,2,3,4),tmp8,shape(1,2,4,3),0.0,tmp5);

                  //now add bottom environment to tmp5
                  DArray<7> tmp7;
                  Gemm(CblasNoTrans,CblasTrans,1.0,env.gb(row-1)[col],tmp5,0.0,tmp7);

                  //now contract left and right for N_eff construction
                  tmp8.clear();
                  Contract(1.0,LI7,shape(0,1,6),tmp7,shape(3,4,0),0.0,tmp8);

                  N_eff.clear();
                  Permute(tmp8,shape(2,0,4,6,3,1,5,7),N_eff);

                  // (2) construct right hand side: work with b_R instead of RI7
                  tmp8.clear();
                  Contract(1.0,peps(row+1,col+1),shape(3,4),b_R,shape(3,1),0.0,tmp8);

                  //and add right operator
                  DArray<6> tmp6;
                  Contract(1.0,rop,shape(1,2,4,5),tmp8,shape(1,2,4,3),0.0,tmp6);

                  //now add bottom environment to tmp6
                  tmp8.clear();
                  Gemm(CblasNoTrans,CblasTrans,1.0,env.gb(row-1)[col],tmp6,0.0,tmp8);

                  //finally add left operator
                  DArray<8> tmp8bis;
                  Contract(1.0,lop,shape(3,4,5),tmp8,shape(4,1,6),0.0,tmp8bis);

                  //contract left and righ
                  tmp5.clear();
                  Contract(1.0,LI7,shape(0,1,2,4,6),tmp8bis,shape(5,6,1,0,3),0.0,tmp5);

                  rhs.clear();
                  Permute(tmp5,shape(1,0,3,4,2),rhs);

               }
               else{//right site of diagonal gate, so site (row+1,col+1) environment

                  // (1) construct N_eff

                  //add bottom peps to LI7
                  DArray<8> tmp8;
                  Contract(1.0,LI7,shape(5,3),peps(row,col),shape(0,1),0.0,tmp8);

                  //and again
                  DArray<7> tmp7;
                  Contract(1.0,tmp8,shape(3,2,5),peps(row,col),shape(0,1,2),0.0,tmp7);

                  //now add bottom environment
                  DArray<5> tmp5;
                  Contract(1.0,tmp7,shape(2,5,3),env.gb(row-1)[col],shape(0,1,2),0.0,tmp5);

                  //contract left and right for N_eff construction
                  DArray<6> tmp6;
                  Contract(1.0,tmp5,shape(3,2,4),RI7,shape(4,5,6),0.0,tmp6);

                  DArray<6> tmp6bis;
                  Permute(tmp6,shape(0,4,2,1,5,3),tmp6bis);

                  int DR = peps(row+1,col+1).shape(4);

                  N_eff = tmp6bis.reshape_clear( shape(D,1,D,DR,D,1,D,DR) );

                  // (2) construct right hand side

                  //add left operator to intermediate
                  DArray<8> tmp8bis;
                  Contract(1.0,tmp8,shape(3,2,5),lop,shape(0,1,2),0.0,tmp8bis);

                  //now add bottom environment
                  tmp6.clear();
                  Contract(1.0,tmp8bis,shape(2,6,3),env.gb(row-1)[col],shape(0,1,2),0.0,tmp6);

                  //and add right operator
                  tmp8.clear();
                  Contract(1.0,tmp6,shape(0,3),rop,shape(0,3),0.0,tmp8);

                  //finally contract left with right
                  tmp5.clear();
                  Contract(1.0,tmp8,shape(7,6,2,1,3),b_R,shape(0,2,4,5,6),0.0,tmp5);

                  rhs.clear();
                  Permute(tmp5,shape(0,1,4,3,2),rhs);

               }

            }

         }

      }

   /**
    * construct the single-site effective environment and right hand side needed for the linear system any gate direction specified by 'dir'
    * @param dir vertical, horizontal,diagonal lurd or diagonal ldru update
    * @param row the row index of the bottom site
    * @param col column index of the vertical column
    * @param peps, full PEPS object before update
    * @param N_eff output object, contains N_eff on output
    * @param rhs output object, contains N_eff on output
    * @param LO Left environment contraction
    * @param RO Right environment contraction
    * @param LI8 left intermediate object
    * @param RI8 right intermediate object
    * @param b_L left intermediate object
    * @param b_R right intermediate object
    * @param left boolean flag for PEPS with left or right operator acted upon
    */
   template<>
      void construct_lin_sys(const PROP_DIR &dir,int row,int col,PEPS<double> &peps,const DArray<6> &lop,
            
            const DArray<6> &rop,DArray<8> &N_eff,DArray<5> &rhs, const DArray<6> &LO, const DArray<6> &RO,
            
            const DArray<8> &LI8,const DArray<8> &RI8,const DArray<8> &b_L,const DArray<8> &b_R,bool left){

         if(dir == VERTICAL){

            if(col == 0){

               if(left){//top site environment

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
               else{//bottom site

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

            }
            else if(col < Lx - 1){//col != 0

               if(left){//top site environment

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
               else{//bottom site

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

            }
            else{//col == Lx - 1

               if(left){//top site environment

                  // (1) calculate N_eff

                  //add top to intermediate
                  DArray<9> tmp9;
                  Contract(1.0,LI8,shape(1,7),peps(row+1,col),shape(0,1),0.0,tmp9);

                  //and another
                  DArray<8> tmp8;
                  Contract(1.0,tmp9,shape(0,5,6),peps(row+1,col),shape(0,1,2),0.0,tmp8);

                  N_eff.clear();
                  Permute(tmp8,shape(0,6,2,7,1,4,3,5),N_eff);

                  // (2) right hand side

                  //add right operator
                  DArray<9> tmp9bis;
                  Contract(1.0,tmp9,shape(0,5,6),rop,shape(0,1,2),0.0,tmp9bis);

                  //and left operator
                  DArray<5> tmp5;
                  Contract(1.0,tmp9bis,shape(0,7,6,2,8),lop,shape(0,1,3,4,5),0.0,tmp5);

                  rhs.clear();
                  Permute(tmp5,shape(0,2,1,3,4),rhs);

               }
               else{//bottom site

                  // (1) calculate N_eff

                  //add bottom peps  to intermediate
                  DArray<9> tmp9;
                  Contract(1.0,LI8,shape(3,5),peps(row,col),shape(0,3),0.0,tmp9);

                  //and another
                  DArray<8> tmp8;
                  Contract(1.0,tmp9,shape(2,7,3),peps(row,col),shape(0,2,3),0.0,tmp8);

                  N_eff.clear();
                  Permute(tmp8,shape(0,2,6,7,1,3,4,5),N_eff);

                  // (2) right hand side

                  //add left operator to intermediate
                  DArray<9> tmp9bis;
                  Contract(1.0,tmp9,shape(2,7,3),lop,shape(0,2,4),0.0,tmp9bis);

                  //and right operator
                  rhs.clear();
                  Contract(1.0,tmp9bis,shape(0,2,7,6,8),rop,shape(0,1,3,4,5),0.0,rhs);

               }

            }

         }//end VERTICAL
         else if(dir == HORIZONTAL){

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

         }//end horizontal
         else if(dir == DIAGONAL_LURD){

            if(left){//left site of gate, so site (row+1,col) environment

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

               //and top bottom environment on col
               tmp8.clear();
               Gemm(CblasNoTrans,CblasNoTrans,1.0,env.gt(row)[col],tmp6,0.0,tmp8);

               //contract with left hand side
               DArray<8> tmp8bis;
               Contract(1.0,LI8,shape(0,5,6,7),tmp8,shape(0,6,5,7),0.0,tmp8bis);

               N_eff.clear();
               Permute(tmp8bis,shape(0,4,2,6,1,5,3,7),N_eff);

               // (2) construct right hand side

               //add right operator
               DArray<9> tmp9bis;
               Contract(1.0,tmp9,shape(3,7,4),rop,shape(1,2,5),0.0,tmp9bis);

               //add bottom environment
               DArray<7> tmp7;
               Contract(1.0,tmp9bis,shape(8,5,3),env.gb(row-1)[col+1],shape(1,2,3),0.0,tmp7);

               //and top environment
               tmp9.clear();
               Gemm(CblasNoTrans,CblasNoTrans,1.0,env.gt(row)[col],tmp7,0.0,tmp9);

               //now add left operator
               tmp9bis.clear();
               Contract(1.0,tmp9,shape(1,7,3),lop,shape(1,3,5),0.0,tmp9bis);

               //attach b_L to right side
               DArray<5> tmp5;
               Contract(1.0,b_L,shape(0,1,3,5,6,7),tmp9bis,shape(0,6,8,4,3,5),0.0,tmp5);

               rhs.clear();
               Permute(tmp5,shape(0,2,1,3,4),rhs);

            }
            else{//right site of gate, so site (row,col+1) environment

               //(1) constsruct N_eff

               //add left-up peps to LI8
               DArray<9> tmp9;
               Contract(1.0,LI8,shape(2,4),peps(row+1,col),shape(0,3),0.0,tmp9);

               //and another peps
               DArray<8> tmp8;
               Contract(1.0,tmp9,shape(1,7,2),peps(row+1,col),shape(0,2,3),0.0,tmp8);

               //now contract with top environment
               DArray<6> tmp6;
               Contract(1.0,tmp8,shape(0,6,4),env.gt(row)[col],shape(0,1,2),0.0,tmp6);

               //and next bottom environment
               tmp8.clear();
               Contract(1.0,tmp6,shape(2),env.gb(row-1)[col+1],shape(0),0.0,tmp8);

               //finally contract with RI8
               DArray<8> tmp8bis;
               Contract(1.0,tmp8,shape(4,3,2,7),RI8,shape(0,1,2,7),0.0,tmp8bis);

               N_eff.clear();
               Permute(tmp8bis,shape(0,4,2,6,1,5,3,7),N_eff);

               // (2) construct right hand side

               //add left-up peps to b_L
               tmp9.clear();
               Contract(1.0,b_L,shape(2,4),peps(row+1,col),shape(0,3),0.0,tmp9);

               //add left operator to intermediate
               DArray<9> tmp9bis;
               Contract(1.0,tmp9,shape(1,7,2),lop,shape(0,2,4),0.0,tmp9bis);

               //now contract with top environment
               DArray<7> tmp7;
               Contract(1.0,tmp9bis,shape(0,6,4),env.gt(row)[col],shape(0,1,2),0.0,tmp7);

               //and next bottom environment
               tmp9.clear();
               Contract(1.0,tmp7,shape(2),env.gb(row-1)[col+1],shape(0),0.0,tmp9);

               //finally add right operator
               tmp9bis.clear();
               Contract(1.0,tmp9,shape(0,3,6),rop,shape(0,3,4),0.0,tmp9bis);

               //and contract with RI8 (same as b_R for lurd) to make right hand side
               DArray<5> tmp5;
               Contract(1.0,tmp9bis,shape(3,2,1,6,8,5),RI8,shape(0,1,2,3,5,7),0.0,tmp5);

               rhs.clear();
               Permute(tmp5,shape(0,3,1,4,2),rhs);

            }

         }
         else{//diagonal lrdu

            if(left){//left site of gate, so site (row,col) environment

               //(1) construct N_eff

               //add upper right peps to intermediate RI8
               DArray<9> tmp9;
               Contract(1.0,peps(row+1,col+1),shape(3,4),RI8,shape(4,2),0.0,tmp9);

               //and another
               DArray<8> tmp8;
               Contract(1.0,peps(row+1,col+1),shape(2,3,4),tmp9,shape(2,5,4),0.0,tmp8);

               //add top environment
               DArray<6> tmp6;
               Contract(1.0,env.gt(row)[col+1],shape(1,2,3),tmp8,shape(1,3,4),0.0,tmp6);

               //add bottom envirnoment
               tmp8.clear();
               Gemm(CblasNoTrans,CblasTrans,1.0,env.gb(row-1)[col],tmp6,0.0,tmp8);

               //now contract LI8 with tmp8
               DArray<8> tmp8bis;
               Contract(1.0,LI8,shape(0,1,2,7),tmp8,shape(3,4,5,0),0.0,tmp8bis);

               N_eff.clear();
               Permute(tmp8bis,shape(2,0,4,6,3,1,5,7),N_eff);

               // (2) construct right hand side

               //add upper right peps to intermediate b_R
               tmp9.clear();
               Contract(1.0,peps(row+1,col+1),shape(3,4),b_R,shape(4,2),0.0,tmp9);

               //and right operator to intermediate
               DArray<9> tmp9bis;
               Contract(1.0,rop,shape(2,4,5),tmp9,shape(2,5,4),0.0,tmp9bis);

               //add top environment
               DArray<7> tmp7;
               Contract(1.0,env.gt(row)[col+1],shape(1,2,3),tmp9bis,shape(1,4,5),0.0,tmp7);

               //add bottom envirnoment
               tmp9.clear();
               Gemm(CblasNoTrans,CblasTrans,1.0,env.gb(row-1)[col],tmp7,0.0,tmp9);

               //finally add left operator
               tmp9bis.clear();
               Contract(1.0,lop,shape(3,4,5),tmp9,shape(5,1,7),0.0,tmp9bis);

               //now attach to LI8 (which is b_L)
               DArray<5> tmp5;
               Contract(1.0,LI8,shape(0,1,2,3,5,7),tmp9bis,shape(5,6,7,1,0,3),0.0,tmp5);

               rhs.clear();
               Permute(tmp5,shape(1,0,3,4,2),rhs);

            }
            else{//right site of gate, so site (row+1,col+1) environment

               //(1) construct N_eff

               //add bottom left peps to intermediate LI8
               DArray<9> tmp9;
               Contract(1.0,LI8,shape(6,4),peps(row,col),shape(0,1),0.0,tmp9);

               //and another
               DArray<8> tmp8;
               Contract(1.0,tmp9,shape(4,3,6),peps(row,col),shape(0,1,2),0.0,tmp8);

               //add bottom environment
               DArray<6> tmp6;
               Contract(1.0,tmp8,shape(3,6,4),env.gb(row-1)[col],shape(0,1,2),0.0,tmp6);

               //add next top envirnoment
               tmp8.clear();
               Gemm(CblasTrans,CblasNoTrans,1.0,tmp6,env.gt(row)[col+1],0.0,tmp8);

               //now attach RI8 to tmp8 
               DArray<8> tmp8bis;
               Contract(1.0,tmp8,shape(7,3,2,4),RI8,shape(0,5,6,7),0.0,tmp8bis);

               N_eff.clear();
               Permute(tmp8bis,shape(0,2,6,4,1,3,7,5),N_eff);

               // (2) construct right hand side

               //add left operator to tmp9
               DArray<9> tmp9bis;
               Contract(1.0,tmp9,shape(4,3,6),lop,shape(0,1,2),0.0,tmp9bis);

               //add bottom environment
               DArray<7> tmp7;
               Contract(1.0,tmp9bis,shape(3,7,4),env.gb(row-1)[col],shape(0,1,2),0.0,tmp7);

               //add next top envirnoment
               tmp9.clear();
               Gemm(CblasTrans,CblasNoTrans,1.0,tmp7,env.gt(row)[col+1],0.0,tmp9);

               //add right operator
               tmp9bis.clear();
               Contract(1.0,tmp9,shape(0,6,3),rop,shape(0,1,3),0.0,tmp9bis);

               //connect tmp9bis with b_R to construct right hand side of equation
               DArray<5> tmp5;
               Contract(1.0,tmp9bis,shape(5,8,7,2,1,3),b_R,shape(0,1,3,5,6,7),0.0,tmp5);

               rhs.clear();
               Permute(tmp5,shape(0,1,4,3,2),rhs);

            }

         }

      }

   /**
    * first guess/ initialization of the peps pair by performing an SVD
    * @param dir vertical, horizontal,diagonal lurd or diagonal ldru update
    * @param row left bottom site row index
    * @param col left bottom site column index
    * @param lop peps acted onto with left trotter operator
    * @param rop peps acted onto with right trotter operator
    * @param peps full PEPS object, not const! relevant elements are changed
    */ 
   void initialize(const PROP_DIR &dir,int row,int col,const DArray<6> &lop,const DArray<6> &rop,PEPS<double> &peps){

      if(dir == VERTICAL){//row --> row+1

         DArray<8> tmp8;
         Contract(1.0,lop,shape(1,3),rop,shape(4,3),0.0,tmp8);

         //svd the fucker
         DArray<5> UL;//left unitary
         DArray<5> VR;//right unitary

         DArray<1> S;
         Gesvd ('S','S', tmp8, S,UL,VR,D);

         //take the square root of the sv's
         for(int i = 0;i < S.size();++i)
            S(i) = sqrt(S(i));

         //and multiply it left and right to the tensors
         Dimm(S,VR);
         Dimm(UL,S);

         //permute back to the peps
         Permute(UL,shape(0,4,1,2,3),peps(row,col));
         Permute(VR,shape(1,2,3,0,4),peps(row+1,col));

      }
      else if(dir == HORIZONTAL){//col --> col + 1

         DArray<8> tmp8;
         Contract(1.0,lop,shape(3,5),rop,shape(3,0),0.0,tmp8);

         //svd the fucker
         DArray<1> S;
         Gesvd ('S','S', tmp8, S,peps(row,col),peps(row,col+1),D);

         //take the square root of the sv's
         for(int i = 0;i < S.size();++i)
            S(i) = sqrt(S(i));

         //and multiply it left and right to the tensors
         Dimm(S,peps(row,col+1));
         Dimm(peps(row,col),S);

      }
      else if(dir == DIAGONAL_LURD){//(row+1,col) --> (row,col+1)

         //make a three-site object: connect lop with peps(row,col)
         DArray<9> tmp9;
         Contract(1.0,lop,shape(4),peps(row,col),shape(1),0.0,tmp9);

         //attach right operator to 
         DArray<11> tmp11;
         Contract(1.0,tmp9,shape(3,8),rop,shape(3,0),0.0,tmp11);

         //first split up in 2 - 1 part
         DArray<8> tmp8;//left unitary: 2-site part

         DArray<1> S;
         Gesvd ('S','S', tmp11, S,tmp8,peps(row,col+1),D);

         //take the square root of the sv's
         for(int i = 0;i < S.size();++i)
            S(i) = sqrt(S(i));

         //and multiply it left and right to the tensors
         Dimm(S,peps(row,col+1));
         Dimm(tmp8,S);

         //now just SVD the two-site part
         DArray<5> UL;//left unitary
         DArray<5> VR;//right unitary

         Gesvd ('S','S', tmp8, S,UL,VR,D);

         //take the square root of the sv's
         for(int i = 0;i < S.size();++i)
            S(i) = sqrt(S(i));

         //and multiply it left and right to the tensors
         Dimm(S,VR);
         Dimm(UL,S);

         //permute back to the peps
         Permute(UL,shape(0,1,2,4,3),peps(row+1,col));
         Permute(VR,shape(1,0,2,3,4),peps(row,col));

      }
      else{//diagonal LDRU

         //make a three-site object: connect lop with peps(row,col+1) (a.k.a. mop)
         DArray<9> tmp9;
         Contract(1.0,lop,shape(5),peps(row,col+1),shape(0),0.0,tmp9);

         //attach right operator to 
         DArray<11> tmp11;
         Contract(1.0,tmp9,shape(3,5),rop,shape(3,4),0.0,tmp11);

         //first split up in 1 site - 2 site part
         DArray<8> tmp8;//right unitary: 2-site part

         DArray<1> S;
         Gesvd ('S','S', tmp11, S,peps(row,col),tmp8,D);

         //take the square root of the sv's
         for(int i = 0;i < S.size();++i)
            S(i) = sqrt(S(i));

         //and multiply it left and right to the tensors
         Dimm(S,tmp8);
         Dimm(peps(row,col),S);

         //now just SVD the two-site part
         DArray<5> UL;//left unitary
         DArray<5> VR;//right unitary

         Gesvd ('S','S', tmp8, S,UL,VR,D);

         //take the square root of the sv's
         for(int i = 0;i < S.size();++i)
            S(i) = sqrt(S(i));

         //and multiply it left and right to the tensors
         Dimm(S,VR);
         Dimm(UL,S);

         //permute back to the peps
         Permute(UL,shape(0,4,1,2,3),peps(row,col+1));
         Permute(VR,shape(1,2,3,0,4),peps(row+1,col+1));

      }

   }

   /**
    * restore peps, i.e. put on equal footing after update is over, fix gauge
    * @param dir vertical, horizontal,diagonal lurd or diagonal ldru update
    * @param row left bottom site row index
    * @param col left bottom site column index
    * @param peps full PEPS object, not const! relevant elements are changed
    */ 
   void equilibrate(const PROP_DIR &dir,int row,int col,PEPS<double> &peps){

      if(dir == VERTICAL){

         DArray<8> tmp8;
         Contract(1.0,peps(row,col),shape(1),peps(row+1,col),shape(3),0.0,tmp8);

         //svd the fucker
         DArray<5> UL;//left unitary
         DArray<5> VR;//right unitary

         DArray<1> S;
         Gesvd ('S','S', tmp8, S,UL,VR,D);

         //take the square root of the sv's
         for(int i = 0;i < S.size();++i)
            S(i) = sqrt(S(i));

         //and multiply it left and right to the tensors
         Dimm(S,VR);
         Dimm(UL,S);

         //permute back to the peps
         Permute(UL,shape(0,4,1,2,3),peps(row,col));
         Permute(VR,shape(1,2,3,0,4),peps(row+1,col));

      }
      else if(dir == HORIZONTAL){

         DArray<8> tmp8;
         Contract(1.0,peps(row,col),shape(4),peps(row,col+1),shape(0),0.0,tmp8);

         //svd the fucker
         DArray<1> S;
         Gesvd ('S','S', tmp8, S,peps(row,col),peps(row,col+1),D);

         //take the square root of the sv's
         for(int i = 0;i < S.size();++i)
            S(i) = sqrt(S(i));

         //and multiply it left and right to the tensors
         Dimm(S,peps(row,col+1));
         Dimm(peps(row,col),S);

      }
      else if(dir == DIAGONAL_LURD){

         //make a three-site object: connect upper and lower peps
         DArray<8> tmp8;
         Contract(1.0,peps(row+1,col),shape(3),peps(row,col),shape(1),0.0,tmp8);

         //attach right peps to it
         DArray<11> tmp11;
         Contract(1.0,tmp8,shape(7),peps(row,col+1),shape(0),0.0,tmp11);

         //first split up in 2 - 1 part
         DArray<5> tmp5;

         DArray<1> S;
         Gesvd ('S','S', tmp11, S,tmp8,peps(row,col+1),D);

         //take the square root of the sv's
         for(int i = 0;i < S.size();++i)
            S(i) = sqrt(S(i));

         //and multiply it left and right to the tensors
         Dimm(S,peps(row,col+1));
         Dimm(tmp8,S);

         //now just SVD the two-site part
         DArray<5> UL;//left unitary
         DArray<5> VR;//right unitary

         Gesvd ('S','S', tmp8, S,UL,VR,D);

         //take the square root of the sv's
         for(int i = 0;i < S.size();++i)
            S(i) = sqrt(S(i));

         //and multiply it left and right to the tensors
         Dimm(S,VR);
         Dimm(UL,S);

         //permute back to the peps
         Permute(UL,shape(0,1,2,4,3),peps(row+1,col));
         Permute(VR,shape(1,0,2,3,4),peps(row,col));

      }
      else{//diagonal ldru

         //make a three-site object: connect peps(row,col) with peps(row,col+1) (a.k.a. mop)
         DArray<8> tmp8;
         Contract(1.0,peps(row,col),shape(4),peps(row,col+1),shape(0),0.0,tmp8);

         //attach peps(row+1,col+1) to it
         DArray<11> tmp11;
         Contract(1.0,tmp8,shape(4),peps(row+1,col+1),shape(3),0.0,tmp11);

         //first split up in 1 site - 2 site part

         DArray<1> S;
         Gesvd ('S','S', tmp11, S,peps(row,col),tmp8,D);

         //take the square root of the sv's
         for(int i = 0;i < S.size();++i)
            S(i) = sqrt(S(i));

         //and multiply it left and right to the tensors
         Dimm(S,tmp8);
         Dimm(peps(row,col),S);

         //now just SVD the two-site part
         DArray<5> UL;//left unitary
         DArray<5> VR;//right unitary

         Gesvd ('S','S', tmp8, S,UL,VR,D);

         //take the square root of the sv's
         for(int i = 0;i < S.size();++i)
            S(i) = sqrt(S(i));

         //and multiply it left and right to the tensors
         Dimm(S,VR);
         Dimm(UL,S);

         //permute back to the peps
         Permute(UL,shape(0,4,1,2,3),peps(row,col+1));
         Permute(VR,shape(1,2,3,0,4),peps(row+1,col+1));

      }

   }

   /**
    * function that calculates intermediate objects that do not change during the sweeping update.
    * By precalculating them a lot of work is avoided
    * @param dir vertical, horizontal,diagonal lurd or diagonal ldru update
    * @param row row index of the bottom peps of the vertical pair
    * @param col column index
    * @param peps full PEPS object
    * @param mop 'middle' peps, connecting the two peps to be updated for diagonal gates
    * @param L left contracted environment around the pair
    * @param R right contracted environment around the pair
    * @param LI7 Left intermediate object to be constructed on output
    * @param RI7 Right intermediate object to be constructed on output
    * @param b_L left intermediate object
    * @param b_R right intermediate object
    */
   template<>
      void construct_intermediate(const PROP_DIR &dir,int row,int col,const PEPS<double> &peps,const DArray<5> &mop,

            const DArray<5> &L,const DArray<5> &R,DArray<7> &LI7,DArray<7> &RI7,DArray<7> &b_L,DArray<7> &b_R){

         if(dir == VERTICAL){

            if(row == 0){

               if(col == 0)
                  Gemm(CblasNoTrans,CblasNoTrans,1.0,env.gt(0)[col],R,0.0,RI7);
               else if(col < Lx - 1)
                  Gemm(CblasTrans,CblasNoTrans,1.0,L,env.gt(0)[col],0.0,LI7);
               else
                  Gemm(CblasTrans,CblasNoTrans,1.0,L,env.gt(0)[col],0.0,LI7);

            }
            else{//row == Lx - 2

               if(col == 0)
                  Gemm(CblasNoTrans,CblasTrans,1.0,env.gb(Lx - 3)[col],R,0.0,RI7);
               else if(col < Lx - 1)
                  Gemm(CblasNoTrans,CblasTrans,1.0,env.gb(Lx - 3)[col],R,0.0,RI7);
               else
                  Gemm(CblasNoTrans,CblasNoTrans,1.0,L,env.gb(Lx - 3)[col],0.0,LI7);

            }

         }
         else if(dir == HORIZONTAL){

            if(row == 0){

               if(col == 0){

                  //right intermediate
                  Gemm(CblasNoTrans,CblasNoTrans,1.0,env.gt(0)[col+1],R,0.0,RI7);

                  DArray<8> tmp8;
                  Contract(1.0,RI7,shape(1,3),peps(row+1,col+1),shape(1,4),0.0,tmp8);

                  DArray<7> tmp7;
                  Contract(1.0,tmp8,shape(1,6,2),peps(row+1,col+1),shape(1,2,4),0.0,tmp7);

                  Permute(tmp7,shape(0,3,5,4,6,1,2),RI7);

                  //left intermediate
                  DArray<5> tmp5;
                  Gemm(CblasTrans,CblasNoTrans,1.0,env.gt(0)[col],peps(row+1,col),0.0,tmp5);

                  DArray<6> tmp6;
                  Contract(1.0,tmp5,shape(0,2),peps(row+1,col),shape(1,2),0.0,tmp6);

                  DArray<6> tmp6bis;
                  Permute(tmp6,shape(0,2,5,1,4,3),tmp6bis);

                  LI7 = tmp6bis.reshape_clear( shape(env.gt(0)[col].shape(3),D,D,D,D,1,1) );

               }
               else if(col < Lx - 2){

                  //right
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

               }
               else{//col == Lx - 2

                  //right
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

               }

            }
            else if(row == Ly - 2){

               if(col == 0){

                  //add top peps to right
                  DArray<8> tmp8;
                  Contract(1.0,peps(row+1,col+1),shape(4),R,shape(0),0.0,tmp8);

                  //and another
                  DArray<7> tmp7;
                  Contract(1.0,peps(row+1,col+1),shape(1,2,4),tmp8,shape(1,2,4),0.0,tmp7);

                  RI7.clear();
                  Permute(tmp7,shape(2,0,3,1,4,5,6),RI7);

                  //left
                  DArray<4> tmp4;
                  Gemm(CblasTrans,CblasNoTrans,1.0,peps(row+1,col),peps(row+1,col),0.0,tmp4);

                  DArray<4> tmp4bis;
                  Permute(tmp4,shape(1,3,0,2),tmp4bis);

                  LI7 = tmp4bis.reshape_clear( shape(D,D,D,D,1,1,1) );

               }
               else if(col < Lx - 2){

                  //add top peps to right
                  DArray<8> tmp8;
                  Contract(1.0,peps(row+1,col+1),shape(4),R,shape(0),0.0,tmp8);

                  //and another
                  DArray<7> tmp7;
                  Contract(1.0,peps(row+1,col+1),shape(1,2,4),tmp8,shape(1,2,4),0.0,tmp7);

                  RI7.clear();
                  Permute(tmp7,shape(2,0,3,1,4,5,6),RI7);

                  //left

                  //add top peps to left
                  tmp8.clear();
                  Contract(1.0,L,shape(0),peps(row+1,col),shape(0),0.0,tmp8);

                  tmp7.clear();
                  Contract(1.0,tmp8,shape(0,4,5),peps(row+1,col),shape(0,1,2),0.0,tmp7);

                  LI7.clear();
                  Permute(tmp7,shape(4,6,3,5,0,1,2),LI7);

               }
               else{//col == Lx - 2

                  //right
                  DArray<4> tmp4;
                  Contract(1.0,peps(row+1,col+1),shape(1,2,4),peps(row+1,col+1),shape(1,2,4),0.0,tmp4);

                  DArray<4> tmp4bis;
                  Permute(tmp4,shape(0,2,1,3),tmp4bis);

                  RI7 = tmp4bis.reshape_clear( shape(D,D,D,D,1,1,1) );

                  //left

                  //add top peps to left
                  DArray<8> tmp8;
                  Contract(1.0,L,shape(0),peps(row+1,col),shape(0),0.0,tmp8);

                  DArray<7> tmp7;
                  Contract(1.0,tmp8,shape(0,4,5),peps(row+1,col),shape(0,1,2),0.0,tmp7);

                  Permute(tmp7,shape(4,6,3,5,0,1,2),LI7);

               }

            }
            else{//row == Ly -1 

               if(col == 0){

                  //right

                  //add bottom env to right
                  Gemm(CblasNoTrans,CblasTrans,1.0,env.gb(Ly-3)[col+1],R,0.0,RI7);

                  //add bottom peps to intermediate
                  DArray<8> tmp8;
                  Contract(1.0,peps(row-1,col+1),shape(3,4),RI7,shape(2,6),0.0,tmp8);

                  //add second bottom peps
                  DArray<7> tmp7;
                  Contract(1.0,peps(row-1,col+1),shape(2,3,4),tmp8,shape(2,4,7),0.0,tmp7);

                  RI7.clear();
                  Permute(tmp7,shape(5,6,1,3,0,2,4),RI7);

                  //left

                  //add bottom peps to bottom env
                  DArray<5> tmp5;
                  Contract(1.0,peps(row-1,col),shape(0,3),env.gb(Ly-3)[col],shape(0,1),0.0,tmp5);

                  //and again
                  DArray<6> tmp6;
                  Contract(1.0,peps(row-1,col),shape(2,3),tmp5,shape(1,3),0.0,tmp6);

                  DArray<6> tmp6bis;
                  Permute(tmp6,shape(0,3,1,4,2,5),tmp6bis);

                  LI7 = tmp6bis.reshape_clear( shape(1,1,D,D,D,D,env.gb(Ly-3)[col].shape(3)) );

               }
               else if(col < Lx - 2){

                  //right

                  //add bottom env to right
                  Gemm(CblasNoTrans,CblasTrans,1.0,env.gb(Ly-3)[col+1],R,0.0,RI7);

                  //add bottom peps to intermediate
                  DArray<8> tmp8;
                  Contract(1.0,peps(row-1,col+1),shape(3,4),RI7,shape(2,6),0.0,tmp8);

                  //add second bottom peps
                  DArray<7> tmp7;
                  Contract(1.0,peps(row-1,col+1),shape(2,3,4),tmp8,shape(2,4,7),0.0,tmp7);

                  RI7.clear();
                  Permute(tmp7,shape(5,6,1,3,0,2,4),RI7);

                  //left

                  //add bottom env to left
                  Gemm(CblasNoTrans,CblasNoTrans,1.0,L,env.gb(Ly-3)[col],0.0,LI7);

                  //add bottom peps to intermediate
                  tmp8.clear();
                  Contract(1.0,LI7,shape(2,4),peps(row-1,col),shape(0,3),0.0,tmp8);

                  //add second bottom peps
                  tmp7.clear();
                  Contract(1.0,tmp8,shape(2,6,3),peps(row-1,col),shape(0,2,3),0.0,tmp7);

                  LI7.clear();
                  Permute(tmp7,shape(0,1,3,5,4,6,2),LI7);

               }
               else{//col == Lx - 2

                  //right

                  //attach bottom peps to bottom env
                  DArray<5> tmp5;
                  Contract(1.0,peps(row-1,col+1),shape(3,4),env.gb(Ly-3)[col+1],shape(2,3),0.0,tmp5);

                  //and another
                  DArray<6> tmp6;
                  Contract(1.0,peps(row-1,col+1),shape(2,3),tmp5,shape(2,4),0.0,tmp6);

                  DArray<6> tmp6bis;
                  Permute(tmp6,shape(2,1,4,0,3,5),tmp6bis);

                  RI7 = tmp6bis.reshape_clear( shape(1,1,D,D,D,D,env.gb(Ly-3)[col+1].shape(0)) );

                  //left

                  //add bottom env to left
                  Gemm(CblasNoTrans,CblasNoTrans,1.0,L,env.gb(Ly-3)[col],0.0,LI7);

                  //add bottom peps to intermediate
                  DArray<8> tmp8;
                  Contract(1.0,LI7,shape(2,4),peps(row-1,col),shape(0,3),0.0,tmp8);

                  //add second bottom peps
                  DArray<7> tmp7;
                  Contract(1.0,tmp8,shape(2,6,3),peps(row-1,col),shape(0,2,3),0.0,tmp7);

                  LI7.clear();
                  Permute(tmp7,shape(0,1,3,5,4,6,2),LI7);

               }

            }

         }
         else if(dir == DIAGONAL_LURD){

            if(row == 0){

               if(col == 0){

                  //create left and right intermediary operators: right

                  //attach top environment to right side
                  DArray<7> tmp7;
                  Contract(1.0,env.gt(row)[col + 1],shape(3),R,shape(0),0.0,tmp7);

                  //add upper right peps to it
                  DArray<8> tmp8;
                  Contract(1.0,tmp7,shape(1,3),peps(row+1,col+1),shape(1,4),0.0,tmp8);

                  //and another
                  tmp7.clear();
                  Contract(1.0,tmp8,shape(1,6,2),peps(row+1,col+1),shape(1,2,4),0.0,tmp7);

                  Permute(tmp7,shape(0,3,5,4,6,1,2),RI7);

                  //b_R is just equal to RI7 for LURD! so leave empty

                  //left: connect bottom left peps to itself
                  DArray<4> tmp4;
                  Contract(1.0,peps(row,col),shape(0,2,3),peps(row,col),shape(0,2,3),0.0,tmp4);

                  DArray<4> tmp4bis;
                  Permute(tmp4,shape(0,2,1,3),tmp4bis);

                  LI7 = tmp4bis.reshape_clear( shape(1,1,1,D,D,D,D) );

                  //for b_L contract mop with peps(row,col)
                  Contract(1.0,mop,shape(0,2,3),peps(row,col),shape(0,2,3),0.0,tmp4);

                  Permute(tmp4,shape(0,2,1,3),tmp4bis);

                  b_L = tmp4bis.reshape_clear( shape(1,1,1,D,D,D,D) );

               }
               else if(col < Lx - 2){

                  //create left and right intermediary operators: right

                  //attach top environment to right side
                  DArray<7> tmp7;
                  Contract(1.0,env.gt(row)[col + 1],shape(3),R,shape(0),0.0,tmp7);

                  //add upper right peps to it
                  DArray<8> tmp8;
                  Contract(1.0,tmp7,shape(1,3),peps(row+1,col+1),shape(1,4),0.0,tmp8);

                  //and another
                  tmp7.clear();
                  Contract(1.0,tmp8,shape(1,6,2),peps(row+1,col+1),shape(1,2,4),0.0,tmp7);

                  Permute(tmp7,shape(0,3,5,4,6,1,2),RI7);

                  //b_R is just equal to RI7 for LURD! so leave empty

                  //left: connect bottom left peps to L
                  tmp8.clear();
                  Contract(1.0,L,shape(4),peps(row,col),shape(0),0.0,tmp8);

                  //add another peps to make LI7
                  tmp7.clear();
                  Contract(1.0,tmp8,shape(3,5,6),peps(row,col),shape(0,2,3),0.0,tmp7);

                  LI7.clear();
                  Permute(tmp7,shape(0,1,2,5,3,6,4),LI7);

                  //for b_L contract tmp8 with mop
                  tmp7.clear();
                  Contract(1.0,tmp8,shape(3,5,6),mop,shape(0,2,3),0.0,tmp7);

                  b_L.clear();
                  Permute(tmp7,shape(0,1,2,5,3,6,4),b_L);

               }
               else{//col == Lx - 2

                  //create left and right intermediary operators: right

                  //attach top environment to upper peps
                  DArray<5> tmp5;
                  Contract(1.0,env.gt(row)[col + 1],shape(2,3),peps(row+1,col+1),shape(1,4),0.0,tmp5);

                  //add another one to it to form RI7
                  DArray<6> tmp6;
                  Contract(1.0,tmp5,shape(1,3),peps(row+1,col+1),shape(1,2),0.0,tmp6);

                  DArray<6> tmp6bis;
                  Permute(tmp6,shape(0,3,1,4,2,5),tmp6bis);

                  RI7 = tmp6bis.reshape_clear( shape(env.gt(row)[col+1].shape(0),D,D,D,D,1,1) );

                  //b_R is just equal to RI7 for LURD! so leave empty

                  //left: connect bottom left peps to L
                  DArray<8> tmp8;
                  Contract(1.0,L,shape(4),peps(row,col),shape(0),0.0,tmp8);

                  //add another peps to make LI7
                  DArray<7> tmp7;
                  Contract(1.0,tmp8,shape(3,5,6),peps(row,col),shape(0,2,3),0.0,tmp7);

                  LI7.clear();
                  Permute(tmp7,shape(0,1,2,5,3,6,4),LI7);

                  //for b_L contract tmp8 with mop
                  tmp7.clear();
                  Contract(1.0,tmp8,shape(3,5,6),mop,shape(0,2,3),0.0,tmp7);

                  b_L.clear();
                  Permute(tmp7,shape(0,1,2,5,3,6,4),b_L);

               }

            }
            else{// row == Ly - 2

               if(col == 0){

                  //add top peps to right
                  DArray<8> tmp8;
                  Contract(1.0,peps(row+1,col+1),shape(4),R,shape(0),0.0,tmp8);

                  //and another
                  DArray<7> tmp7;
                  Contract(1.0,peps(row+1,col+1),shape(1,2,4),tmp8,shape(1,2,4),0.0,tmp7);

                  RI7.clear();
                  Permute(tmp7,shape(2,0,3,1,4,5,6),RI7);

                  //b_R is equal to RI7 for LURD

                  //left

                  //add bottom peps to bottom env
                  DArray<5> tmp5;
                  Contract(1.0,peps(row,col),shape(0,3),env.gb(row-1)[col],shape(0,2),0.0,tmp5);

                  //and again
                  DArray<6> tmp6;
                  Contract(1.0,peps(row,col),shape(2,3),tmp5,shape(1,3),0.0,tmp6);

                  DArray<6> tmp6bis;
                  Permute(tmp6,shape(0,1,3,2,4,5),tmp6bis);

                  LI7 = tmp6bis.reshape_clear( shape(1,1,D,D,D,D,env.gb(Ly-3)[col].shape(3)) );

                  //add mop to tmp5 to construct b_L
                  Contract(1.0,mop,shape(2,3),tmp5,shape(1,3),0.0,tmp6);

                  Permute(tmp6,shape(0,1,3,2,4,5),tmp6bis);

                  b_L = tmp6bis.reshape_clear( shape(1,1,D,D,D,D,env.gb(Ly-3)[col].shape(3)) );

               }
               else if(col < Lx - 2){

                  //right

                  //add top peps to right
                  DArray<8> tmp8;
                  Contract(1.0,peps(row+1,col+1),shape(4),R,shape(0),0.0,tmp8);

                  //and another
                  DArray<7> tmp7;
                  Contract(1.0,peps(row+1,col+1),shape(1,2,4),tmp8,shape(1,2,4),0.0,tmp7);

                  RI7.clear();
                  Permute(tmp7,shape(2,0,3,1,4,5,6),RI7);

                  //b_R is equal to RI7 for LURD

                  //left

                  //add bottom env to left
                  Gemm(CblasNoTrans,CblasNoTrans,1.0,L,env.gb(Ly-3)[col],0.0,LI7);

                  //add bottom peps to intermediate
                  tmp8.clear();
                  Contract(1.0,LI7,shape(3,5),peps(row,col),shape(0,3),0.0,tmp8);

                  //add second bottom peps
                  tmp7.clear();
                  Contract(1.0,tmp8,shape(2,6,3),peps(row,col),shape(0,2,3),0.0,tmp7);

                  LI7.clear();
                  Permute(tmp7,shape(0,1,5,3,6,4,2),LI7);

                  //for b_L construction, add mop to tmp8
                  Contract(1.0,tmp8,shape(2,6,3),mop,shape(0,2,3),0.0,tmp7);

                  Permute(tmp7,shape(0,1,5,3,6,4,2),b_L);

               }
               else{//col == Lx - 2

                  //right
                  DArray<4> tmp4;
                  Contract(1.0,peps(row+1,col+1),shape(1,2,4),peps(row+1,col+1),shape(1,2,4),0.0,tmp4);

                  DArray<4> tmp4bis;
                  Permute(tmp4,shape(0,2,1,3),tmp4bis);

                  RI7 = tmp4bis.reshape_clear( shape(D,D,D,D,1,1,1) );

                  //no b_R

                  //left

                  //add bottom env to left
                  LI7.clear();
                  Gemm(CblasNoTrans,CblasNoTrans,1.0,L,env.gb(Ly-3)[col],0.0,LI7);

                  //add bottom peps to intermediate
                  DArray<8> tmp8;
                  Contract(1.0,LI7,shape(3,5),peps(row,col),shape(0,3),0.0,tmp8);

                  //add second bottom peps
                  DArray<7> tmp7;
                  Contract(1.0,tmp8,shape(2,6,3),peps(row,col),shape(0,2,3),0.0,tmp7);

                  LI7.clear();
                  Permute(tmp7,shape(0,1,5,3,6,4,2),LI7);

                  //for b_L construction, add mop to tmp8
                  Contract(1.0,tmp8,shape(2,6,3),mop,shape(0,2,3),0.0,tmp7);

                  Permute(tmp7,shape(0,1,5,3,6,4,2),b_L);

               }

            }

         }
         else{//DIAGONAL LDRU

            if(row == 0){

               if(col == 0){

                  //create left and right intermediary operators: right

                  //add bottom-right peps to right side
                  DArray<8> tmp8;
                  Contract(1.0,peps(row,col+1),shape(4),R,shape(4),0.0,tmp8);

                  //and another one for RI7
                  DArray<7> tmp7;
                  Contract(1.0,peps(row,col+1),shape(2,3,4),tmp8,shape(2,3,7),0.0,tmp7);

                  //set in the right order
                  Permute(tmp7,shape(4,5,6,1,3,0,2),RI7);

                  //now add mop to tmp8 for b_R construction
                  Contract(1.0,mop,shape(2,3,4),tmp8,shape(2,3,7),0.0,tmp7);

                  //set in the right order
                  Permute(tmp7,shape(4,5,6,1,3,0,2),b_R);

                  //left: connect top environment to top right peps
                  DArray<5> tmp5;
                  Contract(1.0,env.gt(row)[col],shape(0,1),peps(row+1,col),shape(0,1),0.0,tmp5);

                  DArray<6> tmp6;
                  Contract(1.0,tmp5,shape(0,2),peps(row+1,col),shape(1,2),0.0,tmp6);

                  DArray<6> tmp6bis;
                  Permute(tmp6,shape(0,2,5,1,4,3),tmp6bis);

                  LI7 = tmp6bis.reshape_clear( shape(env.gt(row)[col].shape(3),D,D,D,D,1,1) );

               }
               else if(col < Lx - 2){

                  //create left and right intermediary operators: right

                  //add bottom-right peps to right side
                  DArray<8> tmp8;
                  Contract(1.0,peps(row,col+1),shape(4),R,shape(4),0.0,tmp8);

                  //and another one for RI7
                  DArray<7> tmp7;
                  Contract(1.0,peps(row,col+1),shape(2,3,4),tmp8,shape(2,3,7),0.0,tmp7);

                  //set in the right order
                  Permute(tmp7,shape(4,5,6,1,3,0,2),RI7);

                  //now add mop to tmp8 for b_R construction
                  Contract(1.0,mop,shape(2,3,4),tmp8,shape(2,3,7),0.0,tmp7);

                  //set in the right order
                  Permute(tmp7,shape(4,5,6,1,3,0,2),b_R);

                  //left: connect top environment to L
                  tmp7.clear();
                  Gemm(CblasTrans,CblasNoTrans,1.0,L,env.gt(row)[col],0.0,tmp7);

                  //add upper peps
                  tmp8.clear();
                  Contract(1.0,tmp7,shape(0,4),peps(row+1,col),shape(0,1),0.0,tmp8);

                  tmp7.clear();
                  Contract(1.0,tmp8,shape(0,3,5),peps(row+1,col),shape(0,1,2),0.0,tmp7);

                  LI7.clear();
                  Permute(tmp7,shape(2,4,6,3,5,0,1),LI7);

                  //b_L is identical to LI7 for LDRU, so leave empty

               }
               else{//col == Lx - 2

                  //create left and right intermediary operators: right

                  //connect bottom right peps with itsself
                  DArray<4> tmp4;
                  Gemm(CblasNoTrans,CblasTrans,1.0,peps(row,col+1),peps(row,col+1),0.0,tmp4);

                  DArray<4> tmp4bis;
                  Permute(tmp4,shape(1,3,0,2),tmp4bis);

                  RI7 = tmp4bis.reshape_clear( shape(1,1,1,D,D,D,D) );

                  //connect mop with bottom right peps for b_R
                  Gemm(CblasNoTrans,CblasTrans,1.0,mop,peps(row,col+1),0.0,tmp4);

                  Permute(tmp4,shape(1,3,0,2),tmp4bis);

                  b_R = tmp4bis.reshape_clear( shape(1,1,1,D,D,D,D) );

                  //left: connect top environment to L
                  DArray<7> tmp7;
                  Gemm(CblasTrans,CblasNoTrans,1.0,L,env.gt(row)[col],0.0,tmp7);

                  //add upper peps
                  DArray<8> tmp8;
                  Contract(1.0,tmp7,shape(0,4),peps(row+1,col),shape(0,1),0.0,tmp8);

                  tmp7.clear();
                  Contract(1.0,tmp8,shape(0,3,5),peps(row+1,col),shape(0,1,2),0.0,tmp7);

                  LI7.clear();
                  Permute(tmp7,shape(2,4,6,3,5,0,1),LI7);

                  //b_L is identical to LI7 for LDRU, so leave empty

               }

            }
            else{// row == Ly - 2)

               if(col == 0){

                  //right

                  //add bottom env to right
                  Gemm(CblasNoTrans,CblasTrans,1.0,env.gb(Ly-3)[col+1],R,0.0,RI7);

                  //add bottom peps to intermediate
                  DArray<8> tmp8;
                  Contract(1.0,peps(row,col+1),shape(3,4),RI7,shape(2,6),0.0,tmp8);

                  //add second bottom peps
                  DArray<7> tmp7;
                  Contract(1.0,peps(row,col+1),shape(2,3,4),tmp8,shape(2,4,7),0.0,tmp7);

                  RI7.clear();
                  Permute(tmp7,shape(5,6,1,3,0,2,4),RI7);

                  //for b_R construction add mop to tmp8
                  Contract(1.0,mop,shape(2,3,4),tmp8,shape(2,4,7),0.0,tmp7);

                  b_R.clear();
                  Permute(tmp7,shape(5,6,1,3,0,2,4),b_R);

                  //left, b_L is equal to LI7
                  DArray<4> tmp4;
                  Gemm(CblasTrans,CblasNoTrans,1.0,peps(row+1,col),peps(row+1,col),0.0,tmp4);

                  DArray<4> tmp4bis;
                  Permute(tmp4,shape(1,3,0,2),tmp4bis);

                  LI7 = tmp4bis.reshape_clear( shape(D,D,D,D,1,1,1) );

               }
               else if(col < Lx - 2){

                  //right

                  //add bottom env to right
                  Gemm(CblasNoTrans,CblasTrans,1.0,env.gb(Ly-3)[col+1],R,0.0,RI7);

                  //add bottom peps to intermediate
                  DArray<8> tmp8;
                  Contract(1.0,peps(row,col+1),shape(3,4),RI7,shape(2,6),0.0,tmp8);

                  //add second bottom peps
                  DArray<7> tmp7;
                  Contract(1.0,peps(row,col+1),shape(2,3,4),tmp8,shape(2,4,7),0.0,tmp7);

                  RI7.clear();
                  Permute(tmp7,shape(5,6,1,3,0,2,4),RI7);

                  //for b_R construction add mop to tmp8
                  Contract(1.0,mop,shape(2,3,4),tmp8,shape(2,4,7),0.0,tmp7);

                  b_R.clear();
                  Permute(tmp7,shape(5,6,1,3,0,2,4),b_R);

                  //left

                  //add top peps to left
                  tmp8.clear();
                  Contract(1.0,L,shape(0),peps(row+1,col),shape(0),0.0,tmp8);

                  tmp7.clear();
                  Contract(1.0,tmp8,shape(0,4,5),peps(row+1,col),shape(0,1,2),0.0,tmp7);

                  LI7.clear();
                  Permute(tmp7,shape(4,6,3,5,0,1,2),LI7);

               }
               else{//col == Lx - 2

                  //right

                  //attach bottom peps to bottom env
                  DArray<5> tmp5;
                  Contract(1.0,peps(row,col+1),shape(3,4),env.gb(Ly-3)[col+1],shape(2,3),0.0,tmp5);

                  //and another
                  DArray<6> tmp6;
                  Contract(1.0,peps(row,col+1),shape(2,3),tmp5,shape(2,4),0.0,tmp6);

                  DArray<6> tmp6bis;
                  Permute(tmp6,shape(2,1,4,0,3,5),tmp6bis);

                  RI7 = tmp6bis.reshape_clear( shape(1,1,D,D,D,D,env.gb(Ly-3)[col+1].shape(0)) );

                  //for b_R add mop to tmp5
                  Contract(1.0,mop,shape(2,3),tmp5,shape(2,4),0.0,tmp6);

                  Permute(tmp6,shape(2,1,4,0,3,5),tmp6bis);

                  b_R = tmp6bis.reshape_clear( shape(1,1,D,D,D,D,env.gb(Ly-3)[col+1].shape(0)) );

                  //left

                  //add top peps to left
                  DArray<8> tmp8;
                  Contract(1.0,L,shape(0),peps(row+1,col),shape(0),0.0,tmp8);

                  DArray<7> tmp7;
                  Contract(1.0,tmp8,shape(0,4,5),peps(row+1,col),shape(0,1,2),0.0,tmp7);

                  LI7.clear();
                  Permute(tmp7,shape(4,6,3,5,0,1,2),LI7);

               }

            }

         }

      }

   /**
    * function that calculates intermediate objects that do not change during the sweeping update.
    * By precalculating them a lot of work is avoided
    * @param dir vertical, horizontal,diagonal lurd or diagonal ldru update
    * @param row row index of the bottom peps of the vertical pair
    * @param col column index
    * @param peps full PEPS object
    * @param mop 'middle' peps, connecting the two peps to be updated for diagonal gates
    * @param LO left contracted environment around the pair
    * @param RO right contracted environment around the pair
    * @param LI8 Left intermediate object to be constructed on output
    * @param RI8 Right intermediate object to be constructed on output
    * @param b_L left intermediate object
    * @param b_R right intermediate object
    */
   template<>
      void construct_intermediate(const PROP_DIR &dir,int row,int col,const PEPS<double> &peps,const DArray<5> &mop,

            const DArray<6> &LO,const DArray<6> &RO,DArray<8> &LI8,DArray<8> &RI8,DArray<8> &b_L,DArray<8> &b_R){

         if(dir == VERTICAL){

            if(col == 0){

               DArray<8> tmp8;
               Gemm(CblasNoTrans,CblasNoTrans,1.0,env.gt(row)[col],RO,0.0,tmp8);

               //inefficient but it's on the side, so it doesn't matter
               Contract(1.0,tmp8,shape(0,7),env.gb(row-1)[col],shape(0,3),0.0,RI8);

            } 
            else if(col < Lx - 1){

               //right
               Gemm(CblasNoTrans,CblasTrans,1.0,env.gb(row-1)[col],RO,0.0,RI8);

               //left
               Gemm(CblasTrans,CblasNoTrans,1.0,LO,env.gt(row)[col],0.0,LI8);

            }
            else{

               //first add bottom to left
               Gemm(CblasNoTrans,CblasNoTrans,1.0,LO,env.gb(row-1)[Lx-1],0.0,LI8);

               //then top to construct LI8
               DArray<10> tmp10;
               Gemm(CblasTrans,CblasNoTrans,1.0,LI8,env.gt(row)[Lx-1],0.0,tmp10);

               LI8 = tmp10.reshape_clear(shape(D,D,D,D,D,D,D,D));

            }

         }
         else if(dir == HORIZONTAL){

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

            }
            else if(col < Lx - 2){

               //create left and right intermediary operators: right
               Gemm(CblasNoTrans,CblasNoTrans,1.0,env.gt(row)[col+1],RO,0.0,RI8);

               DArray<9> tmp9;
               Contract(1.0,RI8,shape(1,3),peps(row+1,col+1),shape(1,4),0.0,tmp9);

               DArray<8> tmp8;
               Contract(1.0,tmp9,shape(1,7,2),peps(row+1,col+1),shape(1,2,4),0.0,tmp8);

               RI8.clear();
               Permute(tmp8,shape(0,4,6,5,7,1,2,3),RI8);

               //left
               Gemm(CblasTrans,CblasNoTrans,1.0,LO,env.gt(row)[col],0.0,LI8);

               tmp9.clear();
               Contract(1.0,LI8,shape(0,5),peps(row+1,col),shape(0,1),0.0,tmp9);

               tmp8.clear();
               Contract(1.0,tmp9,shape(0,4,6),peps(row+1,col),shape(0,1,2),0.0,tmp8);

               LI8.clear();
               Permute(tmp8,shape(3,5,7,4,6,0,1,2),LI8);

            }
            else{//col == Lx - 2

               //create left and right intermediary operators: right
               DArray<5> tmp5;
               Contract(1.0,env.gt(row)[Lx - 1],shape(2,3),peps(row+1,Lx - 1),shape(1,4),0.0,tmp5);

               DArray<6> tmp6;
               Contract(1.0,tmp5,shape(1,3),peps(row+1,Lx - 1),shape(1,2),0.0,tmp6);

               DArray<6> tmp6bis;
               Permute(tmp6,shape(0,3,1,4,2,5),tmp6bis);

               RI8 = tmp6bis.reshape_clear(shape(env.gt(row)[Lx-1].shape(0),D,D,D,D,1,1,1));

               //left
               Gemm(CblasTrans,CblasNoTrans,1.0,LO,env.gt(row)[col],0.0,LI8);

               DArray<9> tmp9;
               Contract(1.0,LI8,shape(0,5),peps(row+1,col),shape(0,1),0.0,tmp9);

               DArray<8> tmp8;
               Contract(1.0,tmp9,shape(0,4,6),peps(row+1,col),shape(0,1,2),0.0,tmp8);

               LI8.clear();
               Permute(tmp8,shape(3,5,7,4,6,0,1,2),LI8);

            }

         }
         else if(dir == DIAGONAL_LURD){

            if(col == 0){

               //create left and right intermediary operators: right

               //attach top environment to right
               Gemm(CblasNoTrans,CblasNoTrans,1.0,env.gt(row)[col+1],RO,0.0,RI8);

               //and upper peps
               DArray<9> tmp9;
               Contract(1.0,RI8,shape(1,3),peps(row+1,col+1),shape(1,4),0.0,tmp9);

               //again to create RI8
               DArray<8> tmp8;
               Contract(1.0,tmp9,shape(1,7,2),peps(row+1,col+1),shape(1,2,4),0.0,tmp8);

               RI8.clear();
               Permute(tmp8,shape(0,4,6,5,7,1,2,3),RI8);

               //b_R is equal to RI8, for lurd, so leave empty

               //left

               //attach bottom environment to bottom peps
               DArray<7> tmp7;
               Contract(1.0,env.gb(row-1)[col],shape(2),peps(row,col),shape(3),0.0,tmp7);

               //add another bottom peps for LI8 construction
               DArray<6> tmp6;
               Contract(1.0,tmp7,shape(0,5,1),peps(row,col),shape(0,2,3),0.0,tmp6);

               DArray<6> tmp6bis;
               Permute(tmp6,shape(1,4,2,5,3,0),tmp6bis);

               LI8 = tmp6bis.reshape_clear( shape(1,1,1,D,D,D,D,env.gb(row-1)[col].shape(3)) );

               //add mop to tmp9 for b_L construction
               Contract(1.0,tmp7,shape(0,5,1),mop,shape(0,2,3),0.0,tmp6);

               Permute(tmp6,shape(1,4,2,5,3,0),tmp6bis);

               b_L = tmp6bis.reshape_clear( shape(1,1,1,D,D,D,D,env.gb(row-1)[col].shape(3)) );

            }
            else if(col < Lx - 2){

               //create left and right intermediary operators: right

               //attach top environment to right
               Gemm(CblasNoTrans,CblasNoTrans,1.0,env.gt(row)[col+1],RO,0.0,RI8);

               //and upper peps
               DArray<9> tmp9;
               Contract(1.0,RI8,shape(1,3),peps(row+1,col+1),shape(1,4),0.0,tmp9);

               //again to create RI8
               DArray<8> tmp8;
               Contract(1.0,tmp9,shape(1,7,2),peps(row+1,col+1),shape(1,2,4),0.0,tmp8);

               RI8.clear();
               Permute(tmp8,shape(0,4,6,5,7,1,2,3),RI8);

               //b_R is equal to RI8, for lurd, so leave empty

               //left

               //attach bottom environment to LO
               Gemm(CblasNoTrans,CblasNoTrans,1.0,LO,env.gb(row-1)[col],0.0,LI8);

               //attach bottom-left peps
               tmp9.clear();
               Contract(1.0,LI8,shape(4,6),peps(row,col),shape(0,3),0.0,tmp9);

               //add another one for LI8 construction
               tmp8.clear();
               Contract(1.0,tmp9,shape(3,7,4),peps(row,col),shape(0,2,3),0.0,tmp8);

               LI8.clear();
               Permute(tmp8,shape(0,1,2,6,4,7,5,3),LI8);

               //attach mop to tmp9 for b_L construction
               Contract(1.0,tmp9,shape(3,7,4),mop,shape(0,2,3),0.0,tmp8);

               b_L.clear();
               Permute(tmp8,shape(0,1,2,6,4,7,5,3),b_L);

            }
            else{//col == Lx - 2

               //create left and right intermediary operators: right

               //attach top environment to upper-right peps
               DArray<5> tmp5;
               Contract(1.0,env.gt(row)[col+1],shape(2,3),peps(row+1,col+1),shape(1,4),0.0,tmp5);

               //and another one for RI8 construction
               DArray<6> tmp6;
               Contract(1.0,tmp5,shape(1,3),peps(row+1,col+1),shape(1,2),0.0,tmp6);

               DArray<6> tmp6bis;
               Permute(tmp6,shape(0,3,1,4,2,5),tmp6bis);

               RI8 = tmp6bis.reshape_clear( shape(env.gt(row)[col+1].shape(0),D,D,D,D,1,1,1 ) );

               //b_R is equal to RI8, for lurd, so leave empty

               //left

               //attach bottom environment to LO
               Gemm(CblasNoTrans,CblasNoTrans,1.0,LO,env.gb(row-1)[col],0.0,LI8);

               //attach bottom-left peps
               DArray<9> tmp9;
               Contract(1.0,LI8,shape(4,6),peps(row,col),shape(0,3),0.0,tmp9);

               //add another one for LI8 construction
               DArray<8> tmp8;
               Contract(1.0,tmp9,shape(3,7,4),peps(row,col),shape(0,2,3),0.0,tmp8);

               LI8.clear();
               Permute(tmp8,shape(0,1,2,6,4,7,5,3),LI8);

               //attach mop to tmp9 for b_L construction
               Contract(1.0,tmp9,shape(3,7,4),mop,shape(0,2,3),0.0,tmp8);

               b_L.clear();
               Permute(tmp8,shape(0,1,2,6,4,7,5,3),b_L);

            }

         }
         else{///dir == diagonal_ldru

            if(col == 0){

               //create left and right intermediary operators: right

               //add bottom environemtn to RO
               Gemm(CblasNoTrans,CblasTrans,1.0,env.gb(row-1)[col+1],RO,0.0,RI8);

               //attach bottom-right peps
               DArray<9> tmp9;
               Contract(1.0,peps(row,col+1),shape(3,4),RI8,shape(2,7),0.0,tmp9);

               //and another to construct RI8
               DArray<8> tmp8;
               Contract(1.0,peps(row,col+1),shape(2,3,4),tmp9,shape(2,4,8),0.0,tmp8);

               RI8.clear();
               Permute(tmp8,shape(5,6,7,1,3,0,2,4),RI8);

               //construct b_R by adding mop to tmp9
               Contract(1.0,mop,shape(2,3,4),tmp9,shape(2,4,8),0.0,tmp8);

               Permute(tmp8,shape(5,6,7,1,3,0,2,4),b_R);

               //left: add top env to upper left peps
               DArray<5> tmp5;
               Gemm(CblasTrans,CblasNoTrans,1.0,env.gt(row)[col],peps(row+1,col),0.0,tmp5);

               DArray<6> tmp6;
               Contract(1.0,tmp5,shape(0,2),peps(row+1,col),shape(1,2),0.0,tmp6);

               DArray<6> tmp6bis;
               Permute(tmp6,shape(3,0,2,5,1,4),tmp6bis);

               LI8 = tmp6bis.reshape_clear( shape(env.gt(row)[col].shape(3),D,D,D,D,1,1,1) );

               //b_L is equal to LI8 for ldru, so leave empty

            }
            else if(col < Lx - 2){

               //create left and right intermediary operators: right

               //add bottom environemtn to RO
               Gemm(CblasNoTrans,CblasTrans,1.0,env.gb(row-1)[col+1],RO,0.0,RI8);

               //attach bottom-right peps
               DArray<9> tmp9;
               Contract(1.0,peps(row,col+1),shape(3,4),RI8,shape(2,7),0.0,tmp9);

               //and another to construct RI8
               DArray<8> tmp8;
               Contract(1.0,peps(row,col+1),shape(2,3,4),tmp9,shape(2,4,8),0.0,tmp8);

               RI8.clear();
               Permute(tmp8,shape(5,6,7,1,3,0,2,4),RI8);

               //construct b_R by adding mop to tmp9
               Contract(1.0,mop,shape(2,3,4),tmp9,shape(2,4,8),0.0,tmp8);

               Permute(tmp8,shape(5,6,7,1,3,0,2,4),b_R);

               //left: add top env to LO
               Gemm(CblasTrans,CblasNoTrans,1.0,LO,env.gt(row)[col],0.0,LI8);

               //add upper left peps to intermediate
               tmp9.clear();
               Contract(1.0,LI8,shape(0,5),peps(row+1,col),shape(0,1),0.0,tmp9);

               //and another one to construct LI8
               tmp8.clear();
               Contract(1.0,tmp9,shape(0,4,6),peps(row+1,col),shape(0,1,2),0.0,tmp8);

               //permute to right order
               LI8.clear();
               Permute(tmp8,shape(3,5,7,4,6,0,1,2),LI8);

               //b_L is equal to LI8 for ldru, so leave empty


            }
            else{//col == Lx - 2

               //create left and right intermediary operators: right

               //add lower right peps to bottom environment
               DArray<5> tmp5;
               Gemm(CblasNoTrans,CblasTrans,1.0,peps(row,col+1),env.gb(row-1)[col+1],0.0,tmp5);

               //and another
               DArray<6> tmp6;
               Contract(1.0,peps(row,col+1),shape(2,3),tmp5,shape(2,4),0.0,tmp6);

               DArray<6> tmp6bis;
               Permute(tmp6,shape(2,1,4,0,3,5),tmp6bis);

               RI8 = tmp6bis.reshape_clear( shape(1,1,1,D,D,D,D,env.gb(row-1)[col+1].shape(0)) );

               //add mop to tmp5 for b_R construction
               Contract(1.0,mop,shape(2,3),tmp5,shape(2,4),0.0,tmp6);

               Permute(tmp6,shape(2,1,4,0,3,5),tmp6bis);

               b_R = tmp6bis.reshape_clear( shape(1,1,1,D,D,D,D,env.gb(row-1)[col+1].shape(0)) );

               //left: add top env to LO
               Gemm(CblasTrans,CblasNoTrans,1.0,LO,env.gt(row)[col],0.0,LI8);

               //add upper left peps to intermediate
               DArray<9> tmp9;
               Contract(1.0,LI8,shape(0,5),peps(row+1,col),shape(0,1),0.0,tmp9);

               //and another one to construct LI8
               DArray<8> tmp8;
               Contract(1.0,tmp9,shape(0,4,6),peps(row+1,col),shape(0,1,2),0.0,tmp8);

               //permute to right order
               LI8.clear();
               Permute(tmp8,shape(3,5,7,4,6,0,1,2),LI8);

               //b_L is equal to LI8 for ldru, so leave empty

            }

         }

      }

   /**
    * quasi canonicalize the environment of the sites to be updated, for stability reaons in the program.
    * @param dir vertical, horizontal,diagonal lurd or diagonal ldru update
    * @param row the row index of the bottom site
    * @param col column index of the vertical column
    * @param peps, full PEPS object before update
    * @param L Left environment contraction
    * @param R Right environment contraction
    * @param LI7 left intermediate object
    * @param RI7 right intermediate object
    */
   template<>
      void canonicalize(const PROP_DIR &dir,int row,int col,PEPS<double> &peps,const DArray<5> &L, const DArray<5> &R,
            
            const DArray<7> &LI7,const DArray<7> &RI7){

         if(dir == VERTICAL){

         }

      }

   /**
    * quasi canonicalize the environment of the sites to be updated, for stability reaons in the program.
    * @param dir vertical, horizontal,diagonal lurd or diagonal ldru update
    * @param row the row index of the bottom site
    * @param col column index of the vertical column
    * @param peps, full PEPS object before update
    * @param L Left environment contraction
    * @param R Right environment contraction
    * @param LI8 left intermediate object
    * @param RI8 right intermediate object
    */
   template<>
      void canonicalize(const PROP_DIR &dir,int row,int col,PEPS<double> &peps,const DArray<6> &L, const DArray<6> &R, const DArray<8> &LI8,const DArray<8> &RI8){

         if(dir == VERTICAL){

         }

      }

  /**
    * construct the single-site effective environment around a site to be updated
    * @param dir vertical, horizontal,diagonal lurd or diagonal ldru update
    * @param row the row index of the bottom site
    * @param col column index of the vertical column
    * @param peps, full PEPS object before update
    * @param N_eff output object, contains N_eff on output
    * @param L Left environment contraction
    * @param R Right environment contraction
    * @param LI7 left intermediate object
    * @param RI7 right intermediate object
    * @param left boolean flag for peps with left operator or right operator acting on it
    */
   template<>
      void calc_N_eff(const PROP_DIR &dir,int row,int col,PEPS<double> &peps, DArray<8> &N_eff,

            const DArray<5> &L, const DArray<5> &R, const DArray<7> &LI7,const DArray<7> &RI7, bool left){

         if(dir == VERTICAL){

            if(row == 0){

               if(col == 0){

                  if(left){//bottom site of vertical, so site (row,col) environment

                     //paste top peps to right intermediate
                     DArray<6> tmp6;
                     Contract(1.0,RI7,shape(0,2,4),peps(row+1,col),shape(0,1,4),0.0,tmp6);

                     //add another top peps for N_eff
                     DArray<5> tmp5;
                     Contract(1.0,tmp6,shape(0,4,1),peps(row+1,col),shape(1,2,4),0.0,tmp5);

                     DArray<5> tmp5bis;
                     Permute(tmp5,shape(3,4,0,2,1),tmp5bis);

                     N_eff = tmp5bis.reshape_clear( shape(1,D,1,D,1,D,1,D) );

                  }
                  else{//top site (row,col)

                     //paste bottom peps to right intermediate
                     DArray<10> tmp10;
                     Gemm(CblasNoTrans,CblasTrans,1.0,RI7,peps(row,col),0.0,tmp10);

                     //construct right hand side, attach operator to tmp10:
                     DArray<7> tmp7 = tmp10.reshape_clear( shape(D,D,D,D,D,D,d) );

                     //another bottom peps to this one
                     DArray<8> tmp8;
                     Contract(1.0,tmp7,shape(6,4),peps(row,col),shape(2,4),0.0,tmp8);

                     Permute(tmp8,shape(5,0,6,2,7,1,4,3),N_eff);

                  }

               }

            }

         }

      }

   /**
    * construct the single-site effective environment around a site to be updated
    * @param dir vertical, horizontal,diagonal lurd or diagonal ldru update
    * @param row the row index of the bottom site
    * @param col column index of the vertical column
    * @param peps, full PEPS object before update
    * @param N_eff output object, contains N_eff on output
    * @param LO Left environment contraction
    * @param RO Right environment contraction
    * @param LI8 left intermediate object
    * @param RI8 right intermediate object
    * @param left boolean flag for peps with left operator or right operator acting on it
    */
   template<>
      void calc_N_eff(const PROP_DIR &dir,int row,int col,PEPS<double> &peps, DArray<8> &N_eff,

            const DArray<6> &LO, const DArray<6> &RO, const DArray<8> &LI8,const DArray<8> &RI8, bool left){

      }

   /**
    * construct the single-site effective environment and right hand side needed for 
    * the linear system of any gate direction specified by 'dir'
    * @param dir vertical, horizontal,diagonal lurd or diagonal ldru update
    * @param row the row index of the bottom site
    * @param col column index of the vertical column
    * @param peps, full PEPS object before update
    * @param rhs output object, contains N_eff on output
    * @param L Left environment contraction
    * @param R Right environment contraction
    * @param LI7 left intermediate object
    * @param RI7 right intermediate object
    * @param b_L left intermediate object
    * @param b_R right intermediate object
    * @param left boolean flag for peps with left operator or right operator acting on it
    */
   template<>
      void calc_rhs(const PROP_DIR &dir,int row,int col,PEPS<double> &peps,const DArray<6> &lop,const DArray<6> &rop,

            DArray<5> &rhs, const DArray<5> &L, const DArray<5> &R, const DArray<7> &LI7,const DArray<7> &RI7,

            const DArray<7> &b_L,const DArray<7> &b_R,bool left){

      }

   /**
    * construct the single-site effective environment and right hand side needed for 
    * the linear system of any gate direction specified by 'dir'
    * @param dir vertical, horizontal,diagonal lurd or diagonal ldru update
    * @param row the row index of the bottom site
    * @param col column index of the vertical column
    * @param peps, full PEPS object before update
    * @param rhs output object, contains N_eff on output
    * @param L Left environment contraction
    * @param R Right environment contraction
    * @param LI7 left intermediate object
    * @param RI7 right intermediate object
    * @param b_L left intermediate object
    * @param b_R right intermediate object
    * @param left boolean flag for peps with left operator or right operator acting on it
    */
   template<>
      void calc_rhs(const PROP_DIR &dir,int row,int col,PEPS<double> &peps,const DArray<6> &lop,const DArray<6> &rop,

            DArray<5> &rhs, const DArray<6> &L, const DArray<6> &R, const DArray<8> &LI8,const DArray<8> &RI8,

            const DArray<8> &b_L,const DArray<8> &b_R,bool left){

      }
