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

         //construct for left and right operators acting the correct peps elements
         DArray<6> lop;
         DArray<6> rop;

         if(dir == VERTICAL){// (row,col) --> (row+1,col)

            Contract(1.0,peps(row,col),shape(i,j,k,l,m),global::trot.gLO_n(),shape(k,o,n),0.0,lop,shape(i,j,n,o,l,m));
            Contract(1.0,peps(row+1,col),shape(i,j,k,l,m),global::trot.gRO_n(),shape(k,o,n),0.0,rop,shape(i,j,n,o,l,m));

         }
         else if(dir == HORIZONTAL){// (row,col) --> (row,col+1)

            Contract(1.0,peps(row,col),shape(i,j,k,l,m),global::trot.gLO_n(),shape(k,o,n),0.0,lop,shape(i,j,n,o,l,m));
            Contract(1.0,peps(row,col+1),shape(i,j,k,l,m),global::trot.gRO_n(),shape(k,o,n),0.0,rop,shape(i,j,n,o,l,m));

         }
         else if(dir == DIAGONAL_LURD){//(row+1,col) --> (row,col+1)

            Contract(1.0,peps(row+1,col),shape(i,j,k,l,m),global::trot.gLO_nn(),shape(k,o,n),0.0,lop,shape(i,j,n,o,l,m));
            Contract(1.0,peps(row,col+1),shape(i,j,k,l,m),global::trot.gRO_nn(),shape(k,o,n),0.0,rop,shape(i,j,n,o,l,m));

         }
         else{//(row,col) --> (row+1,col+1)

            Contract(1.0,peps(row,col),shape(i,j,k,l,m),global::trot.gLO_nn(),shape(k,o,n),0.0,lop,shape(i,j,n,o,l,m));
            Contract(1.0,peps(row+1,col+1),shape(i,j,k,l,m),global::trot.gRO_nn(),shape(k,o,n),0.0,rop,shape(i,j,n,o,l,m));

         }

         // --- (a) --- initial guess:
         initialize(dir,row,col,lop,rop,peps); 

         // --- (b) --- create intermediary object using during N_eff construction, doesn't change during update
         construct_intermediate(dir,row,col,peps,L,R,LI,RI);

         // --- (c) --- sweeping update
         sweep(dir,row,col,peps,lop,rop,L,R,LI,RI,n_iter);

         // --- (d) --- set top and bottom back on equal footing
         equilibrate(dir,row,col,peps);

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
    * @param LI intermediate object created to simplify N_eff contstruction
    * @param RI intermediate object created to simplify N_eff contstruction
    * @param n_sweeps number of sweeps to execute
    */
   template<size_t M>
      void sweep(const PROP_DIR &dir,int row,int col,PEPS<double> &peps,const DArray<6> &lop,const DArray<6> &rop,

            const DArray<M> &L,const DArray<M> &R,const DArray<M+2> &LI,const DArray<M+2> &RI,int n_sweeps){

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

         while(iter < n_sweeps){

            // --(1)-- 'left' site

            //construct effective environment and right hand side for linear system of top site
            construct_lin_sys(dir,row,col,peps,lop,rop,N_eff,rhs,L,R,LI,RI,true);

            //solve the system
            solve(N_eff,rhs);

            //update 'left' peps
            Permute(rhs,shape(0,1,4,2,3),peps(l_row,l_col));

            // --(2)-- 'right' site

            //construct effective environment and right hand side for linear system of bottom site
            construct_lin_sys(dir,row,col,peps,lop,rop,N_eff,rhs,L,R,LI,RI,false);

            //solve the system
            solve(N_eff,rhs);

            //update 'right' peps
            Permute(rhs,shape(0,1,4,2,3),peps(r_row,r_col));

            //repeat until converged
            ++iter;

         }

      }



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
         update(VERTICAL,0,col,peps,L,R[col],n_sweeps); 

         // --- (2) update the horizontal pair on column 'col'-'col+1' ---
         update(HORIZONTAL,0,col,peps,L,R[col+1],n_sweeps); 

         // --- (3) update diagonal LU-RD
         //         update_diagonal_lurd(0,col,peps,L,R[col+1],n_sweeps); 

         // --- (4) update diagonal LD-RU
         // todo

         contractions::update_L('b',col,peps,L);

      }

      //one last vertical update
      update(VERTICAL,0,Lx-1,peps,L,R[Lx-2],n_sweeps); 

      //update the bottom row for the new peps
      env.gb(0).fill('b',peps);

      // ---------------------------------------------------//
      // --- !!! (2) the middle rows (1 -> Ly-2) (2) !!! ---// 
      // ---------------------------------------------------//

      //renormalized operators for the middle sites
      vector< DArray<6> > RO(Lx - 1);
      DArray<6> LO;

      for(int row = 1;row < Ly-2;++row){

         //first create right renormalized operator
         contractions::init_ro(row,peps,RO);

         for(int col = 0;col < Lx - 1;++col){

            // --- (1) update vertical pair on column 'col', with lowest site on row 'row'
            update(VERTICAL,row,col,peps,LO,RO[col],n_sweeps); 

            // --- (2) update the horizontal pair on column 'col'-'col+1' ---
            update(HORIZONTAL,row,col,peps,LO,RO[col+1],n_sweeps); 

            // --- (3) update diagonal LU-RD
            // todo

            // --- (4) update diagonal LD-RU
            // todo

            //first construct a double layer object for the newly updated bottom 
            contractions::update_L(row,col,peps,LO);

         }

         //finally, last vertical gate
         update(VERTICAL,row,Lx-1,peps,LO,RO[Lx-2],n_sweeps); 

         //finally update the 'bottom' environment for the row
         env.add_layer('b',row,peps);

      }

      // ----------------------------------------------------//
      // --- !!! (3) the top two rows (Ly-2,Ly-1) (3) !!! ---// 
      // ----------------------------------------------------//

      //make the right operators
      contractions::init_ro('t',peps,R);

      for(int col = 0;col < Lx - 1;++col){

         // --- (1) update vertical pair on column 'col' on upper two rows
         update(VERTICAL,Ly-2,col,peps,L,R[col],n_sweeps); 

         // --- (2a) update the horizontal pair on row 'row' and colums 'col'-'col+1' ---
         update(HORIZONTAL,Ly-2,col,peps,L,R[col+1],n_sweeps); 

         // --- (2b) update the horizontal pair on row 'row+1' and colums 'col'-'col+1' ---
         update(HORIZONTAL,Ly-1,col,peps,L,R[col+1],n_sweeps); 

         // --- (3) update diagonal LU-RD
         // todo

         // --- (4) update diagonal LD-RU
         // todo

         contractions::update_L('t',col,peps,L);

      }

      //finally the very last vertical gate
      update(VERTICAL,Ly-2,Lx-1,peps,L,R[Lx-2],n_sweeps); 

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
    * construct the single-site effective environment and right hand side needed for the linear system of any gate direction specified by 'dir'
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
    * @param left boolean flag for peps with left operator or right operator acting on it
    */
   template<>
      void construct_lin_sys(const PROP_DIR &dir,int row,int col,PEPS<double> &peps,const DArray<6> &lop,const DArray<6> &rop,DArray<8> &N_eff,DArray<5> &rhs,

            const DArray<5> &L, const DArray<5> &R, const DArray<7> &LI7,const DArray<7> &RI7,bool left){

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

                     //right hand side
                     DArray<6> tmp6bis;
                     Contract(1.0,tmp6,shape(0,4,1),rop,shape(1,2,5),0.0,tmp6bis);

                     DArray<4> tmp4;
                     Contract(1.0,tmp6bis,shape(3,5,4,0),lop,shape(0,1,3,5),0.0,tmp4);

                     DArray<4> tmp4bis;
                     Permute(tmp4,shape(1,3,0,2),tmp4bis);

                     rhs = tmp4bis.reshape_clear( shape(1,D,1,D,d) );
                  
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

                  //add left operator
                  tmp8bis.clear();
                  Contract(1.0,lop,shape(1,3,5),tmp8,shape(1,7,3),0.0,tmp8bis);

                  //and close with LI7
                  tmp5.clear();
                  Contract(1.0,LI7,shape(0,1,3,5,6),tmp8bis,shape(3,0,2,7,6),0.0,tmp5);

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

                  N_eff = tmp6bis.reshape_clear( shape(D,D,1,D,D,D,1,D) );

                  // (2) construct right hand side

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
            else{//bottom horizontal peps of topmost update

               if(left){//left site of horizontal gate, so site (row,col) environment

               }
               else{//right site of horizontal gate, so site (row+1,col) environment

               }

            }

         }//close DIAGONAL_LURD
         else{//DIAGONAL LDRU

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
    * @param left boolean flag for PEPS with left or right operator acted upon
    */
   template<>
      void construct_lin_sys(const PROP_DIR &dir,int row,int col,PEPS<double> &peps,const DArray<6> &lop,const DArray<6> &rop,DArray<8> &N_eff,DArray<5> &rhs,

            const DArray<6> &LO, const DArray<6> &RO, const DArray<8> &LI8,const DArray<8> &RI8,bool left){

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

         }//end diagonal lurd
         else{//diagonal lrdu

         }

      }

   /**
    * evaluate the cost function of the linear system for two vertically connected PEPS: -2 <\Psi|\Psi'> + <\Psi'|\Psi'> where \Psi full PEPS with operator and \Psi is old PEPS
    * for top or bottom rows, i.e. with L and R of order 5 and intermediates of order 7
    * @param dir vertical, horizontal,diagonal lurd or diagonal ldru update
    * @param row , the row index of the bottom site
    * @param col column index of the vertical column
    * @param peps, full PEPS object before update
    * @param L Left environment contraction
    * @param R Right environment contraction
    * @param LI7 left intermediate object
    * @param RI7 left intermediate object
    */
   template<>
      double cost_function(const PROP_DIR &dir,int row,int col,PEPS<double> &peps,const DArray<6> &lop,const DArray<6> &rop,

            const DArray<5> &L,const DArray<5> &R,const DArray<7> &LI7,const DArray<7> &RI7){

         if(dir == VERTICAL){

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
            else{//row == Lx - 2 

               if(col == 0){

                  // (1) construct N_eff

                  //paste bottom peps to right intermediate
                  DArray<8> tmp8;
                  Contract(1.0,peps(row,col),shape(3,4),RI7,shape(2,6),0.0,tmp8);

                  //and another bottom peps to tmp8
                  DArray<7> tmp7;
                  Contract(1.0,peps(row,col),shape(2,3,4),tmp8,shape(2,4,7),0.0,tmp7);

                  //upper peps
                  DArray<8> tmp8bis;
                  Contract(1.0,peps(row+1,col),shape(3,4),tmp7,shape(3,6),0.0,tmp8bis);

                  DArray<5> tmp5 = tmp8bis.reshape_clear( shape(1,1,d,D,D) );

                  double val = Dot(tmp5,peps(row+1,col));

                  //(2) right hand side:

                  //add left operator to tmp8
                  DArray<6> tmp6;
                  Contract(1.0,lop,shape(0,2,4,5),tmp8,shape(0,2,4,7),0.0,tmp6);

                  //add right operator
                  DArray<6> tmp6bis;
                  Contract(1.0,rop,shape(3,4,5),tmp6,shape(1,0,4),0.0,tmp6bis);

                  val -= 2.0 * blas::dot(tmp6bis.size(), tmp6bis.data(), 1, peps(row+1,col).data(), 1);

                  return val;

               }
               else if(col < Lx - 1){

                  // (1) construct N_eff

                  //paste bottom peps to right intermediate
                  DArray<8> tmp8;
                  Contract(1.0,peps(row,col),shape(3,4),RI7,shape(2,6),0.0,tmp8);

                  //and another bottom peps to tmp8
                  DArray<7> tmp7;
                  Contract(1.0,peps(row,col),shape(2,3,4),tmp8,shape(2,4,7),0.0,tmp7);

                  //add top peps
                  DArray<8> tmp8bis;
                  Contract(1.0,peps(row+1,col),shape(3,4),tmp7,shape(3,6),0.0,tmp8bis);

                  //final top peps
                  DArray<5> tmp5;
                  Contract(1.0,peps(row+1,col),shape(1,2,3,4),tmp8bis,shape(1,2,4,7),0.0,tmp5);

                  //contact with left hand side
                  double val = Dot(tmp5,L);

                  //(2) right hand side:

                  //add left operator to tmp8
                  tmp8bis.clear();
                  Contract(1.0,lop,shape(2,4,5),tmp8,shape(2,4,7),0.0,tmp8bis);

                  //add right operator
                  tmp8.clear();
                  Contract(1.0,rop,shape(3,4,5),tmp8bis,shape(2,1,6),0.0,tmp8);

                  //add top peps
                  tmp5.clear();
                  Contract(1.0,peps(row+1,col),shape(1,2,3,4),tmp8,shape(1,2,5,7),0.0,tmp5);

                  DArray<5> tmp5bis;
                  Permute(tmp5,shape(1,0,2,3,4),tmp5bis);

                  //and contract with left hand side
                  val -= 2.0 * Dot(tmp5bis,L);

                  return val;

               }
               else{//col == Lx - 1

                  // (1) construct N_eff

                  //paste bottom peps to left intermediate
                  DArray<6> tmp6;
                  Contract(1.0,peps(row,col),shape(0,3,4),LI7,shape(3,5,6),0.0,tmp6);

                  //and another
                  DArray<5> tmp5;
                  Contract(1.0,peps(row,col),shape(0,2,3),tmp6,shape(4,1,5),0.0,tmp5);

                  //add top peps to it
                  DArray<4> tmp4;
                  Contract(1.0,peps(row+1,col),shape(0,3,4),tmp5,shape(4,2,1),0.0,tmp4);

                  DArray<4> tmp4bis;
                  Permute(tmp4,shape(3,0,1,2),tmp4bis);

                  tmp5 = tmp4bis.reshape_clear( shape(D,1,d,D,1));

                  double val = Dot(tmp5,peps(row+1,col));

                  //(2) right hand side:

                  //add left operator to tmp6
                  DArray<6> tmp6bis;
                  Contract(1.0,lop,shape(0,2,4),tmp6,shape(4,1,5),0.0,tmp6bis);

                  //add right operator
                  tmp4.clear();
                  Contract(1.0,rop,shape(0,3,4,5),tmp6bis,shape(4,1,0,2),0.0,tmp4);

                  tmp4bis.clear();
                  Permute(tmp4,shape(3,0,1,2),tmp4bis);

                  tmp5 = tmp4bis.reshape_clear( shape(D,1,d,D,1));

                  val -=  2.0 * Dot(tmp5,peps(row+1,col));

                  return val;

               }

            }

         }//end vertical
         else if(dir == HORIZONTAL){

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
            else if(row == Ly - 2){

               //add left peps to LI7
               DArray<8> tmp8_left;
               Contract(1.0,LI7,shape(5,3),peps(row,col),shape(0,1),0.0,tmp8_left);

               //and another peps
               DArray<7> tmp7;
               Contract(1.0,tmp8_left,shape(3,2,5),peps(row,col),shape(0,1,2),0.0,tmp7);

               //add bottom environment
               DArray<5> tmp5_left;
               Contract(1.0,tmp7,shape(2,5,3),env.gb(Ly-3)[col],shape(0,1,2),0.0,tmp5_left);

               //add right peps to intermediate
               DArray<8> tmp8_right;
               Contract(1.0,RI7,shape(3,5),peps(row,col+1),shape(1,4),0.0,tmp8_right);

               //add second peps
               tmp7.clear();
               Contract(1.0,tmp8_right,shape(2,6,3),peps(row,col+1),shape(1,2,4),0.0,tmp7);

               //add bottom environment
               DArray<5> tmp5;
               Contract(1.0,tmp7,shape(6,4,2),env.gb(Ly-3)[col+1],shape(1,2,3),0.0,tmp5);

               double val = Dot(tmp5,tmp5_left);

               //right hand side

               //add left operator
               DArray<8> tmp8bis;
               Contract(1.0,tmp8_left,shape(3,2,5),lop,shape(0,1,2),0.0,tmp8bis);

               //add bottom environment
               DArray<6> tmp6_left;
               Contract(1.0,tmp8bis,shape(2,6,3),env.gb(Ly-3)[col],shape(0,1,2),0.0,tmp6_left);

               //add right operator to tmp8
               tmp8bis.clear();
               Contract(1.0,tmp8_right,shape(2,6,3),rop,shape(1,2,5),0.0,tmp8bis);

               //add bottom environment
               DArray<6> tmp6_right;
               Contract(1.0,env.gb(Ly-3)[col+1],shape(1,2,3),tmp8bis,shape(7,4,2),0.0,tmp6_right);

               DArray<6> tmp6bis;
               Permute(tmp6_right,shape(1,2,3,5,4,0),tmp6bis);

               val -= 2.0 * Dot(tmp6_left,tmp6bis);

               return val;

            }
            else{//row == Ly - 1

               //add right peps to intermediate
               DArray<8> tmp8_right;
               Contract(1.0,peps(row,col+1),shape(3,4),RI7,shape(3,1),0.0,tmp8_right);

               //add second peps
               DArray<5> tmp5_right;
               Contract(1.0,peps(row,col+1),shape(1,2,3,4),tmp8_right,shape(1,2,4,3),0.0,tmp5_right);

               //add left to intermediate
               DArray<8> tmp8_left;
               Contract(1.0,LI7,shape(1,3),peps(row,col),shape(0,3),0.0,tmp8_left);

               //and another
               DArray<5> tmp5_left;
               Contract(1.0,tmp8_left,shape(0,5,6,1),peps(row,col),shape(0,1,2,3),0.0,tmp5_left);

               DArray<5> tmp5bis;
               Permute(tmp5_left,shape(4,3,0,1,2),tmp5bis);

               //norm
               double val = Dot(tmp5bis,tmp5_right);

               //add right operator to tmp8_right
               DArray<6> tmp6_right;
               Contract(1.0,rop,shape(1,2,4,5),tmp8_right,shape(1,2,4,3),0.0,tmp6_right);

               //add left operator
               DArray<6> tmp6_left;
               Contract(1.0,tmp8_left,shape(0,5,6,1),lop,shape(0,1,2,4),0.0,tmp6_left);

               DArray<6> tmp6bis;
               Permute(tmp6_left,shape(5,4,3,0,1,2),tmp6bis);

               val -= 2.0 * Dot(tmp6_right,tmp6bis);

               return val;

            }

         }//end horizontal
         else if(dir == DIAGONAL_LURD){

            if(row == 0){

               // --- (1) calculate overlap of approximated state:

               //right side
               DArray<8> tmp8_right;
               Contract(1.0,RI7,shape(4,6),peps(row,col+1),shape(1,4),0.0,tmp8_right);

               //and again
               DArray<5> tmp5_right;
               Contract(1.0,tmp8_right,shape(3,6,7,4),peps(row,col+1),shape(1,2,3,4),0.0,tmp5_right);

               //left side
               DArray<8> tmp8_left;
               Contract(1.0,LI7,shape(2,4),peps(row+1,col),shape(0,3),0.0,tmp8_left);

               //add another
               DArray<7> tmp7;
               Contract(1.0,tmp8_left,shape(1,6,2),peps(row+1,col),shape(0,2,3),0.0,tmp7);

               //add top environment to intermediate
               DArray<5> tmp5_left;
               Contract(1.0,tmp7,shape(0,5,3),env.gt(row)[col],shape(0,1,2),0.0,tmp5_left);

               DArray<5> tmp5;
               Permute(tmp5_left,shape(4,3,2,1,0),tmp5);

               double val = Dot(tmp5,tmp5_right);

               // --- (2) calculate 'b' part of overlap

               //add right operator to tmp8
               DArray<6> tmp6_right;
               Contract(1.0,tmp8_right,shape(3,6,7,4),rop,shape(1,2,4,5),0.0,tmp6_right);

               //add left operator to tmp8
               DArray<8> tmp8bis;
               Contract(1.0,tmp8_left,shape(1,6,2),lop,shape(0,2,4),0.0,tmp8bis);

               //add top environment to intermediate
               DArray<6> tmp6_left;
               Contract(1.0,tmp8bis,shape(0,5,3),env.gt(row)[col],shape(0,1,2),0.0,tmp6_left);

               DArray<6> tmp6bis;
               Permute(tmp6_left,shape(5,4,2,1,0,3),tmp6bis);

               val -= 2.0 * Dot(tmp6bis,tmp6_right);

               return val;

            }
            else{//row == Ly - 2

               return 0.0;

            }


         }
         else{//diagonal ldru

            return 0.0;

         }

      }

   /**
    * evaluate the cost function of the linear system for two vertically connected PEPS: -2 <\Psi|\Psi'> + <\Psi'|\Psi'> where \Psi full PEPS with operator and \Psi is old PEPS
    * for middle rows, i.e. with R and L order 6 and intermediates of order 8
    * @param dir vertical, horizontal,diagonal lurd or diagonal ldru update
    * @param row , the row index of the bottom site
    * @param col column index of the vertical column
    * @param peps, full PEPS object before update
    * @param LO Left environment contraction
    * @param RO Right environment contraction
    * @param LI8 left intermediate object
    * @param RI8 right intermediate object
    */
   template<>
      double cost_function(const PROP_DIR &dir,int row,int col,PEPS<double> &peps,const DArray<6> &lop,const DArray<6> &rop,

            const DArray<6> &LO,const DArray<6> &RO,const DArray<8> &LI8,const DArray<8> &RI8){

         if(dir == VERTICAL){

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
            else if(col < Lx -  1){//col != 0

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
            else{//col = Lx - 1

               //add bottom peps  to intermediate
               DArray<9> tmp9;
               Contract(1.0,LI8,shape(3,5),peps(row,col),shape(0,3),0.0,tmp9);

               //and another
               DArray<8> tmp8;
               Contract(1.0,tmp9,shape(2,7,3),peps(row,col),shape(0,2,3),0.0,tmp8);

               //upper peps
               DArray<5> tmp5;
               Contract(1.0,tmp8,shape(0,2,6,7),peps(row+1,col),shape(0,1,3,4),0.0,tmp5);

               DArray<5> tmp5bis;
               Permute(tmp5,shape(0,1,4,2,3),tmp5bis);

               double val = Dot(tmp5bis,peps(row+1,col));

               //add left operator to intermediate
               DArray<9> tmp9bis;
               Contract(1.0,tmp9,shape(2,7,3),lop,shape(0,2,4),0.0,tmp9bis);

               //and right operator
               tmp5.clear();
               Contract(1.0,tmp9bis,shape(0,2,7,6,8),rop,shape(0,1,3,4,5),0.0,tmp5);

               tmp5bis.clear();;
               Permute(tmp5,shape(0,1,4,2,3),tmp5bis);

               val -= 2.0 * Dot(tmp5bis,peps(row+1,col));

               return val;

            }

         }//end vertical
         else if(dir == HORIZONTAL){

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

         }//end horizontal
         else if(dir == DIAGONAL_LURD){

            return 0.0;

         }
         else{//diagonal ldru

            return 0.0;

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

      }
      else{//diagonal LDRU

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

      }
      else{//diagonal ldru

      }

   }

   /**
    * function that calculates intermediate objects that do not change during the sweeping update.
    * By precalculating them a lot of work is avoided
    * @param dir vertical, horizontal,diagonal lurd or diagonal ldru update
    * @param row row index of the bottom peps of the vertical pair
    * @param col column index
    * @param peps full PEPS object
    * @param L left contracted environment around the pair
    * @param R right contracted environment around the pair
    * @param LI7 Left intermediate object to be constructed on output
    * @param RI7 Right intermediate object to be constructed on output
    */
   template<>
      void construct_intermediate(const PROP_DIR &dir,int row,int col,const PEPS<double> &peps,const DArray<5> &L,const DArray<5> &R,DArray<7> &LI7,DArray<7> &RI7){

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

                  Permute(tmp7,shape(2,0,3,1,4,5,6),RI7);

                  //left

                  //add top peps to left
                  tmp8.clear();
                  Contract(1.0,L,shape(0),peps(row+1,col),shape(0),0.0,tmp8);

                  tmp7.clear();
                  Contract(1.0,tmp8,shape(0,4,5),peps(row+1,col),shape(0,1,2),0.0,tmp7);

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

                  //left: connect bottom left peps to itsself
                  DArray<4> tmp4;
                  Contract(1.0,peps(row,col),shape(0,2,3),peps(row,col),shape(0,2,3),0.0,tmp4);

                  DArray<4> tmp4bis;
                  Permute(tmp4,shape(0,2,1,3),tmp4bis);

                  LI7 = tmp4bis.reshape_clear( shape(1,1,1,D,D,D,D) );

               }
               else if(col < Lx - 2){

               }
               else{//col == Lx - 2

               }

            }
            else{// row == Ly - 2)

               if(col == 0){

               }
               else if(col < Lx - 2){

               }
               else{//col == Lx - 2

               }

            }

         }
         else{//DIAGONAL LDRU

         }

      }

   /**
    * function that calculates intermediate objects that do not change during the sweeping update.
    * By precalculating them a lot of work is avoided
    * @param dir vertical, horizontal,diagonal lurd or diagonal ldru update
    * @param row row index of the bottom peps of the vertical pair
    * @param col column index
    * @param peps full PEPS object
    * @param LO left contracted environment around the pair
    * @param RO right contracted environment around the pair
    * @param LI8 Left intermediate object to be constructed on output
    * @param RI8 Right intermediate object to be constructed on output
    */
   template<>
      void construct_intermediate(const PROP_DIR &dir,int row,int col,const PEPS<double> &peps,const DArray<6> &LO,const DArray<6> &RO,DArray<8> &LI8,DArray<8> &RI8){

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

         }
         else{///dir == diagonal_ldru

         }

      }

}
