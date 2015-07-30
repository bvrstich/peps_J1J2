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
      void update(const PROP_DIR &dir,int row,int col,PEPS<double> &peps,DArray<M> &L,DArray<M> &R,int n_iter){

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

         // --- (a) --- first construct some intermediates for all the environment calculations
         construct_intermediate(dir,row,col,peps,L,R,LI,RI);

         //the environment 'R' matrices from the QR decompositions of the environment
         std::vector< DArray<2> > R_l(4);
         std::vector< DArray<2> > R_r(4);

         // --- (b) --- canonicalize the environments around the sites to be updated
         canonicalize(dir,row,col,peps,L,R,LI,RI,R_l,R_r);

         if(dir == VERTICAL){// (row,col) --> (row+1,col)

            //left and right operators:
            Contract(1.0,peps(row,col),shape(i,j,k,l,m),global::trot.gLO_n(),shape(k,o,n),0.0,lop,shape(i,j,n,o,l,m));
            Contract(1.0,peps(row+1,col),shape(i,j,k,l,m),global::trot.gRO_n(),shape(k,o,n),0.0,rop,shape(i,j,n,o,l,m));

         }
         else if(dir == HORIZONTAL){// (row,col) --> (row,col+1)

            //left and right operators:
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
/*
         // --- (c) --- initial guess: use SVD to initialize the tensors
         initialize(dir,row,col,lop,rop,peps); 

         //recalculate the intermediates for diagonals
         if(dir == DIAGONAL_LURD || dir == DIAGONAL_LDRU){

            construct_intermediate(dir,row,col,peps,L,R,LI,RI);

            construct_intermediate_rhs(dir,row,col,peps,mop,L,R,b_L,b_R);

         }

#ifdef _DEBUG
         DArray<8> N_eff;
         calc_N_eff(dir,row,col,peps,N_eff,L,R,LI,RI,true);

         DArray<1> eig;
         diagonalize(N_eff,eig);

         cout << endl;
         cout << "after svd" << endl;
         cout << endl;

         cout << "(left)\t" << std::scientific << eig(eig.size() - 1) / eig(0) << endl;

         N_eff.clear();
         calc_N_eff(dir,row,col,peps,N_eff,L,R,LI,RI,false);

         eig.clear();
         diagonalize(N_eff,eig);

         cout << "(right)\t" << std::scientific << eig(eig.size() - 1) / eig(0) << endl;
         cout << endl;
#endif

         // --- (d) --- sweeping update: ALS
         sweep(dir,row,col,peps,lop,rop,L,R,LI,RI,b_L,b_R,n_iter);
*/
         // --- (e) --- restore the tensors, i.e. undo the canonicalization
         restore(dir,row,col,peps,L,R,R_l,R_r);

         // --- (f) --- set top and bottom back on equal footing
//         equilibrate(dir,row,col,peps);

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

         while(iter < n_sweeps){

#ifdef _DEBUG
            cout << iter << "\t" << debug::cost_function(dir,row,col,peps,lop,rop,L,R,LI,RI,b_L,b_R) << endl;
#endif

            // --(1)-- 'left' site

            //construct effective environment and right hand side for linear system of top site
            calc_N_eff(dir,row,col,peps,N_eff,L,R,LI,RI,true);
            regularize(N_eff,reg_const);

            calc_rhs(dir,row,col,peps,lop,rop,rhs,L,R,LI,RI,b_L,b_R,true);

            //solve the system
            solve(N_eff,rhs);

            //update 'left' peps
            Permute(rhs,shape(0,1,4,2,3),peps(l_row,l_col));

#ifdef _DEBUG
            cout << iter << "\t" << debug::cost_function(dir,row,col,peps,lop,rop,L,R,LI,RI,b_L,b_R) << endl;
#endif

            // --(2)-- 'right' site

            //construct effective environment and right hand side for linear system of bottom site
            calc_N_eff(dir,row,col,peps,N_eff,L,R,LI,RI,false);
            regularize(N_eff,reg_const);

            calc_rhs(dir,row,col,peps,lop,rop,rhs,L,R,LI,RI,b_L,b_R,false);

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

      //'canonicalize' top environment
      for(int row = Ly - 1;row > 1;--row)
         shift_row('t',row,peps);

      peps.rescale_tensors(scal_num);

      //calculate top environment
      env.calc('T',peps);

      //containers for the renormalized operators
      vector< DArray<5> > R(Lx);

      //'canonicalize' the right peps for the bottom two rows
      for(int col = Lx - 1;col > 0;--col){

         shift_col('l',0,col,peps);
         shift_col('l',1,col,peps);

      }

      //initialize the right operators for the bottom row
      contractions::init_ro('b',peps,R);

      contractions::rescale_norm('b',peps,R);

      //row == 0
      DArray<5> L(1,1,1,1,1);
      L = 1.0;

      for(int col = 0;col < Lx - 1;++col){

#ifdef _DEBUG
         cout << endl;
         cout << "***************************************" << endl;
         cout << " Update site:\t(0," << col << ")" << endl;
         cout << "***************************************" << endl;
         cout << endl;
#endif

         // --- (1) update the vertical pair on column 'col' ---
         update(VERTICAL,0,col,peps,L,R[col],n_sweeps); 

         // --- (2) update the horizontal pair on column 'col'-'col+1' ---
         update(HORIZONTAL,0,col,peps,L,R[col+1],n_sweeps); 

         // --- (3) update diagonal LU-RD
         update(DIAGONAL_LURD,0,col,peps,L,R[col+1],n_sweeps); 

         // --- (4) update diagonal LD-RU
         update(DIAGONAL_LDRU,0,col,peps,L,R[col+1],n_sweeps); 

         //do a QR decomposition of the updated peps on 'col'
         shift_col('r',0,col,peps);
         shift_col('r',1,col,peps);

         contractions::update_L('b',col,peps,L);

      }

      //one last vertical update
      update(VERTICAL,0,Lx-1,peps,L,R[Lx-1],n_sweeps); 

      //QR the complete row
      shift_row('b',0,peps);
      peps.rescale_tensors(0,scal_num);

      //and make the new bottom environment
      env.gb(0).fill('b',peps);

      //all middle rows:
      //for(int row = 1;row < Lx - 2;++row){
      int row = 1;

      //containers for the renormalized operators
      vector< DArray<6> > RO(Lx);

      //'canonicalize' the right peps for two rows to be updates
      for(int col = Lx - 1;col > 0;--col){

      shift_col('l',row,col,peps);
      shift_col('l',row+1,col,peps);

      }

      //initialize the right operators for the bottom row
      contractions::init_ro(row,peps,RO);

      contractions::rescale_norm(row,peps,RO);

      DArray<6> LO(1,1,1,1,1,1);
      LO = 1.0;

      //for(int col = 0;col < Lx - 1;++col){
      int col = 0;

#ifdef _DEBUG
      cout << endl;
      cout << "***************************************" << endl;
      cout << " Update site:\t(" << row << "," << col << ")" << endl;
      cout << "***************************************" << endl;
      cout << endl;
#endif

      // --- (1) update the vertical pair on column 'col' ---
      update(VERTICAL,row,col,peps,LO,RO[col],n_sweeps); 

      // --- (2) update the horizontal pair on column 'col'-'col+1' ---
      update(HORIZONTAL,row,col,peps,LO,RO[col+1],n_sweeps); 

      // --- (3) update diagonal LU-RD
      update(DIAGONAL_LURD,row,col,peps,LO,RO[col+1],n_sweeps); 

      // --- (4) update diagonal LD-RU
      update(DIAGONAL_LDRU,row,col,peps,LO,RO[col+1],n_sweeps); 

      //do a QR decomposition of the updated peps on 'col'
      shift_col('r',row,col,peps);
      shift_col('r',row+1,col,peps);

      contractions::update_L(row,col,peps,LO);

      //}
      /*
      //one last vertical update
      update(VERTICAL,row,Lx-1,peps,LO,RO[Lx-1],n_sweeps); 

      //QR the complete row
      shift_row('b',row,peps);
      peps.rescale_tensors(row,scal_num);

      //update the environment
      env.add_layer('b',row,peps);

      //}

      // finally row == Lx-2

      //'canonicalize' the right peps for the top two rows
      for(int col = Lx - 1;col > 0;--col){

      shift_col('l',Ly-2,col,peps);
      shift_col('l',Ly-1,col,peps);

      }


      //initialize the right operators for the top two rows
      contractions::init_ro('t',peps,R);

      contractions::rescale_norm('t',peps,R);

      L.resize(shape(1,1,1,1,1));
      L = 1.0;

      for(int col = 0;col < Lx - 1;++col){

#ifdef _DEBUG
cout << endl;
cout << "***************************************" << endl;
cout << " Update site:\t(" << Lx-2 << "," << col << ")" << endl;
cout << "***************************************" << endl;
cout << endl;
#endif

      // --- (1) update the vertical pair on column 'col' ---
      update(VERTICAL,Ly-2,col,peps,L,R[col],n_sweeps); 

      // --- (2a) update the horizontal pair on row Ly-2 column 'col'-'col+1' ---
      update(HORIZONTAL,Ly-2,col,peps,L,R[col+1],n_sweeps); 

      // --- (2b) update the horizontal pair on row Ly-1 column 'col'-'col+1' ---
      update(HORIZONTAL,Ly-1,col,peps,L,R[col+1],n_sweeps); 

      // --- (3) update diagonal LU-RD
      //update(DIAGONAL_LURD,0,col,peps,L,R[col+1],n_sweeps); 

      // --- (4) update diagonal LD-RU
      //update(DIAGONAL_LDRU,0,col,peps,L,R[col+1],n_sweeps); 

      //do a QR decomposition of the updated peps on 'col'
      shift_col('r',Ly-2,col,peps);
      shift_col('r',Ly-1,col,peps);

      contractions::update_L('t',col,peps,L);

      }

      //one last vertical update
      update(VERTICAL,Ly-2,Lx-1,peps,L,R[Lx-1],n_sweeps); 
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

         if(dir == DIAGONAL_LURD){

            if(row == 0){

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

         }//end VERTICAL
         else if(dir == HORIZONTAL){

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
    * @param L left contracted environment around the pair
    * @param R right contracted environment around the pair
    * @param LI7 Left intermediate object to be constructed on output
    * @param RI7 Right intermediate object to be constructed on output
    */
   template<>
      void construct_intermediate(const PROP_DIR &dir,int row,int col,const PEPS<double> &peps,

            const DArray<5> &L,const DArray<5> &R,DArray<7> &LI7,DArray<7> &RI7){

         if(dir == VERTICAL){

            //only LI7
            if(row == 0){

               DArray<7> tmp7;
               Gemm(CblasTrans,CblasNoTrans,1.0,L,env.gt(0)[col],0.0,tmp7);

               Permute(tmp7,shape(6,4,5,0,1,2,3),LI7);

            }
            else{//row == Lx - 2

               Gemm(CblasNoTrans,CblasNoTrans,1.0,L,env.gb(Lx - 3)[col],0.0,LI7);

            }

         }
         else if(dir == HORIZONTAL){

            if(row == 0){

               //right
               RI7.clear();
               Gemm(CblasNoTrans,CblasNoTrans,1.0,env.gt(0)[col+1],R,0.0,RI7);

               DArray<8> tmp8;
               Contract(1.0,RI7,shape(1,3),peps(row+1,col+1),shape(1,4),0.0,tmp8);

               DArray<7> tmp7;
               Contract(1.0,tmp8,shape(1,6,2),peps(row+1,col+1),shape(1,2,4),0.0,tmp7);

               RI7.clear();
               Permute(tmp7,shape(0,3,5,4,6,1,2),RI7);

               //left
               LI7.clear();
               Gemm(CblasTrans,CblasNoTrans,1.0,L,env.gt(0)[col],0.0,LI7);

               tmp8.clear();
               Contract(1.0,LI7,shape(0,4),peps(row+1,col),shape(0,1),0.0,tmp8);

               tmp7.clear();
               Contract(1.0,tmp8,shape(0,3,5),peps(row+1,col),shape(0,1,2),0.0,tmp7);

               LI7.clear();
               Permute(tmp7,shape(2,4,6,3,5,0,1),LI7);

            }
            else if(row == Ly - 2){

               //Left
               Gemm(CblasNoTrans,CblasNoTrans,1.0,L,env.gb(Ly - 3)[col],0.0,LI7);

               //Right
               DArray<7> tmp7;
               Gemm(CblasNoTrans,CblasTrans,1.0,env.gb(Ly - 3)[col+1],R,0.0,tmp7);

               RI7.clear();
               Permute(tmp7,shape(3,4,5,6,1,2,0),RI7);

            }
            else{//row == Ly -1 

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

         }
         else if(dir == DIAGONAL_LURD){

            if(row == 0){

               //first left intermediate:

               //add lower left peps to L
               DArray<8> tmp8;
               Contract(1.0,L,shape(4),peps(row,col),shape(0),0.0,tmp8);

               //and again
               DArray<7> tmp7;
               Contract(1.0,tmp8,shape(3,5,6),peps(row,col),shape(0,2,3),0.0,tmp7);

               LI7.clear();
               Permute(tmp7,shape(0,1,2,5,3,6,4),LI7);

               //right
               RI7.clear();
               Gemm(CblasNoTrans,CblasNoTrans,1.0,env.gt(0)[col+1],R,0.0,RI7);

               //add top right peps
               tmp8.clear();
               Contract(1.0,peps(row+1,col+1),shape(1,4),RI7,shape(1,3),0.0,tmp8);

               //and again
               tmp7.clear();
               Contract(1.0,peps(row+1,col+1),shape(1,2,4),tmp8,shape(4,1,5),0.0,tmp7);

               RI7.clear();
               Permute(tmp7,shape(4,2,0,3,1,5,6),RI7);

            }
            else{// row == Ly - 2

            }

         }
         else{//DIAGONAL LDRU

            if(row == 0){

               //first left intermediate:

               //add top environment to L
               DArray<7> tmp7;
               Contract(1.0,L,shape(0),env.gt(0)[col],shape(0),0.0,tmp7);

               //add left upper peps to tmp7
               DArray<8> tmp8;
               Contract(1.0,tmp7,shape(0,4),peps(row+1,col),shape(0,1),0.0,tmp8);

               //and again
               tmp7.clear();
               Contract(1.0,tmp8,shape(0,3,5),peps(row+1,col),shape(0,1,2),0.0,tmp7);

               LI7.clear();
               Permute(tmp7,shape(2,4,6,3,5,0,1),LI7);

               //right, add right lower peps to R
               tmp8.clear();
               Contract(1.0,peps(row,col+1),shape(4),R,shape(4),0.0,tmp8);

               //and again
               tmp7.clear();
               Contract(1.0,peps(row,col+1),shape(2,3,4),tmp8,shape(2,3,7),0.0,tmp7);

               RI7.clear();
               Permute(tmp7,shape(4,5,6,1,3,0,2),RI7);

            }
            else{// row == Ly - 2

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
    * @param L left contracted environment around the pair
    * @param R right contracted environment around the pair
    * @param b_L Left intermediate object to be constructed on output
    * @param b_R Right intermediate object to be constructed on output
    */
   template<>
      void construct_intermediate_rhs(const PROP_DIR &dir,int row,int col,const PEPS<double> &peps,

            const DArray<5> &mop, const DArray<5> &L,const DArray<5> &R,DArray<7> &b_L,DArray<7> &b_R){

         if(dir == DIAGONAL_LURD){

            if(row == 0){//only b_L here

               //add lower left peps to L
               DArray<8> tmp8;
               Contract(1.0,L,shape(4),peps(row,col),shape(0),0.0,tmp8);

               //add mop on top
               DArray<7> tmp7;
               Contract(1.0,tmp8,shape(3,5,6),mop,shape(0,2,3),0.0,tmp7);

               b_L.clear();
               Permute(tmp7,shape(0,1,2,5,3,6,4),b_L);

            }
            else{// row == Ly - 2

            }

         }
         else{//DIAGONAL LDRU

            if(row == 0){//only b_R

               //right, add right lower peps to R
               DArray<8> tmp8;
               Contract(1.0,peps(row,col+1),shape(4),R,shape(4),0.0,tmp8);

               //and add mop to top
               DArray<7> tmp7;
               Contract(1.0,mop,shape(2,3,4),tmp8,shape(2,3,7),0.0,tmp7);

               b_R.clear();
               Permute(tmp7,shape(4,5,6,1,3,0,2),b_R);

            }
            else{// row == Ly - 2

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
    * @param LO left contracted environment around the pair
    * @param RO right contracted environment around the pair
    * @param LI8 Left intermediate object to be constructed on output
    * @param RI8 Right intermediate object to be constructed on output
    */
   template<>
      void construct_intermediate(const PROP_DIR &dir,int row,int col,const PEPS<double> &peps,

            const DArray<6> &LO,const DArray<6> &RO,DArray<8> &LI8,DArray<8> &RI8){

         if(dir == VERTICAL){

            //right
            DArray<8> tmp8;
            Gemm(CblasNoTrans,CblasTrans,1.0,env.gb(row-1)[col],RO,0.0,tmp8);

            Permute(tmp8,shape(3,4,5,6,7,1,2,0),RI8);

            //left
            tmp8.clear();
            Gemm(CblasTrans,CblasNoTrans,1.0,LO,env.gt(row)[col],0.0,tmp8);

            Permute(tmp8,shape(7,5,6,0,1,2,3,4),LI8);

         }
         else if(dir == HORIZONTAL){

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
         else if(dir == DIAGONAL_LURD){

            //create left and right intermediary operators: right
            RI8.clear();
            Gemm(CblasNoTrans,CblasNoTrans,1.0,env.gt(row)[col+1],RO,0.0,RI8);

            DArray<9> tmp9;
            Contract(1.0,RI8,shape(1,3),peps(row+1,col+1),shape(1,4),0.0,tmp9);

            DArray<8> tmp8;
            Contract(1.0,tmp9,shape(1,7,2),peps(row+1,col+1),shape(1,2,4),0.0,tmp8);

            RI8.clear();
            Permute(tmp8,shape(0,4,6,5,7,1,2,3),RI8);

            //left
            LI8.clear();
            Gemm(CblasNoTrans,CblasNoTrans,1.0,LO,env.gb(row-1)[col],0.0,LI8);

            tmp9.clear();
            Contract(1.0,LI8,shape(3,5),peps(row,col),shape(0,3),0.0,tmp9);

            tmp8.clear();
            Contract(1.0,tmp9,shape(3,7,4),peps(row,col),shape(0,2,3),0.0,tmp8);

            LI8.clear();
            Permute(tmp8,shape(0,1,2,4,6,5,7,3),LI8);

         }
         else{///dir == diagonal_ldru

            //create left and right intermediary operators: right
            RI8.clear();
            Gemm(CblasNoTrans,CblasTrans,1.0,env.gb(row-1)[col+1],RO,0.0,RI8);

            DArray<9> tmp9;
            Contract(1.0,peps(row,col+1),shape(3,4),RI8,shape(2,7),0.0,tmp9);

            DArray<8> tmp8;
            Contract(1.0,peps(row,col+1),shape(2,3,4),tmp9,shape(2,4,8),0.0,tmp8);

            RI8.clear();
            Permute(tmp8,shape(5,6,7,1,3,0,2,4),RI8);

            //left
            Gemm(CblasTrans,CblasNoTrans,1.0,LO,env.gt(row)[col],0.0,LI8);

            tmp9.clear();
            Contract(1.0,LI8,shape(0,5),peps(row+1,col),shape(0,1),0.0,tmp9);

            tmp8.clear();
            Contract(1.0,tmp9,shape(0,4,6),peps(row+1,col),shape(0,1,2),0.0,tmp8);

            LI8.clear();
            Permute(tmp8,shape(3,5,7,4,6,0,1,2),LI8);

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

               if(left){//bottom site

                  //paste top peps to left
                  DArray<8> tmp8;
                  Contract(1.0,LI7,shape(3,1),peps(row+1,col),shape(0,1),0.0,tmp8);

                  //and another: watch out, order is reversed!
                  DArray<7> tmp7;
                  Contract(1.0,tmp8,shape(2,1,5),peps(row+1,col),shape(0,1,2),0.0,tmp7);

                  //now add right side to it
                  DArray<6> tmp6;
                  Contract(1.0,tmp7,shape(0,4,6),R,shape(0,1,2),0.0,tmp6);

                  DArray<6> tmp6bis;
                  Permute(tmp6,shape(0,2,4,1,3,5),tmp6bis);

                  int DL = peps(row,col).shape(0);
                  int DR = peps(row,col).shape(4);

                  N_eff = tmp6bis.reshape_clear(shape(DL,D,1,DR,DL,D,1,DR));

               }
               else{//top site

                  DArray<8> tmp8;
                  Contract(1.0,peps(row,col),shape(4),R,shape(4),0.0,tmp8);

                  DArray<7> tmp7;
                  Contract(1.0,peps(row,col),shape(2,3,4),tmp8,shape(2,3,7),0.0,tmp7);

                  //contract left and right
                  tmp8.clear();
                  Contract(1.0,LI7,shape(0,5,6),tmp7,shape(4,0,2),0.0,tmp8);

                  N_eff.clear();
                  Permute(tmp8,shape(2,0,4,6,3,1,5,7),N_eff);

               }

            }
            else{//row == Lx-2

               if(left){//bottom site

                  //paste top peps to right
                  DArray<8> tmp8;
                  Contract(1.0,peps(row+1,col),shape(4),R,shape(0),0.0,tmp8);

                  //and again
                  DArray<7> tmp7;
                  Contract(1.0,peps(row+1,col),shape(1,2,4),tmp8,shape(1,2,4),0.0,tmp7);

                  //attach left to right
                  tmp8.clear();
                  Contract(1.0,LI7,shape(0,1,6),tmp7,shape(2,0,6),0.0,tmp8);

                  N_eff.clear();
                  Permute(tmp8,shape(0,5,2,6,1,4,3,7),N_eff);

               }
               else{//top site

                  //paste top peps to right
                  DArray<8> tmp8;
                  Contract(1.0,LI7,shape(3,5),peps(row,col),shape(0,3),0.0,tmp8);

                  //and again
                  DArray<7> tmp7;
                  Contract(1.0,tmp8,shape(2,6,3),peps(row,col),shape(0,2,3),0.0,tmp7);

                  //attach left to right
                  DArray<6> tmp6;
                  Contract(1.0,tmp7,shape(6,4,2),R,shape(2,3,4),0.0,tmp6);

                  DArray<6> tmp6bis;
                  Permute(tmp6,shape(0,3,4,1,2,5),tmp6bis);

                  int DL = peps(row+1,col).shape(0);
                  int DR = peps(row+1,col).shape(4);

                  N_eff = tmp6bis.reshape_clear( shape(DL,1,D,DR,DL,1,D,DR) );

               }

            }

         }
         else if(dir == HORIZONTAL){

            if(row == 0){

               if(left){//left site of horizontal gate, so site (row,col) environment

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

               }
               else{//right site of horizontal gate, so site (row+1,col) environment

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

               }

            }
            else if(row == Ly - 2){

               if(left){

                  //fill up the right side
                  DArray<8> tmp8;
                  Contract(1.0,peps(row,col+1),shape(3,4),RI7,shape(5,3),0.0,tmp8);

                  DArray<7> tmp7;
                  Contract(1.0,peps(row,col+1),shape(2,3,4),tmp8,shape(2,6,5),0.0,tmp7);

                  tmp8.clear();
                  Contract(1.0,peps(row+1,col+1),shape(3,4),tmp7,shape(3,5),0.0,tmp8);

                  DArray<5> tmp5;
                  Contract(1.0,peps(row+1,col+1),shape(1,2,3,4),tmp8,shape(1,2,4,6),0.0,tmp5);

                  //add top
                  tmp8.clear();
                  Contract(1.0,peps(row+1,col),shape(4),tmp5,shape(0),0.0,tmp8);

                  tmp7.clear();
                  Contract(1.0,peps(row+1,col),shape(1,2,4),tmp8,shape(1,2,4),0.0,tmp7);

                  //add left side
                  tmp8.clear();
                  Contract(1.0,LI7,shape(0,1,6),tmp7,shape(2,0,6),0.0,tmp8);

                  N_eff.clear();
                  Permute(tmp8,shape(0,5,2,6,1,4,3,7),N_eff);

               }
               else{

                  //fill up the left side
                  DArray<8> tmp8;
                  Contract(1.0,LI7,shape(3,5),peps(row,col),shape(0,3),0.0,tmp8);

                  DArray<7> tmp7;
                  Contract(1.0,tmp8,shape(2,6,3),peps(row,col),shape(0,2,3),0.0,tmp7);

                  tmp8.clear();
                  Contract(1.0,tmp7,shape(1,3),peps(row+1,col),shape(0,3),0.0,tmp8);

                  DArray<5> tmp5;
                  Contract(1.0,tmp8,shape(0,5,6,3),peps(row+1,col),shape(0,1,2,3),0.0,tmp5);

                  //add top of right side
                  tmp8.clear();
                  Contract(1.0,tmp5,shape(4),peps(row+1,col+1),shape(0),0.0,tmp8);

                  tmp7.clear();
                  Contract(1.0,tmp8,shape(3,4,5),peps(row+1,col+1),shape(0,1,2),0.0,tmp7);

                  //add right side
                  tmp8.clear();
                  Contract(1.0,tmp7,shape(4,6,0),RI7,shape(0,1,6),0.0,tmp8);

                  N_eff.clear();
                  Permute(tmp8,shape(1,2,6,4,0,3,7,5),N_eff);

               }

            }
            else{//row = Ly - 1

               if(left){//left site of horizontal gate, so site (row,col) environment

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

               }
               else{//right site of horizontal gate, so site (row+1,col) environment

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

               }

            }

         }
         else if(dir == DIAGONAL_LURD){

            if(row == 0){

               if(left){//left site of diagonal gate, so site (row+1,col) environment

                  //add right peps to intermediate
                  DArray<8> tmp8;
                  Contract(1.0,peps(row,col+1),shape(1,4),RI7,shape(3,5),0.0,tmp8);

                  //add second peps
                  DArray<5> tmp5;
                  Contract(1.0,peps(row,col+1),shape(1,2,3,4),tmp8,shape(6,1,2,7),0.0,tmp5);

                  //add top environment
                  DArray<7> tmp7;
                  Contract(1.0,env.gt(0)[col],shape(3),tmp5,shape(2),0.0,tmp7);

                  //contract with left hand side
                  tmp8.clear();
                  Contract(1.0,LI7,shape(0,5,6),tmp7,shape(0,4,3),0.0,tmp8);

                  N_eff.clear();
                  Permute(tmp8,shape(0,4,2,6,1,5,3,7),N_eff);

               }
               else{//right site of horizontal gate, so site (row+1,col) environment

                  //add left peps to left intermediate
                  DArray<8> tmp8;
                  Contract(1.0,LI7,shape(2,4),peps(row+1,col),shape(0,3),0.0,tmp8);

                  //and again
                  DArray<7> tmp7;
                  Contract(1.0,tmp8,shape(1,6,2),peps(row+1,col),shape(0,2,3),0.0,tmp7);

                  //add top environment
                  DArray<5> tmp5;
                  Contract(1.0,tmp7,shape(0,5,3),env.gt(0)[col],shape(0,1,2),0.0,tmp5);

                  //now add right
                  DArray<6> tmp6;
                  Contract(1.0,tmp5,shape(4,3,2),RI7,shape(0,1,2),0.0,tmp6);

                  DArray<6> tmp6bis;
                  Permute(tmp6,shape(0,2,4,1,3,5),tmp6bis);

                  int DL = peps(row,col+1).shape(0);
                  int DR = peps(row,col+1).shape(4);

                  N_eff = tmp6bis.reshape_clear( shape(DL,D,1,DR,DL,D,1,DR) );

               }

            }
            else{//row == Ly - 2


            }

         }
         else{//last but not least! ---> diagonal LDRU

            if(row == 0){

               if(left){//left site of diagonal gate, so site (row,col) environment

                  //add right peps to intermediate
                  DArray<8> tmp8;
                  Contract(1.0,peps(row+1,col+1),shape(3,4),RI7,shape(4,2),0.0,tmp8);

                  //add second peps
                  DArray<7> tmp7;
                  Contract(1.0,peps(row+1,col+1),shape(2,3,4),tmp8,shape(2,5,4),0.0,tmp7);

                  //add top environment
                  DArray<5> tmp5;
                  Contract(1.0,env.gt(0)[col+1],shape(1,2,3),tmp7,shape(1,3,4),0.0,tmp5);

                  //contract with left hand side
                  DArray<6> tmp6;
                  Contract(1.0,LI7,shape(0,1,2),tmp5,shape(0,1,2),0.0,tmp6);

                  DArray<6> tmp6bis;
                  Permute(tmp6,shape(2,0,4,3,1,5),tmp6bis);

                  int DL = peps(row,col).shape(0);
                  int DR = peps(row,col).shape(4);

                  N_eff = tmp6bis.reshape_clear( shape(DL,D,1,DR,DL,D,1,DR) );

               }
               else{//right site of horizontal gate, so site (row+1,col+1) environment

                  //add left peps to left intermediate
                  DArray<8> tmp8;
                  Contract(1.0,LI7,shape(5,3),peps(row,col),shape(0,1),0.0,tmp8);

                  //and again
                  DArray<5> tmp5;
                  Contract(1.0,tmp8,shape(4,3,5,6),peps(row,col),shape(0,1,2,3),0.0,tmp5);

                  //add top environment
                  DArray<7> tmp7;
                  Contract(1.0,tmp5,shape(0),env.gt(0)[col+1],shape(0),0.0,tmp7);

                  //now add right
                  tmp8.clear();
                  Contract(1.0,tmp7,shape(6,2,3),RI7,shape(0,5,6),0.0,tmp8);

                  N_eff.clear();
                  Permute(tmp8,shape(0,2,6,4,1,3,7,5),N_eff);

               }

            }
            else{//row == Ly-2

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

         if(dir == VERTICAL){

            if(left){//bottom site

               //add upper peps to LI8
               DArray<9> tmp9;
               Contract(1.0,LI8,shape(4,2),peps(row+1,col),shape(0,1),0.0,tmp9);

               //and another
               DArray<8> tmp8;
               Contract(1.0,tmp9,shape(2,1,6),peps(row+1,col),shape(0,1,2),0.0,tmp8);

               //contract with right intermediate
               DArray<8> tmp8bis;
               Contract(1.0,tmp8,shape(0,7,5,3),RI8,shape(0,1,2,7),0.0,tmp8bis);

               N_eff.clear();
               Permute(tmp8bis,shape(0,3,6,4,1,2,7,5),N_eff);

            }
            else{//top site

               //add bottom peps  to intermediate right
               DArray<9> tmp9;
               Contract(1.0,peps(row,col),shape(3,4),RI8,shape(6,4),0.0,tmp9);

               //and another
               DArray<8> tmp8;
               Contract(1.0,peps(row,col),shape(2,3,4),tmp9,shape(2,7,6),0.0,tmp8);

               //now contract with left side
               DArray<8> tmp8bis;
               Contract(1.0,LI8,shape(0,5,6,7),tmp8,shape(4,0,2,7),0.0,tmp8bis);

               N_eff.clear();
               Permute(tmp8bis,shape(2,0,4,6,3,1,5,7),N_eff);

            }

         }
         else if(dir == HORIZONTAL){

            if(left){//left site of horizontal gate, so site (row,col) environment

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

            }
            else{//right site of horizontal gate, so site (row+1,col) environment

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

            }

         }
         else if(dir == DIAGONAL_LURD){

            if(left){//left site of gate, so site (row+1,col) environment

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

            }
            else{//right site of gate, so site (row,col+1) environment

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

            }

         }
         else{//dir == DIAGONAL_LDRU

            if(left){//left site of gate, so site (row,col) environment

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

            }
            else{//right site of gate, so site (row+1,col+1) environment

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

            }

         }

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

         if(dir == VERTICAL){

            if(row == 0){

               if(left){//bottom site

                  //paste top peps to left
                  DArray<8> tmp8;
                  Contract(1.0,LI7,shape(4,2),peps(row+1,col),shape(0,1),0.0,tmp8);

                  //and right operator
                  DArray<8> tmp8bis;
                  Contract(1.0,tmp8,shape(2,1,5),rop,shape(0,1,2),0.0,tmp8bis);

                  //then left operator
                  tmp8.clear();
                  Contract(1.0,tmp8bis,shape(1,6,5),lop,shape(0,1,3),0.0,tmp8);

                  //now add left and right together
                  DArray<5> tmp5;
                  Contract(1.0,tmp8,shape(0,4,3,7),R,shape(0,1,2,3),0.0,tmp5);

                  rhs.clear();
                  Permute(tmp5,shape(0,1,3,4,2),rhs);

               }
               else{//top site

                  DArray<8> tmp8;
                  Contract(1.0,peps(row,col),shape(4),R,shape(4),0.0,tmp8);

                  //add left operator
                  DArray<8> tmp8bis;
                  Contract(1.0,lop,shape(2,4,5),tmp8,shape(2,3,7),0.0,tmp8bis);

                  //and right operator
                  tmp8.clear();
                  Contract(1.0,rop,shape(3,4,5),tmp8bis,shape(2,1,6),0.0,tmp8);

                  //contract left and right
                  DArray<5> tmp5;
                  Contract(1.0,LI7,shape(0,1,3,5,6),tmp8,shape(6,1,0,3,4),0.0,tmp5);

                  rhs.clear();
                  Permute(tmp5,shape(1,0,3,4,2),rhs);

               }

            }
            else{ //row == Lx-2

               if(left){//bottom site

                  //paste right operator to right
                  DArray<9> tmp9;
                  Contract(1.0,rop,shape(5),R,shape(0),0.0,tmp9);

                  //and again
                  DArray<8> tmp8;
                  Contract(1.0,peps(row+1,col),shape(1,2,4),tmp9,shape(1,2,5),0.0,tmp8);

                  //attach left operator
                  DArray<8> tmp8bis;
                  Contract(1.0,lop,shape(1,3,5),tmp8,shape(4,3,5),0.0,tmp8bis);

                  //contract left and right
                  DArray<5> tmp5;
                  Contract(1.0,LI7,shape(0,1,2,4,6),tmp8bis,shape(5,3,0,2,7),0.0,tmp5);

                  rhs.clear();
                  Permute(tmp5,shape(0,3,1,4,2),rhs);

               }
               else{//top site

                  DArray<8> tmp8;
                  Contract(1.0,LI7,shape(3,5),peps(row,col),shape(0,3),0.0,tmp8);

                  //add left operator
                  DArray<8> tmp8bis;
                  Contract(1.0,tmp8,shape(2,6,3),lop,shape(0,2,4),0.0,tmp8bis);

                  //add right operator
                  tmp8.clear();
                  Contract(1.0,tmp8bis,shape(0,6,5),rop,shape(0,3,4),0.0,tmp8);

                  //now contract left and right
                  DArray<5> tmp5;
                  Contract(1.0,tmp8,shape(7,4,3,1),R,shape(0,2,3,4),0.0,tmp5);

                  rhs.clear();
                  Permute(tmp5,shape(0,2,1,4,3),rhs);

               }

            }

         }
         else if(dir == HORIZONTAL){

            if(row == 0){

               if(left){//left site of horizontal gate, so site (row,col) environment

                  //add right peps to intermediate
                  DArray<8> tmp8;
                  Contract(1.0,RI7,shape(4,6),peps(row,col+1),shape(1,4),0.0,tmp8);

                  //add right operator to tmp8
                  DArray<6> tmp6;
                  Contract(1.0,tmp8,shape(3,6,7,4),rop,shape(1,2,4,5),0.0,tmp6);

                  //attach LI7 to right side
                  DArray<7> tmp7;
                  Gemm(CblasTrans,CblasNoTrans,1.0,LI7,tmp6,0.0,tmp7);

                  //now paste left operator in
                  DArray<5> tmp5;
                  Contract(1.0,tmp7,shape(2,0,6,5),lop,shape(0,1,3,5),0.0,tmp5);

                  rhs.clear();
                  Permute(tmp5,shape(1,0,4,2,3),rhs);

               }
               else{//right site of horizontal gate, so site (row+1,col) environment

                  //add left peps to LI7
                  DArray<8> tmp8;
                  Contract(1.0,LI7,shape(6,4),peps(row,col),shape(0,1),0.0,tmp8);

                  //left operator another peps
                  DArray<6> tmp6;
                  Contract(1.0,tmp8,shape(4,3,5,6),lop,shape(0,1,2,4),0.0,tmp6);

                  //contract with RI7
                  DArray<7> tmp7;
                  Gemm(CblasTrans,CblasNoTrans,1.0,tmp6,RI7,0.0,tmp7);

                  //now paste right operator in
                  DArray<5> tmp5;
                  Contract(1.0,tmp7,shape(2,3,1,5),rop,shape(0,1,3,5),0.0,tmp5);

                  rhs.clear();
                  Permute(tmp5,shape(0,1,4,2,3),rhs);

               }

            }
            else if(row == Ly - 2){

               if(left){//left site of horizontal gate, so site (row,col) environment

                  //fill up the right side
                  DArray<8> tmp8;
                  Contract(1.0,peps(row,col+1),shape(3,4),RI7,shape(5,3),0.0,tmp8);

                  //right operator
                  DArray<8> tmp8bis;
                  Contract(1.0,rop,shape(2,4,5),tmp8,shape(2,6,5),0.0,tmp8bis);

                  DArray<9> tmp9;
                  Contract(1.0,peps(row+1,col+1),shape(3,4),tmp8bis,shape(4,6),0.0,tmp9);

                  DArray<6> tmp6;
                  Contract(1.0,peps(row+1,col+1),shape(1,2,3,4),tmp9,shape(1,2,4,7),0.0,tmp6);

                  //add top
                  tmp9.clear();
                  Contract(1.0,peps(row+1,col),shape(4),tmp6,shape(0),0.0,tmp9);

                  tmp8.clear();
                  Contract(1.0,peps(row+1,col),shape(1,2,4),tmp9,shape(1,2,4),0.0,tmp8);

                  //finally left operator
                  tmp8bis.clear();
                  Contract(1.0,lop,shape(1,3,5),tmp8,shape(3,5,4),0.0,tmp8bis);

                  //contract with LI7
                  DArray<5> tmp5;
                  Contract(1.0,LI7,shape(0,1,2,4,6),tmp8bis,shape(5,3,0,2,7),0.0,tmp5);

                  rhs.clear();
                  Permute(tmp5,shape(0,3,1,4,2),rhs);

               }
               else{//right site of horizontal gate, so site (row+1,col) environment

                  //fill up the left side
                  DArray<8> tmp8;
                  Contract(1.0,LI7,shape(3,5),peps(row,col),shape(0,3),0.0,tmp8);

                  //add left operator
                  DArray<8> tmp8bis;
                  Contract(1.0,tmp8,shape(2,6,3),lop,shape(0,2,4),0.0,tmp8bis);

                  DArray<9> tmp9;
                  Contract(1.0,tmp8bis,shape(1,3),peps(row+1,col),shape(0,3),0.0,tmp9);

                  DArray<6> tmp6;
                  Contract(1.0,tmp9,shape(0,6,7,3),peps(row+1,col),shape(0,1,2,3),0.0,tmp6);

                  //add top of right side
                  tmp9.clear();
                  Contract(1.0,tmp6,shape(5),peps(row+1,col+1),shape(0),0.0,tmp9);

                  tmp8.clear();
                  Contract(1.0,tmp9,shape(4,5,6),peps(row+1,col+1),shape(0,1,2),0.0,tmp8);

                  //add right operator
                  tmp8bis.clear();
                  Contract(1.0,tmp8,shape(3,4,2),rop,shape(0,1,3),0.0,tmp8bis);

                  //contract with RI7
                  DArray<5> tmp5;
                  Contract(1.0,tmp8bis,shape(2,4,7,6,0),RI7,shape(0,1,2,4,6),0.0,tmp5);

                  rhs.clear();
                  Permute(tmp5,shape(0,1,4,3,2),rhs);

               }

            }
            else{//row == Ly - 1

               if(left){//left site of horizontal gate, so site (row,col) environment

                  //add right peps to intermediate
                  DArray<8> tmp8;
                  Contract(1.0,peps(row,col+1),shape(3,4),RI7,shape(3,1),0.0,tmp8);

                  //add right operator to tmp8
                  DArray<6> tmp6;
                  Contract(1.0,rop,shape(1,2,4,5),tmp8,shape(1,2,4,3),0.0,tmp6);

                  //now paste left operator in
                  tmp8.clear();
                  Contract(1.0,lop,shape(5,3),tmp6,shape(0,1),0.0,tmp8);

                  //contract with left hand side
                  DArray<5> tmp5;
                  Contract(1.0,LI7,shape(0,2,4,5,6),tmp8,shape(0,3,5,6,7),0.0,tmp5);

                  rhs.clear();
                  Permute(tmp5,shape(0,2,1,4,3),rhs);

               }
               else{//right site of horizontal gate, so site (row+1,col) environment

                  //add left to intermediate
                  DArray<8> tmp8;
                  Contract(1.0,LI7,shape(1,3),peps(row,col),shape(0,3),0.0,tmp8);

                  //add left operator
                  DArray<6> tmp6;
                  Contract(1.0,tmp8,shape(0,5,6,1),lop,shape(0,1,2,4),0.0,tmp6);

                  //and right
                  tmp8.clear();
                  Contract(1.0,tmp6,shape(5,4),rop,shape(0,3),0.0,tmp8);

                  //contract with RI7 hand side
                  DArray<5> tmp5;
                  Contract(1.0,tmp8,shape(7,6,0,1,2),RI7,shape(0,2,4,5,6),0.0,tmp5);

                  rhs.clear();
                  Permute(tmp5,shape(0,1,4,3,2),rhs);

               }

            }

         }
         else if(dir == DIAGONAL_LURD){

            if(row == 0){

               if(left){

                  DArray<8> tmp8;
                  Contract(1.0,RI7,shape(4,6),peps(row,col+1),shape(1,4),0.0,tmp8);

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
                  DArray<5> tmp5;
                  Contract(1.0,tmp9,shape(0,5,8,7,3,4),tmp8,shape(0,1,3,7,6,5),0.0,tmp5);

                  rhs.clear();
                  Permute(tmp5,shape(0,3,1,4,2),rhs);

               }
               else{

                  //add upper-left peps to intermediate b_L
                  DArray<8> tmp8;
                  Contract(1.0,b_L,shape(2,4),peps(row+1,col),shape(0,3),0.0,tmp8);

                  //add left operator to tmp8
                  DArray<8> tmp8bis;
                  Contract(1.0,tmp8,shape(1,6,2),lop,shape(0,2,4),0.0,tmp8bis);

                  //add top environment to intermediate
                  DArray<6> tmp6;
                  Contract(1.0,tmp8bis,shape(0,5,3),env.gt(row)[col],shape(0,1,2),0.0,tmp6);

                  //add right operator to RI7
                  DArray<9> tmp9;
                  Contract(1.0,RI7,shape(3,5),rop,shape(1,5),0.0,tmp9);

                  //now contract left and right
                  DArray<5> tmp5;
                  Contract(1.0,tmp6,shape(5,4,3,2,0),tmp9,shape(0,1,7,2,5),0.0,tmp5);

                  rhs.clear();
                  Permute(tmp5,shape(0,1,4,2,3),rhs);

               }

            }

         }
         else{//diagonal ldru

            if(row == 0){

               if(left){

                  //add top right peps to b_R
                  DArray<8> tmp8;
                  Contract(1.0,peps(row+1,col+1),shape(3,4),b_R,shape(4,2),0.0,tmp8);

                  //add right operator to tmp8
                  DArray<8> tmp8bis;
                  Contract(1.0,rop,shape(2,4,5),tmp8,shape(2,5,4),0.0,tmp8bis);

                  //add top environment to intermediate
                  DArray<6> tmp6;
                  Contract(1.0,env.gt(row)[col+1],shape(1,2,3),tmp8bis,shape(1,4,5),0.0,tmp6);

                  //add left operator to LI7
                  DArray<9> tmp9;
                  Contract(1.0,LI7,shape(5,3),lop,shape(0,1),0.0,tmp9);

                  //contract both sides to form right hand side of equation
                  DArray<5> tmp5;
                  Contract(1.0,tmp9,shape(0,1,2,6,8),tmp6,shape(0,1,3,2,4),0.0,tmp5);

                  rhs.clear();
                  Permute(tmp5,shape(1,0,3,4,2),rhs);

               }
               else{

                  DArray<8> tmp8;
                  Contract(1.0,LI7,shape(6,4),peps(row,col),shape(0,1),0.0,tmp8);

                  //add left operator to tmp8
                  DArray<6> tmp6;
                  Contract(1.0,tmp8,shape(4,3,5,6),lop,shape(0,1,2,4),0.0,tmp6);

                  tmp8.clear();
                  Gemm(CblasTrans,CblasNoTrans,1.0,tmp6,env.gt(row)[col+1],0.0,tmp8);

                  //add right operator to b_R
                  DArray<9> tmp9;
                  Contract(1.0,rop,shape(4,5),b_R,shape(3,1),0.0,tmp9);

                  //contract both sides to form right hand side of equation
                  DArray<5> tmp5;
                  Contract(1.0,tmp8,shape(7,5,0,4,3,2),tmp9,shape(4,1,0,7,3,8),0.0,tmp5);

                  rhs.clear();
                  Permute(tmp5,shape(0,1,4,3,2),rhs);

               }

            }

         }

      }

   /**
    * construct the right hand side needed for the linear system of any gate direction specified by 'dir'
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

         if(dir ==  VERTICAL){

            if(left){

               //add upper peps to LI8
               DArray<9> tmp9;
               Contract(1.0,LI8,shape(4,2),peps(row+1,col),shape(0,1),0.0,tmp9);

               //add right operator to intermediate
               DArray<9> tmp9bis;
               Contract(1.0,tmp9,shape(2,1,6),rop,shape(0,1,2),0.0,tmp9bis);

               //next add left operator
               tmp9.clear();
               Contract(1.0,tmp9bis,shape(1,7,6),lop,shape(0,1,3),0.0,tmp9);

               DArray<5> tmp5;
               Contract(1.0,tmp9,shape(0,5,4,8,7,2),RI8,shape(0,1,2,3,5,7),0.0,tmp5);

               rhs.clear();
               Permute(tmp5,shape(0,1,4,3,2),rhs);

            }
            else{

               //add bottom peps  to intermediate right
               DArray<9> tmp9;
               Contract(1.0,peps(row,col),shape(3,4),RI8,shape(6,4),0.0,tmp9);

               //add left operator to intermediate right
               DArray<9> tmp9bis;
               Contract(1.0,lop,shape(2,4,5),tmp9,shape(2,7,6),0.0,tmp9bis);

               //and right operator
               tmp9.clear();
               Contract(1.0,rop,shape(3,4,5),tmp9bis,shape(2,1,6),0.0,tmp9);

               //contract with left hand side
               DArray<5> tmp5;
               Contract(1.0,LI8,shape(0,1,3,5,6,7),tmp9,shape(6,1,0,3,4,8),0.0,tmp5);

               rhs.clear();
               Permute(tmp5,shape(1,0,3,4,2),rhs);

            }

         }
         else if(dir == HORIZONTAL){

            if(left){//left site of horizontal gate, so site (row,col) environment

               //add right peps to intermediate
               DArray<9> tmp9;
               Contract(1.0,RI8,shape(4,6),peps(row,col+1),shape(1,4),0.0,tmp9);

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

               //add left peps to LI8
               DArray<9> tmp9;
               Contract(1.0,LI8,shape(6,4),peps(row,col),shape(0,1),0.0,tmp9);

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
    * @param R_l vector of DArray<2> objects containing the inverse R of the QR decompostion (left site)
    * @param R_r vector of DArray<2> objects containing the inverse R of the QR decompostion (right site)
    */
   template<>
      void canonicalize(const PROP_DIR &dir,int row,int col,PEPS<double> &peps,DArray<5> &L, DArray<5> &R, 

            DArray<7> &LI7,DArray<7> &RI7,std::vector< DArray<2> > &R_l,std::vector< DArray<2> > &R_r){

         // ----------------------------//
         // --- (A) ---- LEFT SITE ---- //
         // ----------------------------//

         //calculate the environment for the left site
         DArray<8> N_eff;
         calc_N_eff(dir,row,col,peps,N_eff,L,R,LI7,RI7,true);

         //eigenvalues
         DArray<1> eig;
         diagonalize(N_eff,eig);

#ifdef _DEBUG
         cout << endl;
         cout << "before canonicalization" << endl;
         cout << endl;
         cout << "(left)\t" << std::scientific << eig(eig.size() - 1) / eig(0) << endl;
#endif

         //are needed for the calculation of the positive approximant: physical dimension is last index (4)
         DArray<5> X;
         get_X(N_eff,eig,X);

         DArray<5> X_copy;

         //set left and right indices
         int lrow = row;
         int rrow = row;

         int lcol = col;
         int rcol = col;

         if(dir == VERTICAL)
            rrow++;
         else if(dir == HORIZONTAL)
            rcol++;
         else if(dir == DIAGONAL_LURD){

            lrow++;
            rcol++;

         }
         else{//DIAGONAL_LDRU

            rrow++;
            rcol++;

         }

         // ----------------------------------------- //
         // --- canoncalize left-site environment --- //
         // ----------------------------------------- //

         if(N_eff.shape(0) > 1){//left

            X_copy = X;

            //QR: watch out, LQ decomposition
            Gelqf(R_l[0],X_copy);

            //add to left side of tensor
            DArray<5> tmp5;
            Contract(1.0,R_l[0],shape(0),peps(lrow,lcol),shape(0),0.0,tmp5);

            peps(lrow,lcol) = std::move(tmp5);

            //add  inverse to environment, for upper and lower layer
            invert(R_l[0]);

            //for horizontal and vertical: change LI7
            if(dir == VERTICAL || dir == HORIZONTAL){

               if(row == 0){

                  DArray<7> tmp7;
                  Contract(1.0,LI7,shape(5),R_l[0],shape(1),0.0,tmp7);

                  //and again
                  Contract(1.0,tmp7,shape(5),R_l[0],shape(1),0.0,LI7);

               }
               else if(row == Ly - 2){

                  DArray<7> tmp7;
                  Contract(1.0,LI7,shape(2),R_l[0],shape(1),0.0,tmp7);

                  //and again
                  LI7.clear();
                  Contract(1.0,tmp7,shape(2),R_l[0],shape(1),0.0,LI7);

                  tmp7.clear();
                  Permute(LI7,shape(0,1,5,6,2,3,4),tmp7);

                  LI7 = std::move(tmp7);

               }
               else{//row == Ly - 1

                  DArray<7> tmp7;
                  Contract(1.0,LI7,shape(0),R_l[0],shape(1),0.0,tmp7);

                  //and again
                  LI7.clear();
                  Contract(1.0,tmp7,shape(0),R_l[0],shape(1),0.0,LI7);

                  tmp7.clear();
                  Permute(LI7,shape(5,6,0,1,2,3,4),tmp7);

                  LI7 = std::move(tmp7);

               }

            }
            else{//diagonal gates have to recalculate LI7 and RI7,

               //add to right side of tensor
               DArray<5> tmp5;
               Contract(1.0,peps(lrow,lcol-1),shape(4),R_l[0],shape(1),0.0,tmp5);

               peps(lrow,lcol-1) = std::move(tmp5);

               //for consistency with right environment construction, multiply with LI7 
               if(dir == DIAGONAL_LURD){

                  //add inverse to L
                  tmp5.clear();
                  Contract(1.0,L,shape(1),R_l[0],shape(1),0.0,tmp5);

                  //and again
                  L.clear();
                  Contract(1.0,tmp5,shape(1),R_l[0],shape(1),0.0,L);

                  Permute(L,shape(0,3,4,1,2),tmp5);

                  L = std::move(tmp5);

                  //and add inverse to LI7
                  DArray<7> tmp7;
                  Contract(1.0,LI7,shape(1),R_l[0],shape(1),0.0,tmp7);

                  //and again
                  Contract(1.0,tmp7,shape(1),R_l[0],shape(1),0.0,LI7);

                  Permute(LI7,shape(0,5,6,1,2,3,4),tmp7);

                  LI7 = std::move(tmp7);

               }
               else{//DIAGONAL_LDRU

                  //add inverse to L
                  tmp5.clear();
                  Contract(1.0,L,shape(3),R_l[0],shape(1),0.0,tmp5);

                  //and again
                  L.clear();
                  Contract(1.0,tmp5,shape(3),R_l[0],shape(1),0.0,L);

                  //and add inverse to LI7
                  DArray<7> tmp7;
                  Contract(1.0,LI7,shape(5),R_l[0],shape(1),0.0,tmp7);

                  //and again
                  Contract(1.0,tmp7,shape(5),R_l[0],shape(1),0.0,LI7);

               }

            }

         }

         if(N_eff.shape(1) > 1 && dir != VERTICAL){//up

            Permute(X,shape(0,2,3,4,1),X_copy);

            //QR
            Geqrf(X_copy,R_l[1]);

            //add to up side of tensor
            DArray<5> tmp5;
            Contract(1.0,peps(lrow,lcol),shape(1),R_l[1],shape(1),0.0,tmp5);

            Permute(tmp5,shape(0,4,1,2,3),peps(lrow,lcol));

            //add  inverse to environment, for upper and lower layer
            invert(R_l[1]);

            if(dir == HORIZONTAL){

               if(row == 0){

                  DArray<7> tmp7;
                  Contract(1.0,LI7,shape(3),R_l[1],shape(0),0.0,tmp7);

                  //and again
                  LI7.clear();
                  Contract(1.0,tmp7,shape(3),R_l[1],shape(0),0.0,LI7);

                  tmp7.clear();
                  Permute(LI7,shape(0,1,2,5,6,3,4),tmp7);

                  LI7 = std::move(tmp7);

               }
               else if(row == Ly - 2){

                  DArray<5> tmp5;
                  Contract(1.0,peps(row+1,col),shape(3),R_l[1],shape(0),0.0,tmp5);

                  Permute(tmp5,shape(0,1,2,4,3),peps(row+1,col));

               }

            }
            else if(dir == DIAGONAL_LURD){//diagonal gates recalculate LI7 and RI7

               //add inverse to top environment
               DArray<4> tmp4;
               Contract(1.0,env.gt(0)[lcol],shape(1),R_l[1],shape(0),0.0,tmp4);

               //and again
               env.gt(0)[lcol].clear();
               Contract(1.0,tmp4,shape(1),R_l[1],shape(0),0.0,env.gt(0)[lcol]);

               tmp4.clear();
               Permute(env.gt(0)[lcol],shape(0,2,3,1),tmp4);

               env.gt(0)[lcol] = std::move(tmp4);

            }
            else{//diagonal gates recalculate LI7 and RI7, multiply the upper tensor for LDRU

               //add inverse to down side of upper tensor
               DArray<5> tmp5;
               Contract(1.0,peps(lrow+1,lcol),shape(3),R_l[1],shape(0),0.0,tmp5);

               Permute(tmp5,shape(0,1,2,4,3),peps(lrow+1,lcol));

               //add to LI7 for consistency
               DArray<7> tmp7;
               Contract(1.0,LI7,shape(3),R_l[1],shape(0),0.0,tmp7);

               //and again
               LI7.clear();
               Contract(1.0,tmp7,shape(3),R_l[1],shape(0),0.0,LI7);

               tmp7.clear();
               Permute(LI7,shape(0,1,2,5,6,3,4),tmp7);

               LI7 = std::move(tmp7);

            }

         }

         if(N_eff.shape(2) > 1){//down

            Permute(X,shape(0,1,3,4,2),X_copy);

            //QR
            Geqrf(X_copy,R_l[2]);

            //add to down side of tensor
            DArray<5> tmp5;
            Contract(1.0,peps(lrow,lcol),shape(3),R_l[2],shape(1),0.0,tmp5);

            Permute(tmp5,shape(0,1,2,4,3),peps(lrow,lcol));

            //add  inverse to environment, for upper and lower layer
            invert(R_l[2]);

            if(dir == VERTICAL){

               DArray<7> tmp7;
               Contract(1.0,LI7,shape(4),R_l[2],shape(0),0.0,tmp7);

               //and again
               LI7.clear();
               Contract(1.0,tmp7,shape(4),R_l[2],shape(0),0.0,LI7);

               tmp7.clear();
               Permute(LI7,shape(0,1,2,3,5,6,4),tmp7);

               LI7 = std::move(tmp7);

            }
            else if(dir == HORIZONTAL){

               if(row == Ly - 2){

                  DArray<7> tmp7;
                  Contract(1.0,LI7,shape(4),R_l[2],shape(0),0.0,tmp7);

                  //and again
                  LI7.clear();
                  Contract(1.0,tmp7,shape(4),R_l[2],shape(0),0.0,LI7);

                  tmp7.clear();
                  Permute(LI7,shape(0,1,2,3,5,6,4),tmp7);

                  LI7 = std::move(tmp7);

               }
               else{//row == Ly - 1

                  DArray<7> tmp7;
                  Contract(1.0,LI7,shape(2),R_l[2],shape(0),0.0,tmp7);

                  //and again
                  LI7.clear();
                  Contract(1.0,tmp7,shape(2),R_l[2],shape(0),0.0,LI7);

                  tmp7.clear();
                  Permute(LI7,shape(0,1,5,6,2,3,4),tmp7);

                  LI7 = std::move(tmp7);

               }

            }
            else if(dir == DIAGONAL_LURD){//diagonal lurd, change the actual tensors

               //add to upper side of lower tensor
               DArray<5> tmp5;
               Contract(1.0,peps(lrow-1,lcol),shape(1),R_l[2],shape(0),0.0,tmp5);

               Permute(tmp5,shape(0,4,1,2,3),peps(lrow-1,lcol));

               DArray<7> tmp7;
               Contract(1.0,LI7,shape(3),R_l[2],shape(0),0.0,tmp7);

               //and again
               LI7.clear();
               Contract(1.0,tmp7,shape(3),R_l[2],shape(0),0.0,LI7);

               tmp7.clear();
               Permute(LI7,shape(0,1,2,5,6,3,4),tmp7);

               LI7 = std::move(tmp7);

            }

         }

         if(N_eff.shape(3) > 1 && dir != HORIZONTAL){//right

            Permute(X,shape(0,1,2,4,3),X_copy);

            //QR
            Geqrf(X_copy,R_l[3]);

            //add to right side of tensor
            DArray<5> tmp5;
            Contract(1.0,peps(lrow,lcol),shape(4),R_l[3],shape(1),0.0,tmp5);

            peps(lrow,lcol) = std::move(tmp5);

            //add  inverse to environment, for upper and lower layer
            invert(R_l[3]);

            if(dir == VERTICAL){

               if(row == 0){

                  tmp5.clear();
                  Contract(1.0,R,shape(3),R_l[3],shape(0),0.0,tmp5);

                  //and again
                  Contract(1.0,tmp5,shape(3),R_l[3],shape(0),0.0,R);

               }
               else{

                  tmp5.clear();
                  Contract(1.0,R,shape(2),R_l[3],shape(0),0.0,tmp5);

                  //and again
                  R.clear();
                  Contract(1.0,tmp5,shape(2),R_l[3],shape(0),0.0,R);

                  tmp5.clear();
                  Permute(R,shape(0,1,3,4,2),tmp5);

                  R = std::move(tmp5);

               }

            }
            else{//diagonal, add inverse to right tensor

               //add to right side of tensor
               DArray<5> tmp5;
               Contract(1.0,R_l[3],shape(0),peps(lrow,lcol+1),shape(0),0.0,tmp5);

               peps(lrow,lcol+1) = std::move(tmp5);

               //change RI7 as well for consistency when calculating right environment
               if(dir == DIAGONAL_LURD){

                  DArray<7> tmp7;
                  Contract(1.0,RI7,shape(1),R_l[3],shape(0),0.0,tmp7);

                  //and again
                  RI7.clear();
                  Contract(1.0,tmp7,shape(1),R_l[3],shape(0),0.0,RI7);

                  tmp7.clear();
                  Permute(RI7,shape(0,5,6,1,2,3,4),tmp7);

                  RI7 = std::move(tmp7);

               }
               else{//diagonal ldru

                  DArray<7> tmp7;
                  Contract(1.0,RI7,shape(5),R_l[3],shape(0),0.0,tmp7);

                  //and again
                  RI7.clear();
                  Contract(1.0,tmp7,shape(5),R_l[3],shape(0),0.0,RI7);

               }

            }

         }

         // -----------------------------//
         // --- (B) ---- RIGHT SITE ---- //
         // -----------------------------//

         //calculate the environment for the right site
         N_eff.clear();
         calc_N_eff(dir,row,col,peps,N_eff,L,R,LI7,RI7,false);

         //eigenvalues
         eig.clear();
         diagonalize(N_eff,eig);

#ifdef _DEBUG
         cout << "(right)\t" << std::scientific << eig(eig.size() - 1) / eig(0) << endl;
#endif

         //are needed for the calculation of the positive approximant: physical dimension is last index (4)
         X.clear();
         get_X(N_eff,eig,X);

         if(N_eff.shape(0) > 1 && dir != HORIZONTAL){//left

            X_copy = X;

            //QR: watch out, LQ decomposition
            Gelqf(R_r[0],X_copy);

            //add to left side of tensor
            DArray<5> tmp5;
            Contract(1.0,R_r[0],shape(0),peps(rrow,rcol),shape(0),0.0,tmp5);

            peps(rrow,rcol) = std::move(tmp5);

            //add  inverse to environment, for upper and lower layer
            invert(R_r[0]);

            if(dir == VERTICAL){

               if(row == 0){

                  DArray<7> tmp7;
                  Contract(1.0,LI7,shape(3),R_r[0],shape(1),0.0,tmp7);

                  //and again
                  Contract(1.0,tmp7,shape(3),R_r[0],shape(1),0.0,LI7);

                  tmp7.clear();
                  Permute(LI7,shape(0,1,2,5,6,3,4),tmp7);

                  LI7 = std::move(tmp7);

               }
               else{

                  DArray<7> tmp7;
                  Contract(1.0,LI7,shape(0),R_r[0],shape(1),0.0,tmp7);

                  //and again
                  LI7.clear();
                  Contract(1.0,tmp7,shape(0),R_r[0],shape(1),0.0,LI7);

                  tmp7.clear();
                  Permute(LI7,shape(5,6,0,1,2,3,4),tmp7);

                  LI7 = std::move(tmp7);

               }

            }
            else{//diagonal

               //add to right side of left tensor
               DArray<5> tmp5;
               Contract(1.0,peps(rrow,rcol-1),shape(4),R_r[0],shape(1),0.0,tmp5);

               peps(rrow,rcol-1) = std::move(tmp5);

            }

         }

         if(N_eff.shape(1) > 1){//up

            Permute(X,shape(0,2,3,4,1),X_copy);

            //QR
            Geqrf(X_copy,R_r[1]);

            //add to up side of tensor
            DArray<5> tmp5;
            Contract(1.0,peps(rrow,rcol),shape(1),R_r[1],shape(1),0.0,tmp5);

            Permute(tmp5,shape(0,4,1,2,3),peps(rrow,rcol));

            //add  inverse to environment, for upper and lower layer
            invert(R_r[1]);

            if(dir == VERTICAL){//only row = 0

               DArray<7> tmp7;
               Contract(1.0,LI7,shape(1),R_r[1],shape(0),0.0,tmp7);

               //and again
               LI7.clear();
               Contract(1.0,tmp7,shape(1),R_r[1],shape(0),0.0,LI7);

               tmp7.clear();
               Permute(LI7,shape(0,5,6,1,2,3,4),tmp7);

               LI7 = std::move(tmp7);

            }
            else if(dir == HORIZONTAL){

               if(row == 0){

                  DArray<7> tmp7;
                  Contract(1.0,RI7,shape(3),R_r[1],shape(0),0.0,tmp7);

                  //and again
                  RI7.clear();
                  Contract(1.0,tmp7,shape(3),R_r[1],shape(0),0.0,RI7);

                  tmp7.clear();
                  Permute(RI7,shape(0,1,2,5,6,3,4),tmp7);

                  RI7 = std::move(tmp7);

               }
               else if(row == Ly - 2){

                  DArray<5> tmp5;
                  Contract(1.0,peps(row+1,col+1),shape(3),R_r[1],shape(0),0.0,tmp5);

                  Permute(tmp5,shape(0,1,2,4,3),peps(row+1,col+1));

               }

            }
            else if(dir == DIAGONAL_LURD){//diagonal, multiply inverse with down side of upper tensor

               //add to up side of tensor
               DArray<5> tmp5;
               Contract(1.0,peps(rrow+1,rcol),shape(3),R_r[1],shape(0),0.0,tmp5);

               Permute(tmp5,shape(0,1,2,4,3),peps(rrow+1,rcol));

            }
            else{//diagonal LDRU, multiply inverse with top environment

               //add inverse to top environment
               DArray<4> tmp4;
               Contract(1.0,env.gt(0)[rcol],shape(1),R_r[1],shape(0),0.0,tmp4);

               //and again
               env.gt(0)[rcol].clear();
               Contract(1.0,tmp4,shape(1),R_r[1],shape(0),0.0,env.gt(0)[rcol]);

               tmp4.clear();
               Permute(env.gt(0)[rcol],shape(0,2,3,1),tmp4);

               env.gt(0)[rcol] = std::move(tmp4);

            }

         }

         if(N_eff.shape(2) > 1 & dir != VERTICAL){//down

            Permute(X,shape(0,1,3,4,2),X_copy);

            //QR
            Geqrf(X_copy,R_r[2]);

            //add to down side of tensor
            DArray<5> tmp5;
            Contract(1.0,peps(rrow,rcol),shape(3),R_r[2],shape(1),0.0,tmp5);

            Permute(tmp5,shape(0,1,2,4,3),peps(rrow,rcol));

            //add  inverse to environment, for upper and lower layer
            invert(R_r[2]);

            if(dir == HORIZONTAL){

               if(row == Ly - 2){

                  DArray<7> tmp7;
                  Contract(1.0,RI7,shape(4),R_r[2],shape(0),0.0,tmp7);

                  //and again
                  RI7.clear();
                  Contract(1.0,tmp7,shape(4),R_r[2],shape(0),0.0,RI7);

                  tmp7.clear();
                  Permute(RI7,shape(0,1,2,3,5,6,4),tmp7);

                  RI7 = std::move(tmp7);

               }
               else{//row == Ly - 1

                  DArray<7> tmp7;
                  Contract(1.0,RI7,shape(2),R_r[2],shape(0),0.0,tmp7);

                  //and again
                  RI7.clear();
                  Contract(1.0,tmp7,shape(2),R_r[2],shape(0),0.0,RI7);

                  tmp7.clear();
                  Permute(RI7,shape(0,1,5,6,2,3,4),tmp7);

                  RI7 = std::move(tmp7);

               }

            }
            else if(dir == DIAGONAL_LURD){//diagonal

               //add inverse to up side of lower tensor
               DArray<5> tmp5;
               Contract(1.0,peps(rrow-1,rcol),shape(1),R_r[2],shape(0),0.0,tmp5);

               Permute(tmp5,shape(0,4,1,2,3),peps(rrow-1,rcol));

            }

         }

         if(N_eff.shape(3) > 1){//right

            Permute(X,shape(0,1,2,4,3),X_copy);

            //QR
            Geqrf(X_copy,R_r[3]);

            //add to right side of tensor
            DArray<5> tmp5;
            Contract(1.0,peps(rrow,rcol),shape(4),R_r[3],shape(1),0.0,tmp5);

            peps(rrow,rcol) = std::move(tmp5);

            //add  inverse to environment, for upper and lower layer
            invert(R_r[3]);

            if(dir == VERTICAL){

               if(row == 0){

                  tmp5.clear();
                  Contract(1.0,R,shape(1),R_r[3],shape(0),0.0,tmp5);

                  //and again
                  R.clear();
                  Contract(1.0,tmp5,shape(1),R_r[3],shape(0),0.0,R);

                  tmp5.clear();
                  Permute(R,shape(0,3,4,1,2),tmp5);

                  R = std::move(tmp5);

               }
               else{

                  tmp5.clear();
                  Contract(1.0,R,shape(0),R_r[3],shape(0),0.0,tmp5);

                  //and again
                  R.clear();
                  Contract(1.0,tmp5,shape(0),R_r[3],shape(0),0.0,R);

                  tmp5.clear();
                  Permute(R,shape(3,4,0,1,2),tmp5);

                  R = std::move(tmp5);

               }

            }
            else if(dir == HORIZONTAL){

               if(row == 0){

                  DArray<7> tmp7;
                  Contract(1.0,RI7,shape(5),R_r[3],shape(0),0.0,tmp7);

                  //and again
                  RI7.clear();
                  Contract(1.0,tmp7,shape(5),R_r[3],shape(0),0.0,RI7);

               }
               else if(row == Ly - 2){

                  DArray<7> tmp7;
                  Contract(1.0,RI7,shape(2),R_r[3],shape(0),0.0,tmp7);

                  //and again
                  RI7.clear();
                  Contract(1.0,tmp7,shape(2),R_r[3],shape(0),0.0,RI7);

                  tmp7.clear();
                  Permute(RI7,shape(0,1,5,6,2,3,4),tmp7);

                  RI7 = std::move(tmp7);

               }
               else{//row == Ly - 1

                  DArray<7> tmp7;
                  Contract(1.0,RI7,shape(0),R_r[3],shape(0),0.0,tmp7);

                  //and again
                  RI7.clear();
                  Contract(1.0,tmp7,shape(0),R_r[3],shape(0),0.0,RI7);

                  tmp7.clear();
                  Permute(RI7,shape(5,6,0,1,2,3,4),tmp7);

                  RI7 = std::move(tmp7);

               }

            }
            else if(dir == DIAGONAL_LURD){//change the R tensor

               DArray<5> tmp5;
               Contract(1.0,R,shape(3),R_r[3],shape(0),0.0,tmp5);

               //and again
               R.clear();
               Contract(1.0,tmp5,shape(3),R_r[3],shape(0),0.0,R);

            }
            else{//diagonal ldru

               DArray<5> tmp5;
               Contract(1.0,R,shape(1),R_r[3],shape(0),0.0,tmp5);

               //and again
               R.clear();
               Contract(1.0,tmp5,shape(1),R_r[3],shape(0),0.0,R);

               tmp5.clear();
               Permute(R,shape(0,3,4,1,2),tmp5);

               R = std::move(tmp5);

            }

         }

      }

   /**
    * quasi canonicalize the environment of the sites to be updated, for stability reasons in the program.
    * @param dir vertical, horizontal,diagonal lurd or diagonal ldru update
    * @param row the row index of the bottom site
    * @param col column index of the vertical column
    * @param peps, full PEPS object before update
    * @param L Left environment contraction
    * @param R Right environment contraction
    * @param LI8 left intermediate object
    * @param RI8 right intermediate object
    * @param R_l vector of DArray<2> objects containing the inverse R of the QR decompostion (left site)
    * @param R_r vector of DArray<2> objects containing the inverse R of the QR decompostion (right site)
    */
   template<>
      void canonicalize(const PROP_DIR &dir,int row,int col,PEPS<double> &peps, DArray<6> &LO, DArray<6> &RO,

            DArray<8> &LI8,DArray<8> &RI8,std::vector< DArray<2> > &R_l,std::vector< DArray<2> > &R_r){

         // ----------------------------//
         // --- (A) ---- LEFT SITE ---- //
         // ----------------------------//

         //calculate the environment for the left site
         DArray<8> N_eff;
         calc_N_eff(dir,row,col,peps,N_eff,LO,RO,LI8,RI8,true);

         //eigenvalues
         DArray<1> eig;
         diagonalize(N_eff,eig);

#ifdef _DEBUG
         cout << endl;
         cout << "before canonicalization" << endl;
         cout << endl;

         cout << "(left)\t" << std::scientific << eig(eig.size() - 1) / eig(0) << endl;
#endif

         //are needed for the calculation of the positive approximant: physical dimension is last index (4)
         DArray<5> X;
         get_X(N_eff,eig,X);

         DArray<5> X_copy;

         //set left and right indices
         int lrow = row;
         int rrow = row;

         int lcol = col;
         int rcol = col;

         if(dir == VERTICAL)
            rrow++;
         else if(dir == HORIZONTAL)
            rcol++;
         else if(dir == DIAGONAL_LURD){

            lrow++;
            rcol++;

         }
         else{//DIAGONAL_LDRU

            rrow++;
            rcol++;

         }

         // ----------------------------------------- //
         // --- canoncalize left-site environment --- //
         // ----------------------------------------- //

         if(N_eff.shape(0) > 1){//left

            X_copy = X;

            //QR: watch out, LQ decomposition
            Gelqf(R_l[0],X_copy);

            //add to tensor
            DArray<5> tmp5;
            Contract(1.0,R_l[0],shape(0),peps(lrow,lcol),shape(0),0.0,tmp5);

            peps(lrow,lcol) = std::move(tmp5);

            //add  inverse to environment, for upper and lower layer
            invert(R_l[0]);

            if(dir == HORIZONTAL || dir == VERTICAL){

               DArray<8> tmp8;
               Contract(1.0,LI8,shape(5),R_l[0],shape(1),0.0,tmp8);

               //and again
               LI8.clear();
               Contract(1.0,tmp8,shape(5),R_l[0],shape(1),0.0,LI8);

               Permute(LI8,shape(0,1,2,3,4,6,7,5),tmp8);

               LI8 = std::move(tmp8);

            }
            else if(dir == DIAGONAL_LURD){//change LO and LI8!

               //LO
               DArray<6> tmp6;
               Contract(1.0,LO,shape(1),R_l[0],shape(1),0.0,tmp6);

               //and again
               LO.clear();
               Contract(1.0,tmp6,shape(1),R_l[0],shape(1),0.0,LO);

               Permute(LO,shape(0,4,5,1,2,3),tmp6);

               LO = std::move(tmp6);

               //LI8
               DArray<8> tmp8;
               Contract(1.0,LI8,shape(1),R_l[0],shape(1),0.0,tmp8);

               //and again
               LI8.clear();
               Contract(1.0,tmp8,shape(1),R_l[0],shape(1),0.0,LI8);

               Permute(LI8,shape(0,6,7,1,2,3,4,5),tmp8);

               LI8 = std::move(tmp8);

            }
            else{//diagonal LDRU

               //LI8
               DArray<8> tmp8;
               Contract(1.0,LI8,shape(5),R_l[0],shape(1),0.0,tmp8);

               //and again
               LI8.clear();
               Contract(1.0,tmp8,shape(5),R_l[0],shape(1),0.0,LI8);

               Permute(LI8,shape(0,1,2,3,4,6,7,5),tmp8);

               LI8 = std::move(tmp8);

            }

         }

         if(N_eff.shape(1) > 1 && dir != VERTICAL){//up

            Permute(X,shape(0,2,3,4,1),X_copy);

            //QR
            Geqrf(X_copy,R_l[1]);

            //add to up side of tensor
            DArray<5> tmp5;
            Contract(1.0,peps(lrow,lcol),shape(1),R_l[1],shape(1),0.0,tmp5);

            Permute(tmp5,shape(0,4,1,2,3),peps(lrow,lcol));

            //add  inverse to environment, for upper and lower layer
            invert(R_l[1]);

            if(dir == HORIZONTAL){

               DArray<8> tmp8;
               Contract(1.0,LI8,shape(3),R_l[1],shape(0),0.0,tmp8);

               //and again
               LI8.clear();
               Contract(1.0,tmp8,shape(3),R_l[1],shape(0),0.0,LI8);

               tmp8.clear();
               Permute(LI8,shape(0,1,2,6,7,3,4,5),tmp8);

               LI8 = std::move(tmp8);

            }
            else if(dir == DIAGONAL_LURD){//change top environment!

               DArray<4> tmp4;
               Contract(1.0,env.gt(row)[col],shape(1),R_l[1],shape(0),0.0,tmp4);

               //and again
               env.gt(row)[col].clear();
               Contract(1.0,tmp4,shape(1),R_l[1],shape(0),0.0,env.gt(row)[col]);

               tmp4.clear();
               Permute(env.gt(row)[col],shape(0,2,3,1),tmp4);

               env.gt(row)[col] = std::move(tmp4);

            }
            else{//diagonal ldru: change LI8 and upper peps

               //LI8
               DArray<8> tmp8;
               Contract(1.0,LI8,shape(3),R_l[1],shape(0),0.0,tmp8);

               //and again
               LI8.clear();
               Contract(1.0,tmp8,shape(3),R_l[1],shape(0),0.0,LI8);

               tmp8.clear();
               Permute(LI8,shape(0,1,2,6,7,3,4,5),tmp8);

               LI8 = std::move(tmp8);

               //upper peps

               //add to up side of tensor
               DArray<5> tmp5;
               Contract(1.0,peps(lrow+1,lcol),shape(3),R_l[1],shape(1),0.0,tmp5);

               Permute(tmp5,shape(0,1,2,4,3),peps(lrow,lcol));

            }

         }

         if(N_eff.shape(2) > 1){//down

            Permute(X,shape(0,1,3,4,2),X_copy);

            //QR
            Geqrf(X_copy,R_l[2]);

            //add to down side of tensor
            DArray<5> tmp5;
            Contract(1.0,peps(lrow,lcol),shape(3),R_l[2],shape(1),0.0,tmp5);

            Permute(tmp5,shape(0,1,2,4,3),peps(lrow,lcol));

            //add  inverse to environment, for upper and lower layer
            invert(R_l[2]);

            if(dir == HORIZONTAL || dir == DIAGONAL_LDRU){

               DArray<4> tmp4;
               Contract(1.0,env.gb(row-1)[col],shape(1),R_l[2],shape(0),0.0,tmp4);

               //and again
               env.gb(row-1)[col].clear();
               Contract(1.0,tmp4,shape(1),R_l[2],shape(0),0.0,env.gb(row-1)[col]);

               tmp4.clear();
               Permute(env.gb(row-1)[col],shape(0,2,3,1),tmp4);

               env.gb(row-1)[col] = std::move(tmp4);

            }
            else if(dir == VERTICAL){

               DArray<8> tmp8;
               Contract(1.0,RI8,shape(5),R_l[2],shape(0),0.0,tmp8);

               //and again
               RI8.clear();
               Contract(1.0,tmp8,shape(5),R_l[2],shape(0),0.0,RI8);

               tmp8.clear();
               Permute(RI8,shape(0,1,2,3,4,6,7,5),tmp8);

               RI8 = std::move(tmp8);

            }
            else{//diagonal lurd, change lower peps and LI8

               //lower peps
               DArray<5> tmp5;
               Contract(1.0,peps(lrow-1,lcol),shape(1),R_l[2],shape(0),0.0,tmp5);

               Permute(tmp5,shape(0,4,1,2,3),peps(lrow-1,lcol));

               //LI8
               DArray<8> tmp8;
               Contract(1.0,LI8,shape(3),R_l[2],shape(0),0.0,tmp8);

               //and again
               LI8.clear();
               Contract(1.0,tmp8,shape(3),R_l[2],shape(0),0.0,LI8);

               tmp8.clear();
               Permute(LI8,shape(0,1,2,6,7,3,4,5),tmp8);

               LI8 = std::move(tmp8);

            }

         }
         
         if(N_eff.shape(3) > 1 && dir != HORIZONTAL){//right

            Permute(X,shape(0,1,2,4,3),X_copy);

            //QR
            Geqrf(X_copy,R_l[3]);

            //add to right side of tensor
            DArray<5> tmp5;
            Contract(1.0,peps(lrow,lcol),shape(4),R_l[3],shape(1),0.0,tmp5);

            peps(lrow,lcol) = std::move(tmp5);

            //add  inverse to environment, for upper and lower layer
            invert(R_l[3]);

            if(dir == VERTICAL){

               DArray<8> tmp8;
               Contract(1.0,RI8,shape(3),R_l[3],shape(0),0.0,tmp8);

               //and again
               RI8.clear();
               Contract(1.0,tmp8,shape(3),R_l[3],shape(0),0.0,RI8);

               tmp8.clear();
               Permute(RI8,shape(0,1,2,6,7,3,4,5),tmp8);

               RI8 = std::move(tmp8);

            }
            else if(dir == DIAGONAL_LURD){

               //change right peps
               DArray<5> tmp5;
               Contract(1.0,R_l[3],shape(0),peps(lrow,lcol+1),shape(0),0.0,tmp5);

               peps(lrow,lcol+1) = std::move(tmp5);

               //change RI8
               DArray<8> tmp8;
               Contract(1.0,RI8,shape(1),R_l[3],shape(0),0.0,tmp8);

               //and again
               RI8.clear();
               Contract(1.0,tmp8,shape(1),R_l[3],shape(0),0.0,RI8);

               tmp8.clear();
               Permute(RI8,shape(0,6,7,1,2,3,4,5),tmp8);

               RI8 = std::move(tmp8);

            }
            else{//diagonal ldru

               //change peps
               DArray<5> tmp5;
               Contract(1.0,R_l[3],shape(0),peps(lrow,lcol+1),shape(0),0.0,tmp5);

               peps(lrow,lcol+1) = std::move(tmp5);

               //change RI8
               DArray<8> tmp8;
               Contract(1.0,RI8,shape(5),R_l[3],shape(0),0.0,tmp8);

               //and again
               RI8.clear();
               Contract(1.0,tmp8,shape(5),R_l[3],shape(0),0.0,RI8);

               tmp8.clear();
               Permute(RI8,shape(0,1,2,3,4,6,7,5),tmp8);

               RI8 = std::move(tmp8);

            }

         }

         // -----------------------------//
         // --- (B) ---- RIGHT SITE ---- //
         // -----------------------------//

         //calculate the environment for the right site
         N_eff.clear();
         calc_N_eff(dir,row,col,peps,N_eff,LO,RO,LI8,RI8,false);

         //eigenvalues
         eig.clear();
         diagonalize(N_eff,eig);

#ifdef _DEBUG
         cout << "(right)\t" << std::scientific << eig(eig.size() - 1) / eig(0) << endl;
#endif

         //are needed for the calculation of the positive approximant: physical dimension is last index (4)
         X.clear();
         get_X(N_eff,eig,X);

         if(N_eff.shape(0) > 1 && dir != HORIZONTAL){//left

            X_copy = X;

            //QR: watch out, LQ decomposition
            Gelqf(R_r[0],X_copy);

            //add to left side of tensor
            DArray<5> tmp5;
            Contract(1.0,R_r[0],shape(0),peps(rrow,rcol),shape(0),0.0,tmp5);

            peps(rrow,rcol) = std::move(tmp5);

            //add  inverse to environment, for upper and lower layer
            invert(R_r[0]);

            if(dir == VERTICAL){

               DArray<8> tmp8;
               Contract(1.0,LI8,shape(3),R_r[0],shape(1),0.0,tmp8);

               //and again
               LI8.clear();
               Contract(1.0,tmp8,shape(3),R_r[0],shape(1),0.0,LI8);

               tmp8.clear();
               Permute(LI8,shape(0,1,2,6,7,3,4,5),tmp8);

               LI8 = std::move(tmp8);

            }
            else{//diagonal lurd or ldru

               DArray<5> tmp5;
               Contract(1.0,peps(rrow,rcol-1),shape(4),R_r[0],shape(1),0.0,tmp5);

               peps(rrow,rcol-1) = std::move(tmp5);

            }

         }
         
         if(N_eff.shape(1) > 1){//up

            Permute(X,shape(0,2,3,4,1),X_copy);

            //QR
            Geqrf(X_copy,R_r[1]);

            //add to up side of tensor
            DArray<5> tmp5;
            Contract(1.0,peps(rrow,rcol),shape(1),R_r[1],shape(1),0.0,tmp5);

            Permute(tmp5,shape(0,4,1,2,3),peps(rrow,rcol));

            //add  inverse to environment, for upper and lower layer
            invert(R_r[1]);

            if(dir == HORIZONTAL){

               DArray<8> tmp8;
               Contract(1.0,RI8,shape(3),R_r[1],shape(0),0.0,tmp8);

               //and again
               RI8.clear();
               Contract(1.0,tmp8,shape(3),R_r[1],shape(0),0.0,RI8);

               tmp8.clear();
               Permute(RI8,shape(0,1,2,6,7,3,4,5),tmp8);

               RI8 = std::move(tmp8);

            }
            else if(dir == VERTICAL){

               DArray<8> tmp8;
               Contract(1.0,LI8,shape(1),R_r[1],shape(0),0.0,tmp8);

               //and again
               LI8.clear();
               Contract(1.0,tmp8,shape(1),R_r[1],shape(0),0.0,LI8);

               tmp8.clear();
               Permute(LI8,shape(0,6,7,1,2,3,4,5),tmp8);

               LI8 = std::move(tmp8);

            }
            else if(dir == DIAGONAL_LURD){//only change upper peps

               DArray<5> tmp5;
               Contract(1.0,peps(rrow+1,rcol),shape(3),R_r[1],shape(0),0.0,tmp5);

               Permute(tmp5,shape(0,1,2,4,3),peps(rrow+1,rcol));

            }
            else{//diagonal ldru: change upper environment

               DArray<4> tmp4;
               Contract(1.0,env.gt(row)[col+1],shape(1),R_r[1],shape(0),0.0,tmp4);

               //and again
               env.gt(row)[col+1].clear();
               Contract(1.0,tmp4,shape(1),R_r[1],shape(0),0.0,env.gt(row)[col+1]);

               tmp4.clear();
               Permute(env.gt(row)[col+1],shape(0,2,3,1),tmp4);

               env.gt(row)[col+1] = std::move(tmp4);

            }

         }
         
         if(N_eff.shape(2) > 1 && dir != VERTICAL){//down

            Permute(X,shape(0,1,3,4,2),X_copy);

            //QR
            Geqrf(X_copy,R_r[2]);

            //add to down side of tensor
            DArray<5> tmp5;
            Contract(1.0,peps(rrow,rcol),shape(3),R_r[2],shape(1),0.0,tmp5);

            Permute(tmp5,shape(0,1,2,4,3),peps(rrow,rcol));

            //add  inverse to environment, for upper and lower layer
            invert(R_r[2]);

            if(dir == HORIZONTAL){

               DArray<4> tmp4;
               Contract(1.0,env.gb(row-1)[col+1],shape(1),R_r[2],shape(0),0.0,tmp4);

               //and again
               env.gb(row-1)[col+1].clear();
               Contract(1.0,tmp4,shape(1),R_r[2],shape(0),0.0,env.gb(row-1)[col+1]);

               tmp4.clear();
               Permute(env.gb(row-1)[col+1],shape(0,2,3,1),tmp4);

               env.gb(row-1)[col+1] = std::move(tmp4);

            }

         }
         /*
            if(N_eff.shape(3) > 1){//right

            Permute(X,shape(0,1,2,4,3),X_copy);

         //QR
         Geqrf(X_copy,R_r[3]);

         //add to right side of tensor
         DArray<5> tmp5;
         Contract(1.0,peps(rrow,rcol),shape(4),R_r[3],shape(1),0.0,tmp5);

         peps(rrow,rcol) = std::move(tmp5);

         //add  inverse to environment, for upper and lower layer
         invert(R_r[3]);

         if(dir == HORIZONTAL){

         DArray<8> tmp8;
         Contract(1.0,RI8,shape(5),R_r[3],shape(0),0.0,tmp8);

         //and again
         RI8.clear();
         Contract(1.0,tmp8,shape(5),R_r[3],shape(0),0.0,RI8);

         Permute(RI8,shape(0,1,2,3,4,6,7,5),tmp8);

         RI8 = std::move(tmp8);

         }
         else if(dir == VERTICAL){

         DArray<8> tmp8;
         Contract(1.0,RI8,shape(1),R_r[3],shape(0),0.0,tmp8);

         //and again
         RI8.clear();
         Contract(1.0,tmp8,shape(1),R_r[3],shape(0),0.0,RI8);

         Permute(RI8,shape(0,6,7,1,2,3,4,5),tmp8);

         RI8 = std::move(tmp8);

         }

         }
         */
      }

   /**
    * diagonalize the effective environment:
    * @param N_eff is the effective environmnt: output eigenvectors
    * @param eig output eigenvalues
    */
   void diagonalize(DArray<8> &N_eff,DArray<1> &eig){

      int n = N_eff.shape(0) * N_eff.shape(1) * N_eff.shape(2) * N_eff.shape(3);

      for(int i = 0;i < n;++i)
         for(int j = i + 1;j < n;++j){

            N_eff.data()[i + n*j] = 0.5 * ( N_eff.data()[i + n*j] + N_eff.data()[j + n*i] );
            N_eff.data()[j + n*i] = N_eff.data()[i + n*j];

         }

      eig.resize(n);
      lapack::syev(CblasRowMajor,'V','U',n,N_eff.data(),n,eig.data());

   }


   /**
    * get the X - peps: X^T X ~ N_eff positive approximant of environment
    * @param N_eff is the effective environmnt: output eigenvectors
    * @param eig output eigenvalues
    */
   void get_X(const DArray<8> &N_eff,const DArray<1> &eig,DArray<5> &X){

      int n = N_eff.shape(0) * N_eff.shape(1) * N_eff.shape(2) * N_eff.shape(3);

      X.resize( shape(N_eff.shape(0),N_eff.shape(1),N_eff.shape(2),N_eff.shape(3),n) );

      //get the square root of the positive approximant:
      for(int i = 0;i < n;++i)
         for(int j = 0;j < n;++j)
            if(eig(j) > 0.0)
               X.data()[i*n + j] = sqrt( eig(j) ) * N_eff.data()[i*n + j];

   }

   /** 
    * wrapper function invert square general matrix DArray<2>.
    * @param A both input as output matrix: on input A, on output A^{-1}
    */
   void invert(DArray<2> &A){

      int *ipiv = new int [A.shape(0)];

      lapack::getrf(CblasRowMajor,A.shape(0),A.shape(1), A.data(), A.shape(1), ipiv);

      lapack::getri(CblasRowMajor,A.shape(0), A.data(), A.shape(1), ipiv);

      delete [] ipiv;

   }

   /** 
    * restore the tensors to their original state, i.e. undo canonicalization
    * @param dir vertical, horizontal,diagonal lurd or diagonal ldru update
    * @param row the row index of the bottom site
    * @param col column index of the vertical column
    * @param peps, full PEPS object before update
    * @param L left environment
    * @param R right environment
    * @param R_l vector of DArray<2> objects containing the inverse R of the QR decompostion (left site)
    * @param R_r vector of DArray<2> objects containing the inverse R of the QR decompostion (right site)
    */
   template<>
      void restore(const PROP_DIR &dir,int row,int col,PEPS<double> &peps, DArray<5> &L,DArray<5> &R,

            std::vector< DArray<2> > &R_l, std::vector< DArray<2> > &R_r){

         //set left and right indices
         int lrow = row;
         int rrow = row;

         int lcol = col;
         int rcol = col;

         if(dir == VERTICAL)
            rrow++;
         else if(dir == HORIZONTAL)
            rcol++;
         else if(dir == DIAGONAL_LURD){

            lrow++;
            rcol++;

         }
         else{//DIAGONAL_LDRU

            rrow++;
            rcol++;

         }

         //--- LEFT SITE ---

         if(peps(lrow,lcol).shape(0) > 1){//left

            //add to left side of tensor
            DArray<5> tmp5;
            Contract(1.0,R_l[0],shape(0),peps(lrow,lcol),shape(0),0.0,tmp5);

            peps(lrow,lcol) = std::move(tmp5);

            if(dir == DIAGONAL_LURD){//diagonal needs more restoring

               //invert back to original
               invert(R_l[0]);

               //add to L
               tmp5.clear();
               Contract(1.0,L,shape(1),R_l[0],shape(1),0.0,tmp5);

               //and again
               L.clear();
               Contract(1.0,tmp5,shape(1),R_l[0],shape(1),0.0,L);

               Permute(L,shape(0,3,4,1,2),tmp5);

               L = std::move(tmp5);

            }
            else if(dir == DIAGONAL_LDRU){

               //invert back to original
               invert(R_l[0]);

               tmp5.clear();
               Contract(1.0,L,shape(3),R_l[0],shape(1),0.0,tmp5);

               //and again
               L.clear();
               Contract(1.0,tmp5,shape(3),R_l[0],shape(1),0.0,L);

            }

         }

         if(peps(rrow,rcol).shape(1) > 1 && dir != VERTICAL ){//up

            //add to up side of tensor
            DArray<5> tmp5;
            Contract(1.0,peps(lrow,lcol),shape(1),R_l[1],shape(1),0.0,tmp5);

            Permute(tmp5,shape(0,4,1,2,3),peps(lrow,lcol));

            if(dir == HORIZONTAL){

               if(row == Ly - 2){

                  //restore the upper tensor
                  invert(R_l[1]);

                  DArray<5> tmp5;
                  Contract(1.0,peps(row+1,col),shape(3),R_l[1],shape(0),0.0,tmp5);

                  Permute(tmp5,shape(0,1,2,4,3),peps(row+1,col));

               }

            }
            else if(dir == DIAGONAL_LURD){//diagonal needs more restoring

               if(row == 0){

                  invert(R_l[1]);

                  DArray<4> tmp4;
                  Contract(1.0,env.gt(0)[lcol],shape(1),R_l[1],shape(0),0.0,tmp4);

                  env.gt(0)[lcol].clear();
                  Contract(1.0,tmp4,shape(1),R_l[1],shape(0),0.0,env.gt(0)[lcol]);

                  Permute(env.gt(0)[lcol],shape(0,2,3,1),tmp4);

                  env.gt(0)[lcol] = std::move(tmp4);

               }
               else{//other stuff, later

               }

            }
            else if(dir == DIAGONAL_LDRU){

               invert(R_l[1]);

               DArray<5> tmp5;
               Contract(1.0,peps(lrow+1,lcol),shape(3),R_l[1],shape(0),0.0,tmp5);

               Permute(tmp5,shape(0,1,2,4,3),peps(lrow+1,lcol));

            }

         }

         if(peps(lrow,lcol).shape(3) > 1){//down

            //add to up side of tensor
            DArray<5> tmp5;
            Contract(1.0,peps(lrow,lcol),shape(3),R_l[2],shape(1),0.0,tmp5);

            Permute(tmp5,shape(0,1,2,4,3),peps(lrow,lcol));

            if(dir == DIAGONAL_LURD){//some more restoration

               invert(R_l[2]);

               DArray<5> tmp5;
               Contract(1.0,peps(lrow-1,lcol),shape(1),R_l[2],shape(0),0.0,tmp5);

               Permute(tmp5,shape(0,4,1,2,3),peps(lrow-1,lcol));

            }

         }

         if(peps(lrow,lcol).shape(4) > 1 && dir != HORIZONTAL){//right

            //add to right side of tensor
            DArray<5> tmp5;
            Contract(1.0,peps(lrow,lcol),shape(4),R_l[3],shape(1),0.0,tmp5);

            peps(lrow,lcol) = std::move(tmp5);

            if(dir == DIAGONAL_LDRU || dir == DIAGONAL_LURD){//some more restoration

               invert(R_l[3]);

               DArray<5> tmp5;
               Contract(1.0,R_l[3],shape(0),peps(lrow,lcol+1),shape(0),0.0,tmp5);

               peps(lrow,lcol+1) = std::move(tmp5);

            }

         }

         //RIGHT SITE
         if(peps(rrow,rcol).shape(0) > 1 && dir != HORIZONTAL){//left

            //add to right side of tensor
            DArray<5> tmp5;
            Contract(1.0,R_r[0],shape(0),peps(rrow,rcol),shape(0),0.0,tmp5);

            peps(rrow,rcol) = std::move(tmp5);

            if(dir == DIAGONAL_LDRU || dir == DIAGONAL_LURD){//some more restoration

               invert(R_r[0]);

               DArray<5> tmp5;
               Contract(1.0,peps(rrow,rcol-1),shape(4),R_r[0],shape(1),0.0,tmp5);

               peps(rrow,rcol-1) = std::move(tmp5);

            }

         }

         if(peps(rrow,rcol).shape(1) > 1){//up

            //add to up side of tensor
            DArray<5> tmp5;
            Contract(1.0,peps(rrow,rcol),shape(1),R_r[1],shape(1),0.0,tmp5);

            Permute(tmp5,shape(0,4,1,2,3),peps(rrow,rcol));

            if(dir == HORIZONTAL){

               if(row == Ly - 2){

                  //restore the upper tensor
                  invert(R_r[1]);

                  DArray<5> tmp5;
                  Contract(1.0,peps(row+1,col+1),shape(3),R_r[1],shape(0),0.0,tmp5);

                  Permute(tmp5,shape(0,1,2,4,3),peps(row+1,col+1));

               }

            }
            else if(dir == DIAGONAL_LURD){

               invert(R_r[1]);

               DArray<5> tmp5;
               Contract(1.0,peps(rrow+1,rcol),shape(3),R_r[1],shape(0),0.0,tmp5);

               Permute(tmp5,shape(0,1,2,4,3),peps(rrow+1,rcol));

            }
            else if(dir == DIAGONAL_LDRU){

               if(row == 0){

                  invert(R_r[1]);

                  DArray<4> tmp4;
                  Contract(1.0,env.gt(0)[rcol],shape(1),R_r[1],shape(0),0.0,tmp4);

                  env.gt(0)[rcol].clear();
                  Contract(1.0,tmp4,shape(1),R_r[1],shape(0),0.0,env.gt(0)[rcol]);

                  Permute(env.gt(0)[rcol],shape(0,2,3,1),tmp4);

                  env.gt(0)[rcol] = std::move(tmp4);

               }
               else{//later

               }

            }

         }

         if(peps(rrow,rcol).shape(3) > 1 && dir != VERTICAL){//down

            //add to up side of tensor
            DArray<5> tmp5;
            Contract(1.0,peps(rrow,rcol),shape(3),R_r[2],shape(1),0.0,tmp5);

            Permute(tmp5,shape(0,1,2,4,3),peps(rrow,rcol));

            if(dir == DIAGONAL_LURD){

               invert(R_r[2]);

               DArray<5> tmp5;
               Contract(1.0,peps(rrow-1,rcol),shape(1),R_r[2],shape(0),0.0,tmp5);

               Permute(tmp5,shape(0,4,1,2,3),peps(rrow-1,rcol));

            }

         }

         if(peps(rrow,rcol).shape(4) > 1){//right

            //add to right side of tensor
            DArray<5> tmp5;
            Contract(1.0,peps(rrow,rcol),shape(4),R_r[3],shape(1),0.0,tmp5);

            peps(rrow,rcol) = std::move(tmp5);

            //set back the R tensors
            if(dir == DIAGONAL_LURD){

               invert(R_r[3]);

               DArray<5> tmp5;
               Contract(1.0,R,shape(3),R_r[3],shape(0),0.0,tmp5);

               //and again
               R.clear();
               Contract(1.0,tmp5,shape(3),R_r[3],shape(0),0.0,R);

            }
            else if(dir == DIAGONAL_LDRU){

               invert(R_r[3]);

               DArray<5> tmp5;
               Contract(1.0,R,shape(1),R_r[3],shape(0),0.0,tmp5);

               //and again
               R.clear();
               Contract(1.0,tmp5,shape(1),R_r[3],shape(0),0.0,R);

               tmp5.clear();
               Permute(R,shape(0,3,4,1,2),tmp5);

               R = std::move(tmp5);

            }

         }

      }


   /** 
    * restore the tensors to their original state, i.e. undo canonicalization
    * @param dir vertical, horizontal,diagonal lurd or diagonal ldru update
    * @param row the row index of the bottom site
    * @param col column index of the vertical column
    * @param LO left environment
    * @param RO right environment
    * @param peps, full PEPS object before update
    * @param R_l vector of DArray<2> objects containing the inverse R of the QR decompostion (left site)
    * @param R_r vector of DArray<2> objects containing the inverse R of the QR decompostion (right site)
    */
   template<>
      void restore(const PROP_DIR &dir,int row,int col,PEPS<double> &peps, DArray<6> &LO,DArray<6> &RO,

            std::vector< DArray<2> > &R_l, std::vector< DArray<2> > &R_r){

         //set left and right indices
         int lrow = row;
         int rrow = row;

         int lcol = col;
         int rcol = col;

         if(dir == VERTICAL)
            rrow++;
         else if(dir == HORIZONTAL)
            rcol++;
         else if(dir == DIAGONAL_LURD){

            lrow++;
            rcol++;

         }
         else{//DIAGONAL_LDRU

            rrow++;
            rcol++;

         }

         //--- LEFT SITE ---

         if(peps(lrow,lcol).shape(0) > 1){//left

            //add to left side of tensor
            DArray<5> tmp5;
            Contract(1.0,R_l[0],shape(0),peps(lrow,lcol),shape(0),0.0,tmp5);

            peps(lrow,lcol) = std::move(tmp5);

            if(dir == DIAGONAL_LDRU){//diagonal needs more restoring

               //invert back to original
               invert(R_l[0]);

               //LO
               DArray<6> tmp6;
               Contract(1.0,LO,shape(1),R_l[0],shape(1),0.0,tmp6);

               //and again
               LO.clear();
               Contract(1.0,tmp6,shape(1),R_l[0],shape(1),0.0,LO);

               Permute(LO,shape(0,4,5,1,2,3),tmp6);

               LO = std::move(tmp6);

            }
            else if(dir == DIAGONAL_LDRU){

               //invert back to original
               invert(R_l[0]);

               //LO
               DArray<6> tmp6;
               Contract(1.0,LO,shape(3),R_l[0],shape(1),0.0,tmp6);

               //and again
               LO.clear();
               Contract(1.0,tmp6,shape(3),R_l[0],shape(1),0.0,LO);

               Permute(LO,shape(0,1,2,4,5,3),tmp6);

               LO = std::move(tmp6);

            }

         }

         if(peps(rrow,rcol).shape(1) > 1 && dir != VERTICAL ){//up

            //add to up side of tensor
            DArray<5> tmp5;
            Contract(1.0,peps(lrow,lcol),shape(1),R_l[1],shape(1),0.0,tmp5);

            Permute(tmp5,shape(0,4,1,2,3),peps(lrow,lcol));

            if(dir == DIAGONAL_LURD){

               //invert back to original
               invert(R_l[1]);

               DArray<4> tmp4;
               Contract(1.0,env.gt(row)[col],shape(1),R_l[1],shape(0),0.0,tmp4);

               //and again
               env.gt(row)[col].clear();
               Contract(1.0,tmp4,shape(1),R_l[1],shape(0),0.0,env.gt(row)[col]);

               tmp4.clear();
               Permute(env.gt(row)[col],shape(0,2,3,1),tmp4);

               env.gt(row)[col] = std::move(tmp4);

            }
            else if(dir == DIAGONAL_LDRU){

               //invert back to original
               invert(R_l[1]);

               //add to up side of tensor
               DArray<5> tmp5;
               Contract(1.0,peps(lrow+1,lcol),shape(3),R_l[1],shape(1),0.0,tmp5);

               Permute(tmp5,shape(0,1,2,4,3),peps(lrow,lcol));

            }

         }

         if(peps(lrow,lcol).shape(3) > 1){//down

            //add to up side of tensor
            DArray<5> tmp5;
            Contract(1.0,peps(lrow,lcol),shape(3),R_l[2],shape(1),0.0,tmp5);

            Permute(tmp5,shape(0,1,2,4,3),peps(lrow,lcol));

            if(dir == HORIZONTAL || dir == DIAGONAL_LDRU){//possibly extra term

               invert(R_l[2]);

               DArray<4> tmp4;
               Contract(1.0,env.gb(row-1)[col],shape(1),R_l[2],shape(0),0.0,tmp4);

               //and again
               env.gb(row-1)[col].clear();
               Contract(1.0,tmp4,shape(1),R_l[2],shape(0),0.0,env.gb(row-1)[col]);

               tmp4.clear();
               Permute(env.gb(row-1)[col],shape(0,2,3,1),tmp4);

               env.gb(row-1)[col] = std::move(tmp4);

            }
            else{//diagonal lurd, some more restoration

               invert(R_l[2]);

               DArray<5> tmp5;
               Contract(1.0,peps(lrow-1,lcol),shape(1),R_l[2],shape(0),0.0,tmp5);

               Permute(tmp5,shape(0,4,1,2,3),peps(lrow-1,lcol));

            }

         }

         if(peps(lrow,lcol).shape(4) > 1 && dir != HORIZONTAL){//right

            //add to right side of tensor
            DArray<5> tmp5;
            Contract(1.0,peps(lrow,lcol),shape(4),R_l[3],shape(1),0.0,tmp5);

            peps(lrow,lcol) = std::move(tmp5);

            if(dir == DIAGONAL_LDRU || dir == DIAGONAL_LURD){//some more restoration

               invert(R_l[3]);

               DArray<5> tmp5;
               Contract(1.0,R_l[3],shape(0),peps(lrow,lcol+1),shape(0),0.0,tmp5);

               peps(lrow,lcol+1) = std::move(tmp5);

            }

         }

         //RIGHT SITE
         if(peps(rrow,rcol).shape(0) > 1 && dir != HORIZONTAL){//left

            //add to right side of tensor
            DArray<5> tmp5;
            Contract(1.0,R_r[0],shape(0),peps(rrow,rcol),shape(0),0.0,tmp5);

            peps(rrow,rcol) = std::move(tmp5);

            if(dir == DIAGONAL_LDRU || dir == DIAGONAL_LURD){//some more restoration

               invert(R_r[0]);

               DArray<5> tmp5;
               Contract(1.0,peps(rrow,rcol-1),shape(4),R_r[0],shape(1),0.0,tmp5);

               peps(rrow,rcol-1) = std::move(tmp5);

            }

         }

         if(peps(rrow,rcol).shape(1) > 1){//up

            //add to up side of tensor
            DArray<5> tmp5;
            Contract(1.0,peps(rrow,rcol),shape(1),R_r[1],shape(1),0.0,tmp5);

            Permute(tmp5,shape(0,4,1,2,3),peps(rrow,rcol));

            if(dir == DIAGONAL_LURD){

               //invert back to original
               invert(R_r[1]);

               DArray<5> tmp5;
               Contract(1.0,peps(rrow+1,rcol),shape(3),R_r[1],shape(0),0.0,tmp5);

               Permute(tmp5,shape(0,1,2,4,3),peps(rrow+1,rcol));

            }
            else if(dir == DIAGONAL_LDRU){

               //invert back to original
               invert(R_r[1]);

               DArray<4> tmp4;
               Contract(1.0,env.gt(row)[col+1],shape(1),R_r[1],shape(0),0.0,tmp4);

               //and again
               env.gt(row)[col+1].clear();
               Contract(1.0,tmp4,shape(1),R_r[1],shape(0),0.0,env.gt(row)[col+1]);

               tmp4.clear();
               Permute(env.gt(row)[col+1],shape(0,2,3,1),tmp4);

               env.gt(row)[col+1] = std::move(tmp4);

            }

         }
         
         if(peps(rrow,rcol).shape(3) > 1 && dir != VERTICAL){//down

            //add to up side of tensor
            DArray<5> tmp5;
            Contract(1.0,peps(rrow,rcol),shape(3),R_r[2],shape(1),0.0,tmp5);

            Permute(tmp5,shape(0,1,2,4,3),peps(rrow,rcol));

            if(dir == HORIZONTAL){//possibly extra term

               invert(R_r[2]);

               DArray<4> tmp4;
               Contract(1.0,env.gb(row-1)[col+1],shape(1),R_r[2],shape(0),0.0,tmp4);

               //and again
               env.gb(row-1)[col+1].clear();
               Contract(1.0,tmp4,shape(1),R_r[2],shape(0),0.0,env.gb(row-1)[col+1]);

               tmp4.clear();
               Permute(env.gb(row-1)[col+1],shape(0,2,3,1),tmp4);

               env.gb(row-1)[col+1] = std::move(tmp4);

            }
            else if(dir == DIAGONAL_LURD || dir == DIAGONAL_LDRU){

               invert(R_r[2]);

               DArray<5> tmp5;
               Contract(1.0,peps(rrow-1,rcol),shape(1),R_r[2],shape(0),0.0,tmp5);

               Permute(tmp5,shape(0,4,1,2,3),peps(rrow-1,rcol));

            }

         }
         /*
            if(peps(rrow,rcol).shape(4) > 1){//right

//add to right side of tensor
DArray<5> tmp5;
Contract(1.0,peps(rrow,rcol),shape(4),R_r[3],shape(1),0.0,tmp5);

peps(rrow,rcol) = std::move(tmp5);

//set back the R tensors
if(dir == DIAGONAL_LURD){

}
else if(dir == DIAGONAL_LDRU){

}

}
*/
}

/**
 * shift the singular values to the next site in the sweep, i.e. do QR and past R to the right
 * @param option shift to left or shift to right
 * @param row index
 * @param col index
 * @param peps object containing all the tensors
 */
void shift_col(char option,int row,int col,PEPS<double> &peps){

   if(option == 'r'){//shift to the right

      //QR
      DArray<2> tmp2;
      Geqrf(peps(row,col),tmp2);

      //add to left side of next tensor
      DArray<5> tmp5;
      Contract(1.0,tmp2,shape(1),peps(row,col+1),shape(0),0.0,tmp5);

      peps(row,col+1) = std::move(tmp5);

   }
   else{//shift to the left

      //LQ
      DArray<2> tmp2;
      Gelqf(tmp2,peps(row,col));

      //add to left side of next tensor
      DArray<5> tmp5;
      Contract(1.0,peps(row,col-1),shape(4),tmp2,shape(0),0.0,tmp5);

      peps(row,col-1) = std::move(tmp5);

   }

}

/**
 * shift the singular values to the next row in the sweep, i.e. do QR and past R to the upper row
 * @param option top or bottom shift?
 * @param row what row are you on
 * @param peps object containing the tensors
 */
void shift_row(char option,int row,PEPS<double> &peps){

   if(option == 'b'){

      DArray<2> tmp2;

      for(int col = 0;col < Ly;++col){

         DArray<5> tmp5;
         Permute(peps(row,col),shape(0,2,3,4,1),tmp5);

         Geqrf(tmp5,tmp2);

         Permute(tmp5,shape(0,4,1,2,3),peps(row,col));

         //add to down side of upper tensor
         tmp5.clear();
         Contract(1.0,tmp2,shape(1),peps(row+1,col),shape(3),0.0,tmp5);

         Permute(tmp5,shape(1,2,3,0,4),peps(row+1,col));

      }

   }
   else{//top: shift the QR downwards

      DArray<2> tmp2;

      for(int col = 0;col < Ly;++col){

         DArray<5> tmp5;
         Permute(peps(row,col),shape(0,1,2,4,3),tmp5);

         Geqrf(tmp5,tmp2);

         Permute(tmp5,shape(0,1,2,4,3),peps(row,col));

         //add to up side of lower tensor
         tmp5.clear();
         Contract(1.0,tmp2,shape(1),peps(row-1,col),shape(1),0.0,tmp5);

         Permute(tmp5,shape(1,0,2,3,4),peps(row-1,col));


      }

   }

}

/**
 * regularize the effective environment by adding a scaled unit matrix to it
 * @param N_eff the environemnt, changed on exit
 * @param num constant by which to mupliply the unit matrix
 */
void regularize(DArray<8> &N_eff,double num){

   int dim = N_eff.shape(0) * N_eff.shape(1) * N_eff.shape(2) * N_eff.shape(3);

   for(int i = 0;i < dim;++i)
      N_eff.data()[i + dim*i] += num;

}

}
