#include <iostream>
#include <iomanip>
#include <fstream>
#include <complex>

using std::cout;
using std::endl;
using std::ostream;
using std::complex;
using std::ofstream;

#include "include.h"

using namespace btas;
using namespace global;

namespace contractions {

  /**
    * update left renormalized operator 
    * @param option 'b'ottom or 't'op
    * @param col is column index
    * @param peps is the PEPS object
    * @param L is the left environment to be updated
    */
   void update_L(char option,int col,const PEPS<double> &peps,DArray<5> &L){

      if(option == 'b'){//row == 0

         DArray<7> tmp7;
         Gemm(CblasTrans,CblasNoTrans,1.0,L,env.gt(0)[col],0.0,tmp7);

         DArray<8> tmp8;
         Contract(1.0,tmp7,shape(0,4),peps(1,col),shape(0,1),0.0,tmp8);

         tmp7.clear();
         Contract(1.0,tmp8,shape(0,3,5),peps(1,col),shape(0,1,2),0.0,tmp7);

         tmp8.clear();
         Contract(1.0,tmp7,shape(0,3),peps(0,col),shape(0,1),0.0,tmp8);

         //and another peps
         L.clear();
         Contract(1.0,tmp8,shape(0,3,5,6),peps(0,col),shape(0,1,2,3),0.0,L);

      }
      else{//top: bottom peps row = Ly - 2

         if(col == 0){

            //first bottom peps to bottom environemnt
            DArray<5> tmp5;
            Contract(1.0,peps(Ly-2,col),shape(0,3),env.gb(Ly - 3)[col],shape(0,2),0.0,tmp5);

            //add another peps to it
            DArray<6> tmp6;
            Contract(1.0,peps(Ly-2,col),shape(2,3),tmp5,shape(1,3),0.0,tmp6);

            //add top peps to it
            DArray<7> tmp7;
            Contract(1.0,peps(Ly-1,col),shape(0,3),tmp6,shape(0,3),0.0,tmp7);

            //one last top peps
            tmp6.clear();
            Contract(1.0,peps(Ly-1,col),shape(1,2,3),tmp7,shape(0,1,3),0.0,tmp6);

            L = tmp6.reshape_clear( shape(D,D,D,D,env.gb(Ly-3)[col].shape(3)) );

         }
         else{//if col != 0

            //add top peps to L
            DArray<8> tmp8;
            Contract(1.0,L,shape(0),peps(Ly-1,col),shape(0),0.0,tmp8);

            //and another top peps
            DArray<7> tmp7;
            Contract(1.0,tmp8,shape(0,4,5),peps(Ly-1,col),shape(0,1,2),0.0,tmp7);

            //now a bottom peps
            tmp8.clear();
            Contract(1.0,tmp7,shape(0,3),peps(Ly-2,col),shape(0,1),0.0,tmp8);

            tmp7.clear();
            Contract(1.0,tmp8,shape(0,3,5),peps(Ly-2,col),shape(0,1,2),0.0,tmp7);

            //finally environment:
            L.clear();
            Contract(1.0,tmp7,shape(0,3,5),env.gb(Ly - 3)[col],shape(0,1,2),0.0,L);

         }

      }

   }

   /** 
    * init the right renormalized operator for the middle rows: 
    * @param row index of the lowest row
    * @param peps The PEPS object
    * @param R vector containing the right operators on exit
    */
   void init_ro(int row,const PEPS<double> &peps,vector< DArray<6> > &RO){

      DArray<5> tmp5;
      DArray<5> tmp5bis;

      DArray<7> tmp7;
      DArray<7> tmp7bis;

      DArray<8> tmp8;
      DArray<8> tmp8bis;

      DArray<9> tmp9;
      DArray<9> tmp9bis;

      //first rightmost site: attach peps Lx - 1 to bottom
      Gemm(CblasNoTrans,CblasTrans,1.0,peps(row,Lx-1),env.gb(row-1)[Lx-1],0.0,tmp5);

      Permute(tmp5,shape(2,4,0,3,1),tmp5bis);

      //add second lowest peps
      int M = peps(row,Lx-1).shape(0) * peps(row,Lx-1).shape(1);
      int N = tmp5bis.shape(2) * tmp5bis.shape(3) * tmp5bis.shape(4);
      int K =  peps(row,Lx-1).shape(2) * peps(row,Lx-1).shape(3);

      tmp5.resize(shape(D,D,D,env.gb(row-1)[Lx-1].shape(0),D));
      blas::gemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, M, N, K, 1.0,peps(row,Lx-1).data(),K,tmp5bis.data(),N,0.0,tmp5.data(),N);

      //add one row up peps
      M = peps(row+1,Lx-1).shape(0) * peps(row+1,Lx-1).shape(1) * peps(row+1,Lx-1).shape(2);
      N = env.gb(row-1)[Lx-1].shape(0) * D * D * D;
      K = peps(row+1,Lx-1).shape(3);

      tmp7.resize( shape(D,D,d,D,D,D,env.gb(row-1)[Lx-1].shape(0)) );
      blas::gemm(CblasRowMajor, CblasNoTrans, CblasTrans, M, N, K, 1.0,peps(row+1,Lx-1).data(),K,tmp5.data(),K,0.0,tmp7.data(),N);

      Permute(tmp7,shape(2,4,0,1,3,5,6),tmp7bis);

      //and another
      M = peps(row+1,Lx-1).shape(0) * peps(row+1,Lx-1).shape(1);
      N = env.gb(row-1)[Lx-1].shape(0) * D * D * D * D;
      K = peps(row+1,Lx-1).shape(2) * peps(row+1,Lx-1).shape(3);

      tmp7.resize( shape(D,D,D,D,D,D,env.gb(row-1)[Lx-1].shape(0)) );
      blas::gemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, M, N, K, 1.0,peps(row+1,Lx-1).data(),K,tmp7bis.data(),N,0.0,tmp7.data(),N);

      tmp7bis.clear();
      Permute(tmp7,shape(1,3,0,2,4,5,6),tmp7bis);

      //finally top environment
      M = env.gt(row)[Lx-1].shape(0);
      N = env.gb(row-1)[Lx-1].shape(0) * D * D * D * D;
      K = env.gt(row)[Lx-1].shape(1) * env.gt(row)[Lx-1].shape(2);

      RO[Lx - 2].resize(env.gt(row)[Lx-1].shape(0),D,D,D,D,env.gb(row-1)[Lx-1].shape(0));
      blas::gemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, M, N, K, 1.0,env.gt(row)[Lx-1].data(),K,tmp7bis.data(),N,0.0,RO[Lx - 2].data(),N);

      //now move from right to left, constructing the rest
      for(int col = Lx - 2;col > 0;--col){

         //first add bottom to right unity
         tmp8.clear();
         Gemm(CblasNoTrans,CblasTrans,1.0,env.gb(row-1)[col],RO[col],0.0,tmp8);

         tmp8bis.clear();
         Permute(tmp8,shape(2,7,0,1,3,4,5,6),tmp8bis);

         //add regular peps on lower site
         tmp9.clear();
         Gemm(CblasNoTrans,CblasNoTrans,1.0,peps(row,col),tmp8bis,0.0,tmp9);

         tmp9bis.clear();
         Permute(tmp9,shape(2,4,8,0,1,3,5,6,7),tmp9bis);

         //and another regular peps on lower site
         tmp8.clear();
         Gemm(CblasNoTrans,CblasNoTrans,1.0,peps(row,col),tmp9bis,0.0,tmp8);

         tmp8bis.clear();
         Permute(tmp8,shape(3,7,1,6,0,2,4,5),tmp8bis);

         //add regular peps on upper site
         tmp9.clear();
         Gemm(CblasNoTrans,CblasNoTrans,1.0,peps(row+1,col),tmp8bis,0.0,tmp9);

         tmp9bis.clear();
         Permute(tmp9,shape(2,3,4,0,1,5,6,7,8),tmp9bis);

         //yet another regular on upper site
         tmp8.clear();
         Gemm(CblasNoTrans,CblasNoTrans,1.0,peps(row+1,col),tmp9bis,0.0,tmp8);

         tmp8bis.clear();
         Permute(tmp8,shape(1,3,7,0,2,4,5,6),tmp8bis);

         //finally top environment for closure
         RO[col - 1].clear();
         Gemm(CblasNoTrans,CblasNoTrans,1.0,env.gt(row)[col],tmp8bis,0.0,RO[col-1]);

      }

   }

   /** 
    * init the right renormalized operator for the top or bottom row
    * @param option == 't'op or 'b'ottom
    * @param R vector containing the right operators on exit
    */
   void init_ro(char option,const PEPS<double> &peps,vector< DArray<5> > &R){

      if(option == 'b'){

         R[Lx-1].resize( shape(1,1,1,1,1) );
         R[Lx - 1] = 1.0;

         int M,N,K;

         DArray<7> tmp7,tmp7bis;
         DArray<8> tmp8,tmp8bis;

         for(int i = Lx - 1;i > 0;--i){

            tmp7.clear();
            Contract(1.0,env.gt(0)[i],shape(3),R[i],shape(0),0.0,tmp7);

            tmp8.clear();
            Contract(1.0,tmp7,shape(1,3),peps(1,i),shape(1,4),0.0,tmp8);

            tmp7.clear();
            Contract(1.0,tmp8,shape(1,6,2),peps(1,i),shape(1,2,4),0.0,tmp7);

            tmp8.clear();
            Contract(1.0,tmp7,shape(4,1),peps(0,i),shape(1,4),0.0,tmp8);

            R[i - 1].clear();
            Contract(1.0,tmp8,shape(4,6,7,1),peps(0,i),shape(1,2,3,4),0.0,R[i-1]);

         }

      }
      else{ //top 2 rows

         //attach bottom environment to lower peps
         DArray<5> tmp5;
         Gemm(CblasNoTrans,CblasTrans,1.0,peps(Ly-2,Lx-1),env.gb(Ly-3)[Lx-1],0.0,tmp5);

         DArray<5> tmp5bis;
         Permute(tmp5,shape(2,4,0,1,3),tmp5bis);

         //another lower regular peps on top
         int M = peps(Ly-2,Lx-1).shape(0) * peps(Ly-2,Lx-1).shape(1);
         int N = tmp5bis.shape(2) * tmp5bis.shape(3) * tmp5bis.shape(4);
         int K = peps(Ly-2,Lx-1).shape(2) * peps(Ly-2,Lx-1).shape(3);

         tmp5.resize(shape( D,D,D,D,env.gb(Ly-3)[Lx-1].shape(0) ) );
         blas::gemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, M, N, K, 1.0,peps(Ly-2,Lx-1).data(),K,tmp5bis.data(),N,0.0,tmp5.data(),N);

         tmp5bis.clear();
         Permute(tmp5,shape(1,3,0,2,4),tmp5bis);

         //finally add top!
         M = env.gt(Ly-3)[Lx-1].shape(0);
         N = tmp5bis.shape(2) * tmp5bis.shape(3) * tmp5bis.shape(4);
         K = env.gt(Ly-3)[Lx-1].shape(1) * env.gt(Ly-3)[Lx-1].shape(2) * env.gt(Ly-3)[Lx-1].shape(3);

         R[Lx-2].resize(shape( D,D,D,D,env.gb(Ly-3)[Lx-1].shape(0) ) );
         blas::gemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, M, N, K, 1.0,env.gt(Ly-3)[Lx-1].data(),K,tmp5bis.data(),N,0.0,R[Lx-2].data(),N);

         DArray<7> tmp7;
         DArray<7> tmp7bis;

         DArray<8> tmp8;
         DArray<8> tmp8bis;

         for(int col = Lx - 2;col > 0;--col){

            tmp7.clear();
            Gemm(CblasNoTrans,CblasTrans,1.0,env.gb(Ly-3)[col],R[col],0.0,tmp7);

            tmp7bis.clear();
            Permute(tmp7,shape(2,6,0,1,3,4,5),tmp7bis);

            tmp8.clear();
            Gemm(CblasNoTrans,CblasNoTrans,1.0,peps(Ly-2,col),tmp7bis,0.0,tmp8);

            tmp8bis.clear();
            Permute(tmp8,shape(2,4,7,0,1,3,5,6),tmp8bis);

            tmp7.clear();
            Gemm(CblasNoTrans,CblasNoTrans,1.0,peps(Ly-2,col),tmp8bis,0.0,tmp7);

            tmp7bis.clear();
            Permute(tmp7,shape(1,3,5,6,0,2,4),tmp7bis);

            M = env.gt(Ly-3)[col].shape(0);
            N = tmp7bis.shape(4) * tmp7bis.shape(5) * tmp7bis.shape(6);
            K = env.gt(Ly-3)[col].shape(1) * env.gt(Ly-3)[col].shape(2) * env.gt(Ly-3)[col].shape(3);

            R[col - 1].resize(shape( D,D,D,D,env.gb(Ly-3)[col].shape(0) ) );
            blas::gemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, M, N, K, 1.0,env.gt(Ly-3)[col].data(),K,tmp7bis.data(),N,0.0,R[col - 1].data(),N);

         }

      }

   }

   /**
    * update left renormalized operator on site (row,col ) for middle sites
    * @param row index of bottom row peps
    * @param col index of leftmost column
    * @param peps the input PEPS object
    * @param LO input old left renormalized operator, output new left renormalized operator
    */
   void update_L(int row,int col,const PEPS<double> &peps,DArray<6> &LO){

      if(col == 0){

         //paste top environment on
         DArray<5> tmp5;
         Gemm(CblasTrans,CblasNoTrans,1.0,env.gt(row)[0],peps(row+1,0),0.0,tmp5);

         DArray<6> tmp6;
         Contract(1.0,tmp5,shape(0,2),peps(row+1,0),shape(1,2),0.0,tmp6);

         DArray<7> tmp7;
         Contract(1.0,tmp6,shape(3,1),peps(row,0),shape(0,1),0.0,tmp7);

         DArray<8> tmp8;
         Contract(1.0,tmp7,shape(2,4),peps(row,0),shape(1,2),0.0,tmp8);

         LO.clear();
         Contract(1.0,tmp8,shape(5,3,6),env.gb(row-1)[0],shape(0,1,2),0.0,LO);

      }
      else{//col != 0

         DArray<8> tmp8;
         Gemm(CblasTrans,CblasNoTrans,1.0,LO,env.gt(row)[col],0.0,tmp8);

         DArray<9> tmp9;
         Contract(1.0,tmp8,shape(0,5),peps(row+1,col),shape(0,1),0.0,tmp9);

         tmp8.clear();
         Contract(1.0,tmp9,shape(0,4,6),peps(row+1,col),shape(0,1,2),0.0,tmp8);

         tmp9.clear();
         Contract(1.0,tmp8,shape(0,4),peps(row,col),shape(0,1),0.0,tmp9);

         tmp8.clear();
         Contract(1.0,tmp9,shape(0,4,6),peps(row,col),shape(0,1,2),0.0,tmp8);

         LO.clear();
         Contract(1.0,tmp8,shape(0,4,6),env.gb(row-1)[col],shape(0,1,2),0.0,LO);

      }

   }

}
