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

         //attach top-row peps to left
         DArray<8> tmp8;
         Contract(1.0,L,shape(0),peps(Ly-1,col),shape(0),0.0,tmp8);

         //and again
         DArray<7> tmp7;
         Contract(1.0,tmp8,shape(0,4,5),peps(Ly-1,col),shape(0,1,2),0.0,tmp7);

         //add second-row PEPS
         tmp8.clear();
         Contract(1.0,tmp7,shape(0,3),peps(Ly-2,col),shape(0,1),0.0,tmp8);

         //and again
         tmp7.clear();
         Contract(1.0,tmp8,shape(0,3,5),peps(Ly-2,col),shape(0,1,2),0.0,tmp7);

         //finally
         L.clear();
         Contract(1.0,tmp7,shape(0,3,5),env.gb(Ly-3)[col],shape(0,1,2),0.0,L);

      }

   }

   /** 
    * init the right renormalized operator for the middle rows: 
    * @param row index of the lowest row
    * @param peps The PEPS object
    * @param R vector containing the right operators on exit
    */
   void init_ro(int row,const PEPS<double> &peps,vector< DArray<6> > &RO){

      DArray<8> tmp8;
      DArray<8> tmp8bis;

      DArray<9> tmp9;
      DArray<9> tmp9bis;

      RO[Lx-1].resize( shape(1,1,1,1,1,1) );
      RO[Lx-1] = 1.0;

      //now move from right to left, constructing the rest
      for(int col = Lx - 1;col > 0;--col){

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
         R[Lx-1] = 1.0;

         DArray<7> tmp7;
         DArray<8> tmp8;

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

         R[Lx-1].resize( shape(1,1,1,1,1) );
         R[Lx-1] = 1.0;

         DArray<7> tmp7;
         DArray<8> tmp8;

         for(int i = Lx - 1;i > 0;--i){

            tmp7.clear();
            Contract(1.0,env.gb(Ly-3)[i],shape(3),R[i],shape(4),0.0,tmp7);

            tmp8.clear();
            Contract(1.0,peps(Ly-2,i),shape(3,4),tmp7,shape(2,6),0.0,tmp8);

            tmp7.clear();
            Contract(1.0,peps(Ly-2,i),shape(2,3,4),tmp8,shape(2,4,7),0.0,tmp7);

            tmp8.clear();
            Contract(1.0,peps(Ly-1,i),shape(3,4),tmp7,shape(3,6),0.0,tmp8);

            R[i-1].clear();
            Contract(1.0,peps(Ly-1,i),shape(1,2,3,4),tmp8,shape(1,2,4,7),0.0,R[i-1]);

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

   /** 
    * init the right renormalized operator for the middle rows: 
    * @param row index of the lowest row
    * @param peps The PEPS object
    * @param R vector containing the right operators on exit
    */
   double rescale_norm(int row,PEPS<double> &peps,vector< DArray<6> > &RO){

      DArray<8> tmp8;
      DArray<8> tmp8bis;

      DArray<9> tmp9;
      DArray<9> tmp9bis;

      DArray<6> tmp6;

      //first add bottom to right unity
      tmp8.clear();
      Gemm(CblasNoTrans,CblasTrans,1.0,env.gb(row-1)[0],RO[0],0.0,tmp8);

      tmp8bis.clear();
      Permute(tmp8,shape(2,7,0,1,3,4,5,6),tmp8bis);

      //add regular peps on lower site
      tmp9.clear();
      Gemm(CblasNoTrans,CblasNoTrans,1.0,peps(row,0),tmp8bis,0.0,tmp9);

      tmp9bis.clear();
      Permute(tmp9,shape(2,4,8,0,1,3,5,6,7),tmp9bis);

      //and another regular peps on lower site
      tmp8.clear();
      Gemm(CblasNoTrans,CblasNoTrans,1.0,peps(row,0),tmp9bis,0.0,tmp8);

      tmp8bis.clear();
      Permute(tmp8,shape(3,7,1,6,0,2,4,5),tmp8bis);

      //add regular peps on upper site
      tmp9.clear();
      Gemm(CblasNoTrans,CblasNoTrans,1.0,peps(row+1,0),tmp8bis,0.0,tmp9);

      tmp9bis.clear();
      Permute(tmp9,shape(2,3,4,0,1,5,6,7,8),tmp9bis);

      //yet another regular on upper site
      tmp8.clear();
      Gemm(CblasNoTrans,CblasNoTrans,1.0,peps(row+1,0),tmp9bis,0.0,tmp8);

      tmp8bis.clear();
      Permute(tmp8,shape(1,3,7,0,2,4,5,6),tmp8bis);

      //finally top environment for closure
      Gemm(CblasNoTrans,CblasNoTrans,1.0,env.gt(row)[0],tmp8bis,0.0,tmp6);

      //---------------------
      //--- rescale stuff ---
      //---------------------
      double full_nrm = tmp6(0,0,0,0,0,0);

      //first bottom environment
      env.gb(row-1).scal(1.0/full_nrm);

      //then the R operators
      double nrm = pow(full_nrm,1.0/(double)Lx);
      double scl = 1.0;

      for(int col = Lx - 2;col >= 0;--col){

         scl *= nrm;
         Scal(1.0/scl,RO[col]);

      }

      //lastly rscale the tensors themselves
      nrm = sqrt(nrm);

      for(int col = 0;col < Lx;++col)
         Scal(1.0/nrm,peps(row-1,col));

      return full_nrm;

   }
   /** 
    * init the right renormalized operator for the middle rows: 
    * @param row index of the lowest row
    * @param peps The PEPS object
    * @param R vector containing the right operators on exit
    */
   double rescale_norm(char option,PEPS<double> &peps,vector< DArray<5> > &R){

      if(option == 'b'){

         DArray<7> tmp7;
         Contract(1.0,env.gt(0)[0],shape(3),R[0],shape(0),0.0,tmp7);

         DArray<8> tmp8;
         Contract(1.0,tmp7,shape(1,3),peps(1,0),shape(1,4),0.0,tmp8);

         tmp7.clear();
         Contract(1.0,tmp8,shape(1,6,2),peps(1,0),shape(1,2,4),0.0,tmp7);

         tmp8.clear();
         Contract(1.0,tmp7,shape(4,1),peps(0,0),shape(1,4),0.0,tmp8);

         DArray<5> tmp5;
         Contract(1.0,tmp8,shape(4,6,7,1),peps(0,0),shape(1,2,3,4),0.0,tmp5);

         //---------------------
         //--- rescale stuff ---
         //---------------------
         double full_nrm = tmp5(0,0,0,0,0);

         double nrm = pow(full_nrm,1.0/(double)Ly);

         double scl = 1.0;

         for(int row = Ly - 3;row >= 0;--row){

            scl *= nrm;
            env.gt(row).scal(1.0/scl);

         }

         //rescale the R tensors
         scl = pow(scl,1.0/(double)Lx);

         nrm = pow(nrm,1.0/(double)Lx);

         nrm = scl*nrm*nrm;

         scl = 1.0;

         for(int col = Ly - 2;col >= 0;--col){

            scl *= nrm;
            Scal(1.0/scl,R[col]);

         }

         //finally rescale all the tensors
         peps.scal(1.0/sqrt(full_nrm));

         return full_nrm;

      }
      else{ //top 2 rows

         DArray<7> tmp7;
         DArray<8> tmp8;

         tmp7.clear();
         Contract(1.0,env.gb(Ly-3)[0],shape(3),R[0],shape(4),0.0,tmp7);

         tmp8.clear();
         Contract(1.0,peps(Ly-2,0),shape(3,4),tmp7,shape(2,6),0.0,tmp8);

         tmp7.clear();
         Contract(1.0,peps(Ly-2,0),shape(2,3,4),tmp8,shape(2,4,7),0.0,tmp7);

         tmp8.clear();
         Contract(1.0,peps(Ly-1,0),shape(3,4),tmp7,shape(3,6),0.0,tmp8);

         DArray<5> tmp5;
         Contract(1.0,peps(Ly-1,0),shape(1,2,3,4),tmp8,shape(1,2,4,7),0.0,tmp5);

         //---------------------
         //--- rescale stuff ---
         //---------------------
         double full_nrm = tmp5(0,0,0,0,0);

         //first bottom environment
         env.gb(Ly-3).scal(1.0/full_nrm);

         //then the R operators
         double nrm = pow(full_nrm,1.0/(double)Lx);
         double scl = 1.0;

         for(int col = Lx - 2;col >= 0;--col){

            scl *= nrm;
            Scal(1.0/scl,R[col]);

         }

         //lastly rescale the tensors themselves
         nrm = sqrt(nrm);

         for(int col = 0;col < Lx;++col)
            Scal(1.0/nrm,peps(Ly-3,col));

         return full_nrm;

      }

   }

}
