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
    * update left renormalized operator on site col 
    * @param option == 't'op ,'b'ottom, 'l'eft or 'r'ight
    * @param rc is row or column index, col for t,b row for r,l
    */
   void update_L(char option,int rc,const PEPS<double> &peps,DArray<3> &L){
/*
      if(option == 'b'){//bottom

         if(rc == 0){

            DArray<7> tmp7;
            Contract(1.0,env.gt(0)[0],shape(1),peps(0,0),shape(1),0.0,tmp7);

            DArray<8> tmp8;
            Contract(1.0,tmp7,shape(1,4),peps(0,0),shape(1,2),0.0,tmp8);

            L = tmp8.reshape_clear( shape(tmp8.shape(1),tmp8.shape(4),tmp8.shape(7)) );

         }
         else if(rc < Lx - 2){

            //construct left renormalized operator for next site: first construct intermediary
            DArray<5> tmp5;
            Contract(1.0,L,shape(0),env.gt(0)[rc],shape(0),0.0,tmp5);

            DArray<6> tmp6;
            Contract(1.0,tmp5,shape(0,2),peps(0,rc),shape(0,1),0.0,tmp6);

            tmp5.clear();
            Contract(1.0,tmp6,shape(0,1,3),peps(0,rc),shape(0,1,2),0.0,tmp5);

            L = tmp5.reshape_clear(shape(env.gt(0)[rc].shape(3),peps(0,rc).shape(4),peps(0,rc).shape(4)));

         }
         else{//nothing, no update necessary

         }

      }
      else if(option == 't'){//top

         if(rc == 0){

            DArray<7> tmp7;
            Contract(1.0,peps(Ly-1,0),shape(3),env.gb(Ly-2)[0],shape(2),0.0,tmp7);

            DArray<8> tmp8;
            Contract(1.0,peps(Ly-1,0),shape(2,3),tmp7,shape(2,5),0.0,tmp8);

            L = tmp8.reshape_clear( shape(tmp8.shape(2),tmp8.shape(5),tmp8.shape(7)) );

         }
         else{

            //construct left renormalized operators for next site: first construct intermediary
            DArray<5> tmp5;
            Contract(1.0,env.gb(Ly-2)[rc],shape(0),L,shape(2),0.0,tmp5);

            DArray<6> tmp6;
            Contract(1.0,peps(Ly-1,rc),shape(3,0),tmp5,shape(1,4),0.0,tmp6);

            //finally construct new unity on the left
            tmp5.clear();
            Contract(1.0,peps(Ly-1,rc),shape(2,3,0),tmp6,shape(1,3,5),0.0,tmp5);

            L = tmp5.reshape_clear(shape(peps(Ly-1,rc).shape(4),peps(Ly-1,rc).shape(4),env.gb(Ly-2)[rc].shape(3)));

         }

      }
      else if(option == 'l'){//left

         if(rc == 0){

            DArray<7> tmp7;
            Contract(1.0,peps(0,0),shape(4),env.gr(0)[0],shape(2),0.0,tmp7);

            DArray<8> tmp8;
            Contract(1.0,peps(0,0),shape(2,4),tmp7,shape(2,5),0.0,tmp8);

            L = tmp8.reshape_clear( shape(tmp8.shape(1),tmp8.shape(4),tmp8.shape(7)) );

         }
         else{

            //construct left renormalized operators for next site: first construct intermediary
            DArray<5> tmp5;
            Contract(1.0,env.gr(0)[rc],shape(0),L,shape(2),0.0,tmp5);

            DArray<6> tmp6;
            Contract(1.0,peps(rc,0),shape(3,4),tmp5,shape(4,1),0.0,tmp6);

            // construct new unity on the left
            tmp5.clear();
            Contract(1.0,peps(rc,0),shape(2,3,4),tmp6,shape(2,5,3),0.0,tmp5);

            L = tmp5.reshape_clear(shape(peps(rc,0).shape(1),peps(rc,0).shape(1),env.gr(0)[rc].shape(3)));

         }

      }
      else{//right

         if(rc == 0){

            DArray<7> tmp7;
            Contract(1.0,env.gl(Lx - 2)[0],shape(1),peps(0,Lx - 1),shape(0),0.0,tmp7);

            DArray<8> tmp8;
            Contract(1.0,tmp7,shape(1,4),peps(0,Lx - 1),shape(0,2),0.0,tmp8);

            L = tmp8.reshape_clear( shape(tmp8.shape(1),tmp8.shape(2),tmp8.shape(5)) );

         }
         else{

            //construct left renormalized operators for next site: first construct intermediary
            DArray<5> tmp5;
            Contract(1.0,L,shape(0),env.gl(Lx - 2)[rc],shape(0),0.0,tmp5);

            DArray<6> tmp6;
            Contract(1.0,tmp5,shape(0,2),peps(rc,Lx - 1),shape(3,0),0.0,tmp6);

            tmp5.clear();
            Contract(1.0,tmp6,shape(0,1,4),peps(rc,Lx-1),shape(3,0,2),0.0,tmp5);

            L = tmp5.reshape_clear(shape(env.gl(Lx - 2)[rc].shape(3),peps(rc,Lx - 1).shape(1),peps(rc,Lx - 1).shape(1)));

         }

      }
  */
   }

   /** 
    * init the right renormalized operator for the middle rows: 
    * @param option 'H'orizontal or 'V'ertical
    * @param rc 'row' index for Horizontal, 'col' index for Vertical
    * @param peps The PEPS object
    * @param R vector containing the right operators on exit
    */
   void init_ro(bool is_local,char option,int rc,const PEPS<double> &peps,vector< DArray<4> > &RO){
/*
      if(option == 'H'){

         DArray<6> tmp6;
         DArray<6> tmp6bis;
         DArray<7> tmp7;
         DArray<8> tmp8;
         DArray<8> tmp8bis;

         int stop;

         if(is_local)
            stop = 0;
         else
            stop = 1;

         //paste top peps 'operators'
         Contract(1.0,env.gt(rc)[Lx - 1],shape(1),peps(rc,Lx-1),shape(1),0.0,tmp7);

         Contract(1.0,tmp7,shape(1,4),peps(rc,Lx-1),shape(1,2),0.0,tmp8);

         Contract(1.0,tmp8,shape(3,6),env.gb(rc-1)[Lx-1],shape(1,2),0.0,tmp8bis);

         //move to a DArray<3> object
         RO[Lx - 2] = tmp8bis.reshape_clear(shape(env.gt(rc)[Lx - 1].shape(0),peps(rc,Lx-1).shape(0),peps(rc,Lx-1).shape(0),env.gb(rc-1)[Lx - 1].shape(0)));

         //now construct the middle operators
         for(int col = Lx - 2;col > stop;--col){

            tmp6.clear();
            Contract(1.0,env.gt(rc)[col],shape(3),RO[col],shape(0),0.0,tmp6);

            tmp7.clear();
            Contract(1.0,tmp6,shape(1,3),peps(rc,col),shape(1,4),0.0,tmp7);

            tmp6.clear();
            Contract(1.0,tmp7,shape(1,2,5),peps(rc,col),shape(1,4,2),0.0,tmp6);

            tmp6bis.clear();
            Permute(tmp6,shape(0,2,4,3,5,1),tmp6bis);

            RO[col - 1].clear();
            Gemm(CblasNoTrans,CblasTrans,1.0,tmp6bis,env.gb(rc - 1)[col],0.0,RO[col - 1]);

         }

      }
      else{//vertical, columns

         DArray<6> tmp6;
         DArray<6> tmp6bis;
         DArray<7> tmp7;
         DArray<8> tmp8;
         DArray<8> tmp8bis;

         int stop;

         if(is_local)
            stop = 0;
         else
            stop = 1;

         Contract(1.0,env.gl(rc - 1)[Ly - 1],shape(1),peps(Ly-1,rc),shape(0),0.0,tmp7);

         Contract(1.0,tmp7,shape(1,4),peps(Ly-1,rc),shape(0,2),0.0,tmp8);

         Contract(1.0,tmp8,shape(4,7),env.gr(rc)[Ly-1],shape(1,2),0.0,tmp8bis);

         //move to a DArray<3> object
         RO[Ly - 2] = tmp8bis.reshape_clear(shape(env.gl(rc - 1)[Ly - 1].shape(0),peps(Ly-1,rc).shape(3),peps(Ly-1,rc).shape(3),env.gr(rc)[Ly - 1].shape(0)));

         //now construct the middle operators
         for(int row = Ly - 2;row > stop;--row){

            tmp6.clear();
            Contract(1.0,env.gl(rc - 1)[row],shape(3),RO[row],shape(0),0.0,tmp6);

            tmp7.clear();
            Contract(1.0,tmp6,shape(1,3),peps(row,rc),shape(0,1),0.0,tmp7);

            tmp6.clear();
            Contract(1.0,tmp7,shape(1,2,4),peps(row,rc),shape(0,1,2),0.0,tmp6);

            tmp6bis.clear();
            Permute(tmp6,shape(0,2,4,3,5,1),tmp6bis);

            RO[row-1].clear();
            Gemm(CblasNoTrans,CblasTrans,1.0,tmp6bis,env.gr(rc)[row],0.0,RO[row - 1]);

         }

      }
  */
   }

   /** 
    * init the right renormalized operator for the top or bottom row, or left or right column
    * @param option == 'l'eft 'r'ight 'top' or 'b'ottom
    * @param R vector containing the right operators on exit
    */
   void init_ro(char option,const PEPS<double> &peps,vector< DArray<5> > &R){

      if(option == 'b'){

         DArray<7> tmp7;
         DArray<7> tmp7bis;

         //first the rightmost operator
         int M = peps(1,Lx-1).shape(0) * peps(1,Lx-1).shape(1) * peps(1,Lx-1).shape(2);
         int N = env.gb(0)[Lx-1].shape(0) * env.gb(0)[Lx-1].shape(1);
         int K = env.gb(0)[Lx-1].shape(2) * env.gb(0)[Lx-1].shape(3);

         DArray<6> tmp6(D,D,d,D,D,D);
         blas::gemm(CblasRowMajor, CblasNoTrans, CblasTrans, M, N, K, 1.0,peps(1,Lx-1).data(),K,env.gb(0)[Lx-1].data(),K,0.0,tmp6.data(),N);

         DArray<6> tmp6bis;
         Permute(tmp6,shape(2,5,0,1,3,4),tmp6bis);

         M = peps(1,Lx-1).shape(0) * peps(1,Lx-1).shape(1);
         N = tmp6bis.shape(2) * tmp6bis.shape(3) * tmp6bis.shape(4) * tmp6bis.shape(5);
         K = peps(1,Lx-1).shape(2) * peps(1,Lx-1).shape(3);

         tmp6.resize( shape(D,D,tmp6bis.shape(2),tmp6bis.shape(3),tmp6bis.shape(4),tmp6bis.shape(5)) );
         blas::gemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, M, N, K, 1.0,peps(1,Lx-1).data(),K,tmp6bis.data(),N,0.0,tmp6.data(),N);

         tmp6bis.clear();
         Permute(tmp6,shape(1,3,0,2,4,5),tmp6bis);

         M = env.gt(0)[Lx - 1].shape(0);
         N = tmp6bis.shape(2) * tmp6bis.shape(3) * tmp6bis.shape(4) * tmp6bis.shape(5);
         K = env.gt(0)[Lx - 1].shape(1) * env.gt(0)[Lx - 1].shape(2);

         R[Lx - 2].resize(shape(M,tmp6bis.shape(2),tmp6bis.shape(3),tmp6bis.shape(4),tmp6bis.shape(5)));
         blas::gemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, M, N, K, 1.0,env.gt(0)[Lx-1].data(),K,tmp6bis.data(),N,0.0,R[Lx - 2].data(),N);

         for(int i = Lx - 2;i > 0;--i){

            M = R[i].shape(0) * R[i].shape(1) * R[i].shape(2);
            N = env.gb(0)[i].shape(0) * env.gb(0)[i].shape(1) * env.gb(0)[i].shape(2);
            K = env.gb(0)[i].shape(3);

            tmp7.resize(shape( R[i].shape(0),D,D,D,D,D,D) );
            blas::gemm(CblasRowMajor, CblasNoTrans, CblasTrans, M, N, K, 1.0,R[i].data(),K,env.gb(0)[i].data(),K,0.0,tmp7.data(),N);

            tmp7bis.clear();
            Permute(tmp7,shape(0,1,3,4,5,6,2),tmp7bis);

            DArray<8> tmp8;
            Gemm(CblasNoTrans,CblasTrans,1.0,tmp7bis,peps(1,i),0.0,tmp8);

            DArray<8> tmp8bis;
            Permute(tmp8,shape(0,5,2,3,6,7,4,1),tmp8bis);

            tmp7.clear();
            Gemm(CblasNoTrans,CblasTrans,1.0,tmp8bis,peps(1,i),0.0,tmp7);

            tmp7bis.clear();
            Permute(tmp7,shape(6,4,0,5,1,2,3),tmp7bis);

            Gemm(CblasNoTrans,CblasNoTrans,1.0,env.gt(0)[i],tmp7bis,0.0,R[i - 1]);

         }

      }
      else if(option == 't'){

      }

   }


   /**
    * update left renormalized operator on site (row,col )
    * @param option 'H'orizonal or 'V'ertical
    * @param row index
    * @param col index
    * @param peps the input PEPS object
    * @param LO input old left renormalized operator, output new left renormalized operator
    */
   void update_L(char option,int row,int col,const PEPS<double> &peps,DArray<4> &LO){
      /*
         if(option == 'H'){

         if(col == 0){

      //paste top environment on
      DArray<7> tmp7;
      Contract(1.0,env.gt(row)[0],shape(1),peps(row,0),shape(1),0.0,tmp7);

      DArray<8> tmp8;
      Contract(1.0,tmp7,shape(1,4),peps(row,0),shape(1,2),0.0,tmp8);

      DArray<8> tmp8bis;
      Contract(1.0,tmp8,shape(3,6),env.gb(row-1)[0],shape(1,2),0.0,tmp8bis);

      //move to a DArray<4> object: order (top-env,(*this)-row,bottom-env)
      LO = tmp8bis.reshape_clear(shape(env.gt(row)[0].shape(3),peps(row,0).shape(4),peps(row,0).shape(4),env.gb(row-1)[0].shape(3)));

      }
      else if(col < Lx - 2){//middle

      //first attach top to left unity
      DArray<6> tmp6;
      Contract(1.0,env.gt(row)[col],shape(0),LO,shape(0),0.0,tmp6);

      //add peps to it, intermediary
      DArray<7> tmp7;
      Contract(1.0,tmp6,shape(3,0),peps(row,col),shape(0,1),0.0,tmp7);

      //finally construct new left unity
      tmp6.clear();
      Contract(1.0,tmp7,shape(2,0,4),peps(row,col),shape(0,1,2),0.0,tmp6);

      DArray<6> tmp6bis;
      Permute(tmp6,shape(0,3,5,1,2,4),tmp6bis);

      LO.clear();
      Gemm(CblasNoTrans,CblasNoTrans,1.0,tmp6bis,env.gb(row - 1)[col],0.0,LO);

      }
      else{

      //no update

      }

      }
      else{//Vertical

      if(row == 0){

      //paste left environment on
      DArray<7> tmp7;
      Contract(1.0,env.gl(col - 1)[0],shape(1),peps(0,col),shape(0),0.0,tmp7);

      DArray<8> tmp8;
      Contract(1.0,tmp7,shape(1,4),peps(0,col),shape(0,2),0.0,tmp8);

      DArray<8> tmp8bis;
      Contract(1.0,tmp8,shape(4,7),env.gr(col)[0],shape(1,2),0.0,tmp8bis);

      //move to a DArray<3> object: order (top-env,(*this)-row,bottom-env)
      LO = tmp8bis.reshape_clear(shape(env.gl(col - 1)[0].shape(3),peps(0,col).shape(1),peps(0,col).shape(1),env.gr(col)[0].shape(3)));

      }
      else if(row < Ly - 2){

      //first attach top to left unity, make intermediary
      DArray<6> tmp6;
      Contract(1.0,env.gl(col - 1)[row],shape(0),LO,shape(0),0.0,tmp6);

      //add peps to it
      DArray<7> tmp7;
      Contract(1.0,tmp6,shape(0,3),peps(row,col),shape(0,3),0.0,tmp7);

      tmp6.clear();
      Contract(1.0,tmp7,shape(0,5,2),peps(row,col),shape(0,2,3),0.0,tmp6);

      DArray<6> tmp6bis;
      Permute(tmp6,shape(0,2,4,1,3,5),tmp6bis);

      LO.clear();
      Gemm(CblasNoTrans,CblasNoTrans,1.0,tmp6bis,env.gr(col)[row],0.0,LO);

   }
      else{

         //nothing

      }

   }

   */
   }

}
