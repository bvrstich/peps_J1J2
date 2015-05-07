#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <complex>
#include <omp.h>

using std::cout;
using std::endl;
using std::vector;
using std::complex;
using std::ofstream;

#include "include.h"

using namespace global;

/** 
 * empty constructor
 */
Environment::Environment(){ }

/** 
 * constructor with allocation
 * @param D_in bond dimension of peps state
 * @param D_aux_in contraction bond dimension
 * @param comp_sweeps_in sets the number of sweeps done for MPO compression
 */
Environment::Environment(int D_in,int D_aux_in,int comp_sweeps_in){

   t.resize(Ly - 2);
   b.resize(Ly - 2);

   D = D_in;
   D_aux = D_aux_in;
   comp_sweeps = comp_sweeps_in;

   //allocate the memory
   
   //bottom
   int tmp = D*D;

   for(int i = 0;i < Ly - 2;++i){

      if(tmp < D_aux){

         b[i] = MPO<double>(Lx,D,tmp);
         tmp *= D*D;

      }
      else
         b[i] = MPO<double>(Lx,D,D_aux);

   }
   
   //top
   tmp = D*D;

   for(int i = Ly - 3;i >= 0;--i){

      if(tmp < D_aux){

         t[i] = MPO<double>(Lx,D,tmp);
         tmp *= D*D;

      }
      else
         t[i] = MPO<double>(Lx,D,D_aux);

   }

}

/** 
 * copy constructor with allocation
 */
Environment::Environment(const Environment &env_copy){

   t = env_copy.gt();
   b = env_copy.gb();

   D = env_copy.gD();
   D_aux = env_copy.gD_aux();

   comp_sweeps = env_copy.gcomp_sweeps();

}

/**
 * empty destructor
 */
Environment::~Environment(){ }

/**
 * construct the enviroment mps's for the input PEPS
 * @param option if 'L' construct full left environment
 *               if 'R' construct full right environment
 *               if 'T' construct full top environment
 *               if 'B' construct full bottom environment
 *               if 'A' construct all environments
 * @param peps input PEPS<double>
 * @param D_aux dimension to which environment will be compressed
 */
void Environment::calc(const char option,PEPS<double> &peps){

   if(option == 'A'){

#pragma omp parallel sections
      {

#pragma omp section
         {

            b[0].fill('b',peps);
            b[0].canonicalize(Right,false);

            for(int i = 1;i < Ly - 2;++i)
               this->add_layer('b',i,peps);
  
         }
#pragma omp section
         {

            t[Ly - 3].fill('t',peps);

            for(int i = Ly - 4;i >= 0;--i)
               this->add_layer('t',i,peps);

         }

      } 

   }
   else if(option == 'B'){

      b[0].fill('b',peps);

      for(int i = 1;i < Ly - 2;++i)
         this->add_layer('b',i,peps);

   }
   else if(option == 'T'){

      t[Ly - 3].fill('t',peps);

      for(int i = Ly - 4;i >= 0;--i)
         this->add_layer('t',i,peps);

   }

}

/**
 * test if the enviroment is correctly contracted
 */
void Environment::test(){

   for(int i = 0;i < Ly - 3;++i)
      cout << i + 2 << "\t" << b[i + 1].dot(t[i]) << endl;

}

/**
 * const version
 * @param row the row index
 * @return the top boundary 'MPO' environment on row 'row'
 */
const MPO<double> &Environment::gt(int row) const {

   return t[row];

}

/**
 * access version
 * @param row the row index
 * @return the top boundary 'MPO' environment on row 'row'
 */
MPO<double> &Environment::gt(int row) {

   return t[row];

}

/**
 * const version
 * @param row the row index
 * @return the bottom boundary 'MPO' environment on row 'row'
 */
const MPO<double> &Environment::gb(int row) const {

   return b[row];

}

/**
 * access version
 * @param row the row index
 * @return the bottom boundary 'MPO' environment on row 'row'
 */
MPO<double> &Environment::gb(int row) {

   return b[row];

}

/**
 * @return the auxiliary bond dimension for the contraction
 **/
int Environment::gD_aux() const {

   return D_aux;

}

/**
 * @return the auxiliary bond dimension for the contraction
 **/
int Environment::gD() const {

   return D;

}

/**
 * @return the number of sweeps performed during compression
 **/
int Environment::gcomp_sweeps() const {

   return comp_sweeps;

}


/**
 * set a new bond dimension
 */
void Environment::sD(int D_in) {

   D = D_in;

}

/**
 * set a new auxiliary bond dimension
 */
void Environment::sD_aux(int D_aux_in) {

   D_aux = D_aux_in;

}

/**
 * @return the full bottom boundary 'MPO'
 */
const vector< MPO<double> > &Environment::gb() const {

   return b;

}

/**
 * @return the full top boundary 'MPO'
 */
const vector< MPO<double> > &Environment::gt() const {

   return t;

}

/**
 * construct the (t or b) environment on row/col 'rc' by adding a the appropriate peps row/col and compressing the boundary MPO
 * @param option 't'op or 'b'ottom
 * @param row row index
 * @param peps the input PEPS<double> object 
 */
void Environment::add_layer(const char option,int row,PEPS<double> &peps){

   if(option == 'b'){

      //initialize using svd
      init_svd(option,row,peps);
      
      std::vector< DArray<4> > R(Lx+1);

      //first construct rightmost operator
      R[Lx].resize(1,1,1,1);
      R[Lx] = 1.0;

      //now move from right to left to construct the rest
      for(int col = Lx - 1;col > 0;--col){

         DArray<6> tmp6;
         Contract(1.0,b[row - 1][col],shape(3),R[col+1],shape(0),0.0,tmp6);

         DArray<7> tmp7;
         Contract(1.0,tmp6,shape(1,3),peps(row,col),shape(3,4),0.0,tmp7);

         tmp6.clear();
         Contract(1.0,tmp7,shape(1,2,6),peps(row,col),shape(3,4,2),0.0,tmp6);

         Contract(1.0,tmp6,shape(3,5,1),b[row][col],shape(1,2,3),0.0,R[col]);

      }

      R[0].resize(1,1,1,1);
      R[0] = 1.0;

      cout << cost_function('b',row,0,peps,R) << endl;
/*
      int iter = 0;

//      while(iter < comp_sweeps){

         //now start sweeping to get the compressed boundary MPO
         DArray<6> tmp6;
         Contract(1.0,b[row - 1][0],shape(3),R[0],shape(0),0.0,tmp6);

         tmp7.clear();
         Contract(1.0,peps(row,0),shape(3,4),tmp6,shape(2,4),0.0,tmp7);

         tmp6.clear();
         Contract(1.0,peps(row,0),shape(2,3,4),tmp7,shape(2,4,5),0.0,tmp6);

         b[row][0] = tmp6.reshape_clear(shape(1,D,D,tmp6.shape(5)));

         //QR
         DArray<2> tmp2;
         Geqrf(b[row][0],tmp2);

         //construct new left operator
         tmp7.clear();
         Contract(1.0,b[row-1][0],shape(1),peps(row,0),shape(3),0.0,tmp7);

         tmp8.clear();
         Contract(1.0,tmp7,shape(1,5),peps(row,0),shape(3,2),0.0,tmp8);

         tmp8bis.clear();
         Contract(1.0,tmp8,shape(3,6),b[row][0],shape(1,2),0.0,tmp8bis);

         R[0] = tmp8bis.reshape_clear(shape(tmp8bis.shape(1),tmp8bis.shape(3),tmp8bis.shape(5),tmp8bis.shape(7)));

         //now for the rest of the rightgoing sweep.
         for(int i = 1;i < Lx-1;++i){

            tmp6.clear();
            Contract(1.0,R[i - 1],shape(0),b[row - 1][i],shape(0),0.0,tmp6);

            tmp7.clear();
            Contract(1.0,tmp6,shape(0,3),peps(row,i),shape(0,3),0.0,tmp7);

            tmp6.clear();
            Contract(1.0,tmp7,shape(0,2,5),peps(row,i),shape(0,3,2),0.0,tmp6);

            DArray<6> tmp6bis;
            Permute(tmp6,shape(0,2,4,1,3,5),tmp6bis);

            Gemm(CblasNoTrans,CblasNoTrans,1.0,tmp6bis,R[i],0.0,b[row][i]);

            //QR
            tmp2.clear();
            Geqrf(b[row][i],tmp2);

            //construct new left operator
            Gemm(CblasTrans,CblasNoTrans,1.0,tmp6bis,b[row][i],0.0,R[i]);

         }

         //rightmost site
         tmp6.clear();
         Contract(1.0,R[Lx - 2],shape(0),b[row - 1][Lx - 1],shape(0),0.0,tmp6);

         tmp7.clear();
         Contract(1.0,tmp6,shape(0,3),peps(row,Lx - 1),shape(0,3),0.0,tmp7);

         tmp6.clear();
         Contract(1.0,tmp7,shape(0,2,5),peps(row,Lx - 1),shape(0,3,2),0.0,tmp6);

         b[row][Lx - 1] = tmp6.reshape_clear(shape(tmp6.shape(0),D,D,1));

         //LQ
         tmp2.clear();
         Gelqf(tmp2,b[row][Lx - 1]);

         //construct new right operator
         tmp7.clear();
         Contract(1.0,b[row - 1][Lx - 1],shape(1),peps(row,Lx - 1),shape(3),0.0,tmp7);

         tmp8.clear();
         Contract(1.0,tmp7,shape(1,5),peps(row,Lx - 1),shape(3,2),0.0,tmp8);

         tmp8bis.clear();
         Contract(1.0,tmp8,shape(3,6),b[row][Lx - 1],shape(1,2),0.0,tmp8bis);

         R[Lx - 2] = tmp8bis.reshape_clear(shape(tmp8bis.shape(0),tmp8bis.shape(2),tmp8bis.shape(4),tmp8bis.shape(6)));

         //back to the beginning with a leftgoing sweep
         for(int i = Lx-2;i > 0;--i){

            tmp6.clear();
            Contract(1.0,b[row - 1][i],shape(3),R[i],shape(0),0.0,tmp6);

            tmp7.clear();
            Contract(1.0,tmp6,shape(1,3),peps(row,i),shape(3,4),0.0,tmp7);

            tmp6.clear();
            Contract(1.0,tmp7,shape(1,2,6),peps(row,i),shape(3,4,2),0.0,tmp6);

            DArray<6> tmp6bis;
            Permute(tmp6,shape(0,2,4,3,5,1),tmp6bis);

            Gemm(CblasTrans,CblasNoTrans,1.0,R[i - 1],tmp6bis,0.0,b[row][i]);

            //LQ
            tmp2.clear();
            Gelqf(tmp2,b[row][i]);

            //construct right operator
            Gemm(CblasNoTrans,CblasTrans,1.0,tmp6bis,b[row][i],0.0,R[i - 1]);

         }

         //multiply the last L matrix with the first matrix:
         DArray<4> tmp4;
         Gemm(CblasNoTrans,CblasNoTrans,1.0,b[row][0],tmp2,0.0,tmp4);

         b[row][0] = std::move(tmp4);

         ++iter;

 //     }

      //redistribute the norm over the chain
      double nrm =  Nrm2(b[row][0]);

      //rescale the first site
      Scal((1.0/nrm), b[row][0]);

      //then multiply the norm over the whole chain
      b[row].scal(nrm);
*/
   }
   else {

      vector< DArray<4> > R(Lx - 1);

      //first construct rightmost operator
      DArray<7> tmp7;
      Contract(1.0,t[row + 1][Lx - 1],shape(1),peps(row+2,Lx - 1),shape(1),0.0,tmp7);

      DArray<8> tmp8;
      Contract(1.0,tmp7,shape(1,4),peps(row+2,Lx - 1),shape(1,2),0.0,tmp8);

      DArray<8> tmp8bis;
      Contract(1.0,tmp8,shape(3,6),t[row][Lx - 1],shape(1,2),0.0,tmp8bis);

      R[Lx - 2] = tmp8bis.reshape_clear(shape(tmp8bis.shape(0),tmp8bis.shape(2),tmp8bis.shape(4),tmp8bis.shape(6)));

      //now move from right to left to construct the rest
      for(int i = Lx - 2;i > 0;--i){

         DArray<6> tmp6;
         Contract(1.0,t[row + 1][i],shape(3),R[i],shape(0),0.0,tmp6);

         tmp7.clear();
         Contract(1.0,tmp6,shape(1,3),peps(row+2,i),shape(1,4),0.0,tmp7);

         tmp6.clear();
         Contract(1.0,tmp7,shape(1,2,5),peps(row+2,i),shape(1,4,2),0.0,tmp6);

         Contract(1.0,tmp6,shape(3,5,1),t[row][i],shape(1,2,3),0.0,R[i - 1]);

      }

      int iter = 0;

      while(iter < comp_sweeps){

         //now start sweeping to get the compressed boundary MPO
         DArray<6> tmp6;
         Contract(1.0,t[row + 1][0],shape(3),R[0],shape(0),0.0,tmp6);

         tmp7.clear();
         Contract(1.0,peps(row+2,0),shape(1,4),tmp6,shape(2,4),0.0,tmp7);

         tmp6.clear();
         Contract(1.0,peps(row+2,0),shape(1,2,4),tmp7,shape(4,1,5),0.0,tmp6);

         t[row][0] = tmp6.reshape_clear(shape(1,D,D,tmp6.shape(5)));

         //QR
         DArray<2> tmp2;
         Geqrf(t[row][0],tmp2);

         //construct new left operator
         tmp7.clear();
         Contract(1.0,t[row+1][0],shape(1),peps(row+2,0),shape(1),0.0,tmp7);

         tmp8.clear();
         Contract(1.0,tmp7,shape(1,4),peps(row+2,0),shape(1,2),0.0,tmp8);

         tmp8bis.clear();
         Contract(1.0,tmp8,shape(3,6),t[row][0],shape(1,2),0.0,tmp8bis);

         R[0] = tmp8bis.reshape_clear(shape(tmp8bis.shape(1),tmp8bis.shape(3),tmp8bis.shape(5),tmp8bis.shape(7)));

         //now for the rest of the rightgoing sweep.
         for(int i = 1;i < Lx-1;++i){

            tmp6.clear();
            Contract(1.0,R[i - 1],shape(0),t[row + 1][i],shape(0),0.0,tmp6);

            tmp7.clear();
            Contract(1.0,tmp6,shape(0,3),peps(row+2,i),shape(0,1),0.0,tmp7);

            tmp6.clear();
            Contract(1.0,tmp7,shape(0,2,4),peps(row+2,i),shape(0,1,2),0.0,tmp6);

            DArray<6> tmp6bis;
            Permute(tmp6,shape(0,2,4,1,3,5),tmp6bis);

            Gemm(CblasNoTrans,CblasNoTrans,1.0,tmp6bis,R[i],0.0,t[row][i]);

            //QR
            tmp2.clear();
            Geqrf(t[row][i],tmp2);

            //construct left operator
            Gemm(CblasTrans,CblasNoTrans,1.0,tmp6bis,t[row][i],0.0,R[i]);

         }

         //rightmost site
         tmp6.clear();
         Contract(1.0,R[Lx - 2],shape(0),t[row + 1][Lx - 1],shape(0),0.0,tmp6);

         tmp7.clear();
         Contract(1.0,tmp6,shape(0,3),peps(row+2,Lx - 1),shape(0,1),0.0,tmp7);

         tmp6.clear();
         Contract(1.0,tmp7,shape(0,2,4),peps(row+2,Lx - 1),shape(0,1,2),0.0,tmp6);

         t[row][Lx - 1] = tmp6.reshape_clear(shape(tmp6.shape(0),D,D,1));

         //LQ
         tmp2.clear();
         Gelqf(tmp2,t[row][Lx - 1]);

         //construct new right operator
         tmp7.clear();
         Contract(1.0,t[row + 1][Lx - 1],shape(1),peps(row+2,Lx - 1),shape(1),0.0,tmp7);

         tmp8.clear();
         Contract(1.0,tmp7,shape(1,4),peps(row+2,Lx - 1),shape(1,2),0.0,tmp8);

         tmp8bis.clear();
         Contract(1.0,tmp8,shape(3,6),t[row][Lx - 1],shape(1,2),0.0,tmp8bis);

         R[Lx - 2] = tmp8bis.reshape_clear(shape(tmp8bis.shape(0),tmp8bis.shape(2),tmp8bis.shape(4),tmp8bis.shape(6)));

         //back to the beginning with a leftgoing sweep
         for(int i = Lx-2;i > 0;--i){

            tmp6.clear();
            Contract(1.0,t[row + 1][i],shape(3),R[i],shape(0),0.0,tmp6);

            tmp7.clear();
            Contract(1.0,tmp6,shape(1,3),peps(row+2,i),shape(1,4),0.0,tmp7);

            tmp6.clear();
            Contract(1.0,tmp7,shape(1,5,2),peps(row+2,i),shape(1,2,4),0.0,tmp6);

            DArray<6> tmp6bis;
            Permute(tmp6,shape(0,2,4,3,5,1),tmp6bis);

            Gemm(CblasTrans,CblasNoTrans,1.0,R[i - 1],tmp6bis,0.0,t[row][i]);

            //LQ
            tmp2.clear();
            Gelqf(tmp2,t[row][i]);

            //construct right operator
            Gemm(CblasNoTrans,CblasTrans,1.0,tmp6bis,t[row][i],0.0,R[i - 1]);

         }

         //multiply the last L matrix with the first matrix:
         DArray<4> tmp4;
         Gemm(CblasNoTrans,CblasNoTrans,1.0,t[row][0],tmp2,0.0,tmp4);

         t[row][0] = std::move(tmp4);

         ++iter;

      }

      //redistribute the norm over the chain
      double nrm =  Nrm2(t[row][0]);

      //rescale the first site
      Scal((1.0/nrm), t[row][0]);

      //then multiply the norm over the whole chain
      t[row].scal(nrm);

   }

}

/**
 * construct the (t or b) environment on row/col 'rc' by adding a the appropriate peps row/col and compressing the boundary MPO
 * @param option 't'op or 'b'ottom
 * @param row row index
 * @param col column index
 * @param peps the input PEPS<double> object 
 */
double Environment::cost_function(const char option,int row,int col,const PEPS<double> &peps,const std::vector< DArray<4> > &R){

   if(option == 'b'){

      //environment of b is completely unitary
      double val = Dot(b[row][col],b[row][col]);

      //add row -1 to right hand side (col + 1)
      DArray<6> tmp6;
      Contract(1.0,b[row - 1][col],shape(3),R[col+1],shape(0),0.0,tmp6);

      DArray<7> tmp7;
      Contract(1.0,tmp6,shape(1,3),peps(row,col),shape(3,4),0.0,tmp7);

      tmp6.clear();
      Contract(1.0,tmp7,shape(6,1,2),peps(row,col),shape(2,3,4),0.0,tmp6);

      DArray<4> tmp4;
      Contract(1.0,tmp6,shape(3,5,1),b[row][col],shape(1,2,3),0.0,tmp4);

      val -= 2.0 * Dot(tmp4,R[col]);

      return val;

   }
   else{

      return 0.0;

   }

}

/**
 * initialize the environment on 'row' by performing an svd-compression on the 'full' environment b[row-1] * peps(row,...) * peps(row,...)
 * output is right canonical, which is needed for the compression algorithm!
 * @param option 'b'ottom or 't'op environment
 * @param row index of the row to be added into the environment
 */
void Environment::init_svd(char option,int row,const PEPS<double> &peps){

   if(option == 'b'){

      //add rightmost peps to right bottom environment
      DArray<5> tmp5;
      Contract(1.0,peps(row,Lx-1),shape(3,4),b[row-1][Lx-1],shape(2,3),0.0,tmp5);

      DArray<6> tmp6;
      Contract(1.0,peps(row,Lx-1),shape(2,3),tmp5,shape(2,4),0.0,tmp6);

      DArray<6> tmp6bis;
      Permute(tmp6,shape(0,3,5,1,4,2),tmp6bis);

      //now svd the large object
      DArray<1> S;
      DArray<4> U;

      Gesvd('S','S',tmp6bis,S,U,b[row][Lx-1],D_aux);

      //paste S to left 
      Dimm(U,S);

      //the rest of the columns
      //for(int col = Lx - 2;col > 0;--col){
     int col = Lx - 2; 
         //add next bottom to 'U' from previous column
         tmp6.clear();
         Contract(1.0,b[row-1][col],shape(3),U,shape(2),0.0,tmp6);

         //add next peps(row,col) to intermediate
         DArray<7> tmp7;
         Contract(1.0,peps(row,col),shape(3,4),tmp6,shape(2,4),0.0,tmp7);

         //and again
         tmp6.clear();
         Contract(1.0,peps(row,col),shape(2,3,4),tmp7,shape(2,4,5),0.0,tmp6);

         //permute!
         tmp6bis.clear();
         Permute(tmp6,shape(0,2,4,1,3,5),tmp6bis);

         //and svd!
         S.clear();
         U.clear();

         Gesvd('S','S',tmp6bis,S,U,b[row][col],D_aux);

         //paste S to left for next iteration
         Dimm(U,S);

      //}
     col = 1; 
      //add next bottom to 'U' from previous column
         tmp6.clear();
         Contract(1.0,b[row-1][col],shape(3),U,shape(2),0.0,tmp6);

         //add next peps(row,col) to intermediate
         tmp7.clear();
         Contract(1.0,peps(row,col),shape(3,4),tmp6,shape(2,4),0.0,tmp7);

         //and again
         tmp6.clear();
         Contract(1.0,peps(row,col),shape(2,3,4),tmp7,shape(2,4,5),0.0,tmp6);

         //permute!
         tmp6bis.clear();
         Permute(tmp6,shape(0,2,4,1,3,5),tmp6bis);

         //and svd!
         S.clear();
         U.clear();

         Gesvd('S','S',tmp6bis,S,U,b[row][col],D_aux);
         cout << b[row][col].shape() << endl;

         //paste S to left for next iteration
         Dimm(U,S);
/*
      //finally first (leftmost) site
      tmp6.clear();
      Contract(1.0,b[row-1][0],shape(3),U,shape(2),0.0,tmp6);

      //add next peps(row,col) to intermediate
      DArray<7> tmp7;
      Contract(1.0,peps(row,0),shape(3,4),tmp6,shape(2,4),0.0,tmp7);

      //and again
      tmp6.clear();
      Contract(1.0,peps(row,0),shape(2,3,4),tmp7,shape(2,4,5),0.0,tmp6);
  */
   }

}
