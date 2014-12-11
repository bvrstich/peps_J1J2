#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <complex>

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

   t.resize(Ly - 1);
   b.resize(Ly - 1);

   r.resize(Lx - 1);
   l.resize(Lx - 1);

   flag_b = false;
   flag_t = false;
   flag_r = false;
   flag_l = false;

   D = D_in;
   D_aux = D_aux_in;
   comp_sweeps = comp_sweeps_in;

   //allocate the memory
   
   //bottom
   int tmp = D*D;

   for(int i = 0;i < Ly - 1;++i){

      if(tmp < D_aux){

         b[i] = MPO<double>(Lx,D,tmp);
         tmp *= D*D;

      }
      else
         b[i] = MPO<double>(Lx,D,D_aux);

   }
   
   //top
   tmp = D*D;

   for(int i = Ly - 2;i >= 0;--i){

      if(tmp < D_aux){

         t[i] = MPO<double>(Lx,D,tmp);
         tmp *= D*D;

      }
      else
         t[i] = MPO<double>(Lx,D,D_aux);

   }

   //left
   tmp = D*D;

   for(int i = 0;i < Lx - 1;++i){

      if(tmp < D_aux){

         l[i] = MPO<double>(Ly,D,tmp);
         tmp *= D*D;

      }
      else
         l[i] = MPO<double>(Ly,D,D_aux);

   }
   
   //finally right
   tmp = D*D;

   for(int i = Lx - 2;i >= 0;--i){

      if(tmp < D_aux){

         r[i] = MPO<double>(Ly,D,tmp);
         tmp *= D*D;

      }
      else
         r[i] = MPO<double>(Ly,D,D_aux);

   }

}

/** 
 * copy constructor with allocation
 */
Environment::Environment(const Environment &env_copy){

   t = env_copy.gt();
   b = env_copy.gb();

   r = env_copy.gr();
   l = env_copy.gl();

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
void Environment::calc(const char option,const PEPS<double> &peps){

   if(option == 'B' || option == 'A'){

      b[0].fill('b',peps);

      for(int i = 1;i < Ly - 1;++i)
         this->add_layer('b',i,peps);

      flag_b = true;

   }

   if(option == 'T' || option == 'A'){

      t[Ly - 2].fill('t',peps);

      for(int i = Ly - 3;i >= 0;--i)
         this->add_layer('t',i,peps);

      flag_t = true;

   }

   if(option == 'R' || option == 'A'){

      r[Lx - 2].fill('r',peps);

      for(int i = Lx - 3;i >= 0;--i)
         this->add_layer('r',i,peps);

      flag_r = true;

   }

   if(option == 'L' || option == 'A'){

      l[0].fill('l',peps);

      for(int i = 1;i < Lx - 1;++i)
         this->add_layer('l',i,peps);

      flag_l = true;

   }

}

/**
 * test if the enviroment is correctly contracted
 */
void Environment::test(){

   cout << endl;
   cout << "FROM BOTTOM TO TOP" << endl;
   cout << endl;
   for(int i = 0;i < Ly - 1;++i)
      cout << i << "\t" << b[i].dot(t[i]) << endl;

   cout << endl;
   cout << "FROM LEFT TO RIGHT" << endl;
   cout << endl;
   for(int i = 0;i < Lx - 1;++i)
      cout << i << "\t" << l[i].dot(r[i]) << endl;
   cout << endl;
}

/**
 * const version
 * @param col the column index
 * @return the right boundary 'MPO' environment on column col
 */
const MPO<double> &Environment::gr(int col) const {

   return r[col];

}

/**
 * @param col the column index: access version
 * @return the right boundary 'MPO' environment on column col
 */
MPO<double> &Environment::gr(int col) {

   return r[col];

}

/**
 * const version
 * @param col the column index
 * @return the left boundary 'MPO' environment on column col
 */
const MPO<double> &Environment::gl(int col) const {

   return l[col];

}

/**
 * @param col the column index: access version
 * @return the left boundary 'MPO' environment on column col
 */
MPO<double> &Environment::gl(int col) {

   return l[col];

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
 * @return the full left boundary 'MPO'
 */
const vector< MPO<double> > &Environment::gl() const {

   return l;

}

/**
 * @return the full right boundary 'MPO'
 */
const vector< MPO<double> > &Environment::gr() const {

   return r;

}

/**
 * construct the (t,b,l or r) environment on row/col 'rc' by adding a the appropriate peps row/col and compressing the boundary MPO
 * @param option 't'op, 'b'ottom, 'l'eft or 'r'ight environment
 * @param rc row or column index
 * @param peps the input PEPS<double> object 
 */
void Environment::add_layer(const char option,int rc,const PEPS<double> &peps){

   if(option == 'b'){

      if(!flag_b)
         b[rc].fill_Random();

      vector< DArray<4> > R(Lx - 1);

      //first construct rightmost operator
      DArray<7> tmp7;
      Contract(1.0,b[rc - 1][Lx - 1],shape(1),peps(rc,Lx - 1),shape(3),0.0,tmp7);

      DArray<8> tmp8;
      Contract(1.0,tmp7,shape(1,5),peps(rc,Lx - 1),shape(3,2),0.0,tmp8);

      DArray<8> tmp8bis;
      Contract(1.0,tmp8,shape(3,6),b[rc][Lx - 1],shape(1,2),0.0,tmp8bis);

      R[Lx - 2] = tmp8bis.reshape_clear(shape(tmp8bis.shape(0),tmp8bis.shape(2),tmp8bis.shape(4),tmp8bis.shape(6)));

      //now move from right to left to construct the rest
      for(int i = Lx - 2;i > 0;--i){

         DArray<6> tmp6;
         Contract(1.0,b[rc - 1][i],shape(3),R[i],shape(0),0.0,tmp6);

         tmp7.clear();
         Contract(1.0,tmp6,shape(1,3),peps(rc,i),shape(3,4),0.0,tmp7);

         tmp6.clear();
         Contract(1.0,tmp7,shape(1,2,6),peps(rc,i),shape(3,4,2),0.0,tmp6);

         Contract(1.0,tmp6,shape(3,5,1),b[rc][i],shape(1,2,3),0.0,R[i - 1]);

      }
      int iter = 0;

      while(iter < comp_sweeps){

         //now start sweeping to get the compressed boundary MPO
         DArray<6> tmp6;
         Contract(1.0,b[rc - 1][0],shape(3),R[0],shape(0),0.0,tmp6);

         tmp7.clear();
         Contract(1.0,peps(rc,0),shape(3,4),tmp6,shape(2,4),0.0,tmp7);

         tmp6.clear();
         Contract(1.0,peps(rc,0),shape(2,3,4),tmp7,shape(2,4,5),0.0,tmp6);

         b[rc][0] = tmp6.reshape_clear(shape(1,D,D,tmp6.shape(5)));

         //QR
         DArray<2> tmp2;
         Geqrf(b[rc][0],tmp2);

         //construct new left operator
         tmp7.clear();
         Contract(1.0,b[rc-1][0],shape(1),peps(rc,0),shape(3),0.0,tmp7);

         tmp8.clear();
         Contract(1.0,tmp7,shape(1,5),peps(rc,0),shape(3,2),0.0,tmp8);

         tmp8bis.clear();
         Contract(1.0,tmp8,shape(3,6),b[rc][0],shape(1,2),0.0,tmp8bis);

         R[0] = tmp8bis.reshape_clear(shape(tmp8bis.shape(1),tmp8bis.shape(3),tmp8bis.shape(5),tmp8bis.shape(7)));

         //now for the rest of the rightgoing sweep.
         for(int i = 1;i < Lx-1;++i){

            tmp6.clear();
            Contract(1.0,R[i - 1],shape(0),b[rc - 1][i],shape(0),0.0,tmp6);

            tmp7.clear();
            Contract(1.0,tmp6,shape(0,3),peps(rc,i),shape(0,3),0.0,tmp7);

            tmp6.clear();
            Contract(1.0,tmp7,shape(0,2,5),peps(rc,i),shape(0,3,2),0.0,tmp6);

            DArray<6> tmp6bis;
            Permute(tmp6,shape(0,2,4,1,3,5),tmp6bis);

            Gemm(CblasNoTrans,CblasNoTrans,1.0,tmp6bis,R[i],0.0,b[rc][i]);

            //QR
            tmp2.clear();
            Geqrf(b[rc][i],tmp2);

            //construct new left operator
            Gemm(CblasTrans,CblasNoTrans,1.0,tmp6bis,b[rc][i],0.0,R[i]);

         }

         //rightmost site
         tmp6.clear();
         Contract(1.0,R[Lx - 2],shape(0),b[rc - 1][Lx - 1],shape(0),0.0,tmp6);

         tmp7.clear();
         Contract(1.0,tmp6,shape(0,3),peps(rc,Lx - 1),shape(0,3),0.0,tmp7);

         tmp6.clear();
         Contract(1.0,tmp7,shape(0,2,5),peps(rc,Lx - 1),shape(0,3,2),0.0,tmp6);

         b[rc][Lx - 1] = tmp6.reshape_clear(shape(tmp6.shape(0),D,D,1));

         //LQ
         tmp2.clear();
         Gelqf(tmp2,b[rc][Lx - 1]);

         //construct new right operator
         tmp7.clear();
         Contract(1.0,b[rc - 1][Lx - 1],shape(1),peps(rc,Lx - 1),shape(3),0.0,tmp7);

         tmp8.clear();
         Contract(1.0,tmp7,shape(1,5),peps(rc,Lx - 1),shape(3,2),0.0,tmp8);

         tmp8bis.clear();
         Contract(1.0,tmp8,shape(3,6),b[rc][Lx - 1],shape(1,2),0.0,tmp8bis);

         R[Lx - 2] = tmp8bis.reshape_clear(shape(tmp8bis.shape(0),tmp8bis.shape(2),tmp8bis.shape(4),tmp8bis.shape(6)));

         //back to the beginning with a leftgoing sweep
         for(int i = Lx-2;i > 0;--i){

            tmp6.clear();
            Contract(1.0,b[rc - 1][i],shape(3),R[i],shape(0),0.0,tmp6);

            tmp7.clear();
            Contract(1.0,tmp6,shape(1,3),peps(rc,i),shape(3,4),0.0,tmp7);

            tmp6.clear();
            Contract(1.0,tmp7,shape(1,2,6),peps(rc,i),shape(3,4,2),0.0,tmp6);

            DArray<6> tmp6bis;
            Permute(tmp6,shape(0,2,4,3,5,1),tmp6bis);

            Gemm(CblasTrans,CblasNoTrans,1.0,R[i - 1],tmp6bis,0.0,b[rc][i]);

            //LQ
            tmp2.clear();
            Gelqf(tmp2,b[rc][i]);

            //construct right operator
            Gemm(CblasNoTrans,CblasTrans,1.0,tmp6bis,b[rc][i],0.0,R[i - 1]);

         }

         //multiply the last L matrix with the first matrix:
         DArray<4> tmp4;
         Gemm(CblasNoTrans,CblasNoTrans,1.0,b[rc][0],tmp2,0.0,tmp4);

         b[rc][0] = std::move(tmp4);

         ++iter;

      }

      //redistribute the norm over the chain
      double nrm =  Nrm2(b[rc][0]);

      //rescale the first site
      Scal((1.0/nrm), b[rc][0]);

      //then multiply the norm over the whole chain
      b[rc].scal(nrm);

   }
   else if(option == 't'){

      if(!flag_t)
         t[rc].fill_Random();

      vector< DArray<4> > R(Lx - 1);

      //first construct rightmost operator
      DArray<7> tmp7;
      Contract(1.0,t[rc + 1][Lx - 1],shape(1),peps(rc+1,Lx - 1),shape(1),0.0,tmp7);

      DArray<8> tmp8;
      Contract(1.0,tmp7,shape(1,4),peps(rc+1,Lx - 1),shape(1,2),0.0,tmp8);

      DArray<8> tmp8bis;
      Contract(1.0,tmp8,shape(3,6),t[rc][Lx - 1],shape(1,2),0.0,tmp8bis);

      R[Lx - 2] = tmp8bis.reshape_clear(shape(tmp8bis.shape(0),tmp8bis.shape(2),tmp8bis.shape(4),tmp8bis.shape(6)));

      //now move from right to left to construct the rest
      for(int i = Lx - 2;i > 0;--i){

         DArray<6> tmp6;
         Contract(1.0,t[rc + 1][i],shape(3),R[i],shape(0),0.0,tmp6);

         tmp7.clear();
         Contract(1.0,tmp6,shape(1,3),peps(rc+1,i),shape(1,4),0.0,tmp7);

         tmp6.clear();
         Contract(1.0,tmp7,shape(1,2,5),peps(rc+1,i),shape(1,4,2),0.0,tmp6);

         Contract(1.0,tmp6,shape(3,5,1),t[rc][i],shape(1,2,3),0.0,R[i - 1]);

      }

      int iter = 0;

      while(iter < comp_sweeps){

         //now start sweeping to get the compressed boundary MPO
         DArray<6> tmp6;
         Contract(1.0,t[rc + 1][0],shape(3),R[0],shape(0),0.0,tmp6);

         tmp7.clear();
         Contract(1.0,peps(rc+1,0),shape(1,4),tmp6,shape(2,4),0.0,tmp7);

         tmp6.clear();
         Contract(1.0,peps(rc+1,0),shape(1,2,4),tmp7,shape(4,1,5),0.0,tmp6);

         t[rc][0] = tmp6.reshape_clear(shape(1,D,D,tmp6.shape(5)));

         //QR
         DArray<2> tmp2;
         Geqrf(t[rc][0],tmp2);

         //construct new left operator
         tmp7.clear();
         Contract(1.0,t[rc+1][0],shape(1),peps(rc+1,0),shape(1),0.0,tmp7);

         tmp8.clear();
         Contract(1.0,tmp7,shape(1,4),peps(rc+1,0),shape(1,2),0.0,tmp8);

         tmp8bis.clear();
         Contract(1.0,tmp8,shape(3,6),t[rc][0],shape(1,2),0.0,tmp8bis);

         R[0] = tmp8bis.reshape_clear(shape(tmp8bis.shape(1),tmp8bis.shape(3),tmp8bis.shape(5),tmp8bis.shape(7)));

         //now for the rest of the rightgoing sweep.
         for(int i = 1;i < Lx-1;++i){

            tmp6.clear();
            Contract(1.0,R[i - 1],shape(0),t[rc + 1][i],shape(0),0.0,tmp6);

            tmp7.clear();
            Contract(1.0,tmp6,shape(0,3),peps(rc+1,i),shape(0,1),0.0,tmp7);

            tmp6.clear();
            Contract(1.0,tmp7,shape(0,2,4),peps(rc+1,i),shape(0,1,2),0.0,tmp6);

            DArray<6> tmp6bis;
            Permute(tmp6,shape(0,2,4,1,3,5),tmp6bis);

            Gemm(CblasNoTrans,CblasNoTrans,1.0,tmp6bis,R[i],0.0,t[rc][i]);

            //QR
            tmp2.clear();
            Geqrf(t[rc][i],tmp2);

            //construct left operator
            Gemm(CblasTrans,CblasNoTrans,1.0,tmp6bis,t[rc][i],0.0,R[i]);

         }

         //rightmost site
         tmp6.clear();
         Contract(1.0,R[Lx - 2],shape(0),t[rc + 1][Lx - 1],shape(0),0.0,tmp6);

         tmp7.clear();
         Contract(1.0,tmp6,shape(0,3),peps(rc+1,Lx - 1),shape(0,1),0.0,tmp7);

         tmp6.clear();
         Contract(1.0,tmp7,shape(0,2,4),peps(rc+1,Lx - 1),shape(0,1,2),0.0,tmp6);

         t[rc][Lx - 1] = tmp6.reshape_clear(shape(tmp6.shape(0),D,D,1));

         //LQ
         tmp2.clear();
         Gelqf(tmp2,t[rc][Lx - 1]);

         //construct new right operator
         tmp7.clear();
         Contract(1.0,t[rc + 1][Lx - 1],shape(1),peps(rc+1,Lx - 1),shape(1),0.0,tmp7);

         tmp8.clear();
         Contract(1.0,tmp7,shape(1,4),peps(rc+1,Lx - 1),shape(1,2),0.0,tmp8);

         tmp8bis.clear();
         Contract(1.0,tmp8,shape(3,6),t[rc][Lx - 1],shape(1,2),0.0,tmp8bis);

         R[Lx - 2] = tmp8bis.reshape_clear(shape(tmp8bis.shape(0),tmp8bis.shape(2),tmp8bis.shape(4),tmp8bis.shape(6)));

         //back to the beginning with a leftgoing sweep
         for(int i = Lx-2;i > 0;--i){

            tmp6.clear();
            Contract(1.0,t[rc + 1][i],shape(3),R[i],shape(0),0.0,tmp6);

            tmp7.clear();
            Contract(1.0,tmp6,shape(1,3),peps(rc+1,i),shape(1,4),0.0,tmp7);

            tmp6.clear();
            Contract(1.0,tmp7,shape(1,5,2),peps(rc+1,i),shape(1,2,4),0.0,tmp6);

            DArray<6> tmp6bis;
            Permute(tmp6,shape(0,2,4,3,5,1),tmp6bis);

            Gemm(CblasTrans,CblasNoTrans,1.0,R[i - 1],tmp6bis,0.0,t[rc][i]);

            //LQ
            tmp2.clear();
            Gelqf(tmp2,t[rc][i]);

            //construct right operator
            Gemm(CblasNoTrans,CblasTrans,1.0,tmp6bis,t[rc][i],0.0,R[i - 1]);

         }

         //multiply the last L matrix with the first matrix:
         DArray<4> tmp4;
         Gemm(CblasNoTrans,CblasNoTrans,1.0,t[rc][0],tmp2,0.0,tmp4);

         t[rc][0] = std::move(tmp4);

         ++iter;

      }

      //redistribute the norm over the chain
      double nrm =  Nrm2(t[rc][0]);

      //rescale the first site
      Scal((1.0/nrm), t[rc][0]);

      //then multiply the norm over the whole chain
      t[rc].scal(nrm);

   }
   else if(option == 'r'){

      if(!flag_r)
         r[rc].fill_Random();

      vector< DArray<4> > R(Ly - 1);

      //first construct rightmost operator
      DArray<7> tmp7;
      Contract(1.0,r[rc + 1][Ly - 1],shape(1),peps(Ly - 1,rc+1),shape(4),0.0,tmp7);

      DArray<8> tmp8;
      Contract(1.0,tmp7,shape(1,5),peps(Ly-1,rc+1),shape(4,2),0.0,tmp8);

      DArray<8> tmp8bis;
      Contract(1.0,tmp8,shape(2,5),r[rc][Ly - 1],shape(1,2),0.0,tmp8bis);

      R[Ly - 2] = tmp8bis.reshape_clear(shape(tmp8bis.shape(0),tmp8bis.shape(3),tmp8bis.shape(5),tmp8bis.shape(6)));

      //now move from right to left to construct the rest
      for(int i = Ly - 2;i > 0;--i){

         DArray<6> tmp6;
         Contract(1.0,r[rc + 1][i],shape(3),R[i],shape(0),0.0,tmp6);

         tmp7.clear();
         Contract(1.0,tmp6,shape(1,3),peps(i,rc+1),shape(4,1),0.0,tmp7);

         tmp6.clear();
         Contract(1.0,tmp7,shape(1,2,5),peps(i,rc+1),shape(4,1,2),0.0,tmp6);

         Contract(1.0,tmp6,shape(2,4,1),r[rc][i],shape(1,2,3),0.0,R[i - 1]);

      }

      int iter = 0;

      while(iter < comp_sweeps){

         //now start sweeping to get the compressed boundary MPO
         DArray<6> tmp6;
         Contract(1.0,r[rc + 1][0],shape(3),R[0],shape(0),0.0,tmp6);

         tmp7.clear();
         Contract(1.0,peps(0,rc+1),shape(1,4),tmp6,shape(4,2),0.0,tmp7);

         tmp6.clear();
         Contract(1.0,peps(0,rc+1),shape(2,4,1),tmp7,shape(1,4,5),0.0,tmp6);

         r[rc][0] = tmp6.reshape_clear(shape(1,D,D,tmp6.shape(5)));

         //QR
         DArray<2> tmp2;
         Geqrf(r[rc][0],tmp2);

         //construct new left operator
         tmp7.clear();
         Contract(1.0,r[rc+1][0],shape(1),peps(0,rc+1),shape(4),0.0,tmp7);

         tmp8.clear();
         Contract(1.0,tmp7,shape(1,5),peps(0,rc+1),shape(4,2),0.0,tmp8);

         tmp8bis.clear();
         Contract(1.0,tmp8,shape(2,5),r[rc][0],shape(1,2),0.0,tmp8bis);

         R[0] = tmp8bis.reshape_clear(shape(tmp8bis.shape(1),tmp8bis.shape(2),tmp8bis.shape(4),tmp8bis.shape(7)));

         //now for the rest of the rightgoing sweep.
         for(int i = 1;i < Ly-1;++i){

            tmp6.clear();
            Contract(1.0,R[i - 1],shape(0),r[rc + 1][i],shape(0),0.0,tmp6);

            tmp7.clear();
            Contract(1.0,tmp6,shape(0,3),peps(i,rc+1),shape(3,4),0.0,tmp7);

            tmp6.clear();
            Contract(1.0,tmp7,shape(0,2,6),peps(i,rc+1),shape(3,4,2),0.0,tmp6);

            DArray<6> tmp6bis;
            Permute(tmp6,shape(0,2,4,1,3,5),tmp6bis);

            Gemm(CblasNoTrans,CblasNoTrans,1.0,tmp6bis,R[i],0.0,r[rc][i]);

            //QR
            tmp2.clear();
            Geqrf(r[rc][i],tmp2);

            //construct left operator
            Gemm(CblasTrans,CblasNoTrans,1.0,tmp6bis,r[rc][i],0.0,R[i]);

         }

         //rightmost site
         tmp6.clear();
         Contract(1.0,R[Ly - 2],shape(0),r[rc + 1][Ly - 1],shape(0),0.0,tmp6);

         tmp7.clear();
         Contract(1.0,tmp6,shape(0,3),peps(Ly - 1,rc+1),shape(3,4),0.0,tmp7);

         tmp6.clear();
         Contract(1.0,tmp7,shape(0,2,6),peps(Ly - 1,rc + 1),shape(3,4,2),0.0,tmp6);

         r[rc][Ly - 1] = tmp6.reshape_clear(shape(tmp6.shape(0),D,D,1));

         //LQ
         tmp2.clear();
         Gelqf(tmp2,r[rc][Ly - 1]);

         //construct new right operator
         tmp7.clear();
         Contract(1.0,r[rc + 1][Ly - 1],shape(1),peps(Ly - 1,rc+1),shape(4),0.0,tmp7);

         tmp8.clear();
         Contract(1.0,tmp7,shape(1,5),peps(Ly - 1,rc+1),shape(4,2),0.0,tmp8);

         tmp8bis.clear();
         Contract(1.0,tmp8,shape(2,5),r[rc][Ly - 1],shape(1,2),0.0,tmp8bis);

         R[Ly - 2] = tmp8bis.reshape_clear(shape(tmp8bis.shape(0),tmp8bis.shape(3),tmp8bis.shape(5),tmp8bis.shape(6)));

         //back to the beginning with a leftgoing sweep
         for(int i = Ly-2;i > 0;--i){

            tmp6.clear();
            Contract(1.0,r[rc + 1][i],shape(3),R[i],shape(0),0.0,tmp6);

            tmp7.clear();
            Contract(1.0,tmp6,shape(1,3),peps(i,rc+1),shape(4,1),0.0,tmp7);

            tmp6.clear();
            Contract(1.0,tmp7,shape(1,2,5),peps(i,rc+1),shape(4,1,2),0.0,tmp6);

            DArray<6> tmp6bis;
            Permute(tmp6,shape(0,3,5,2,4,1),tmp6bis);

            Gemm(CblasTrans,CblasNoTrans,1.0,R[i - 1],tmp6bis,0.0,r[rc][i]);

            //LQ
            tmp2.clear();
            Gelqf(tmp2,r[rc][i]);

            //construct right operator
            Gemm(CblasNoTrans,CblasTrans,1.0,tmp6bis,r[rc][i],0.0,R[i - 1]);

         }

         //multiply the last L matrix with the first matrix:
         DArray<4> tmp4;
         Gemm(CblasNoTrans,CblasNoTrans,1.0,r[rc][0],tmp2,0.0,tmp4);

         r[rc][0] = std::move(tmp4);

         ++iter;

      }

      //redistribute the norm over the chain
      double nrm =  Nrm2(r[rc][0]);

      //rescale the first site
      Scal((1.0/nrm), r[rc][0]);

      //then multiply the norm over the whole chain
      r[rc].scal(nrm);

   }
   else{//left

      if(!flag_l)
         l[rc].fill_Random();

      vector< DArray<4> > R(Ly - 1);

      //first construct rightmost operator
      DArray<7> tmp7;
      Contract(1.0,l[rc - 1][Ly - 1],shape(1),peps(Ly-1,rc),shape(0),0.0,tmp7);

      DArray<8> tmp8;
      Contract(1.0,tmp7,shape(1,4),peps(Ly - 1,rc),shape(0,2),0.0,tmp8);

      DArray<8> tmp8bis;
      Contract(1.0,tmp8,shape(4,7),l[rc][Ly - 1],shape(1,2),0.0,tmp8bis);

      R[Ly - 2] = tmp8bis.reshape_clear(shape(tmp8bis.shape(0),tmp8bis.shape(3),tmp8bis.shape(5),tmp8bis.shape(6)));

      //now move from right to left to construct the rest
      for(int i = Ly - 2;i > 0;--i){

         DArray<6> tmp6;
         Contract(1.0,l[rc - 1][i],shape(3),R[i],shape(0),0.0,tmp6);

         tmp7.clear();
         Contract(1.0,tmp6,shape(1,3),peps(i,rc),shape(0,1),0.0,tmp7);

         tmp6.clear();
         Contract(1.0,tmp7,shape(1,2,4),peps(i,rc),shape(0,1,2),0.0,tmp6);

         Contract(1.0,tmp6,shape(3,5,1),l[rc][i],shape(1,2,3),0.0,R[i - 1]);

      }

      int iter = 0;

      while(iter < comp_sweeps){

         //now start sweeping to get the compressed boundary MPO
         DArray<6> tmp6;
         Contract(1.0,l[rc - 1][0],shape(3),R[0],shape(0),0.0,tmp6);

         tmp7.clear();
         Contract(1.0,peps(0,rc),shape(0,1),tmp6,shape(2,4),0.0,tmp7);

         tmp6.clear();
         Contract(1.0,peps(0,rc),shape(0,1,2),tmp7,shape(4,5,0),0.0,tmp6);

         l[rc][0] = tmp6.reshape_clear(shape(1,D,D,tmp6.shape(5)));

         //QR
         DArray<2> tmp2;
         Geqrf(l[rc][0],tmp2);

         //construct new left operator
         tmp7.clear();
         Contract(1.0,l[rc-1][0],shape(1),peps(0,rc),shape(0),0.0,tmp7);

         tmp8.clear();
         Contract(1.0,tmp7,shape(1,4),peps(0,rc),shape(0,2),0.0,tmp8);

         tmp8bis.clear();
         Contract(1.0,tmp8,shape(4,7),l[rc][0],shape(1,2),0.0,tmp8bis);

         R[0] = tmp8bis.reshape_clear(shape(tmp8bis.shape(1),tmp8bis.shape(2),tmp8bis.shape(4),tmp8bis.shape(7)));

         //now for the rest of the rightgoing sweep.
         for(int i = 1;i < Ly-1;++i){

            tmp6.clear();
            Contract(1.0,R[i - 1],shape(0),l[rc - 1][i],shape(0),0.0,tmp6);

            tmp7.clear();
            Contract(1.0,tmp6,shape(0,3),peps(i,rc),shape(3,0),0.0,tmp7);

            tmp6.clear();
            Contract(1.0,tmp7,shape(0,2,5),peps(i,rc),shape(3,0,2),0.0,tmp6);

            DArray<6> tmp6bis;
            Permute(tmp6,shape(0,3,5,1,2,4),tmp6bis);

            Gemm(CblasNoTrans,CblasNoTrans,1.0,tmp6bis,R[i],0.0,l[rc][i]);

            //QR
            tmp2.clear();
            Geqrf(l[rc][i],tmp2);

            //construct left operator
            Gemm(CblasTrans,CblasNoTrans,1.0,tmp6bis,l[rc][i],0.0,R[i]);

         }

         //rightmost site
         tmp6.clear();
         Contract(1.0,R[Ly - 2],shape(0),l[rc - 1][Ly - 1],shape(0),0.0,tmp6);

         tmp7.clear();
         Contract(1.0,tmp6,shape(0,3),peps(Ly - 1,rc),shape(3,0),0.0,tmp7);

         tmp6.clear();
         Contract(1.0,tmp7,shape(0,2,5),peps(Ly - 1,rc),shape(3,0,2),0.0,tmp6);

         l[rc][Ly - 1] = tmp6.reshape_clear(shape(tmp6.shape(0),D,D,1));

         //LQ
         tmp2.clear();
         Gelqf(tmp2,l[rc][Ly - 1]);

         //construct new right operator
         tmp7.clear();
         Contract(1.0,l[rc - 1][Ly - 1],shape(1),peps(Ly - 1,rc),shape(0),0.0,tmp7);

         tmp8.clear();
         Contract(1.0,tmp7,shape(1,4),peps(Ly - 1,rc),shape(0,2),0.0,tmp8);

         tmp8bis.clear();
         Contract(1.0,tmp8,shape(4,7),l[rc][Ly - 1],shape(1,2),0.0,tmp8bis);

         R[Ly - 2] = tmp8bis.reshape_clear(shape(tmp8bis.shape(0),tmp8bis.shape(3),tmp8bis.shape(5),tmp8bis.shape(6)));

         //back to the beginning with a leftgoing sweep
         for(int i = Ly-2;i > 0;--i){

            tmp6.clear();
            Contract(1.0,l[rc - 1][i],shape(3),R[i],shape(0),0.0,tmp6);

            tmp7.clear();
            Contract(1.0,tmp6,shape(1,3),peps(i,rc),shape(0,1),0.0,tmp7);

            tmp6.clear();
            Contract(1.0,tmp7,shape(1,2,4),peps(i,rc),shape(0,1,2),0.0,tmp6);

            DArray<6> tmp6bis;
            Permute(tmp6,shape(0,2,4,3,5,1),tmp6bis);

            Gemm(CblasTrans,CblasNoTrans,1.0,R[i - 1],tmp6bis,0.0,l[rc][i]);

            //LQ
            tmp2.clear();
            Gelqf(tmp2,l[rc][i]);

            //construct right operator
            Gemm(CblasNoTrans,CblasTrans,1.0,tmp6bis,l[rc][i],0.0,R[i - 1]);

         }

         //multiply the last L matrix with the first matrix:
         DArray<4> tmp4;
         Gemm(CblasNoTrans,CblasNoTrans,1.0,l[rc][0],tmp2,0.0,tmp4);

         l[rc][0] = std::move(tmp4);

         ++iter;

      }

      //redistribute the norm over the chain
      double nrm =  Nrm2(l[rc][0]);

      //rescale the first site
      Scal((1.0/nrm), l[rc][0]);

      //then multiply the norm over the whole chain
      l[rc].scal(nrm);

   }

}
