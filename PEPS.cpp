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
 * construct an empty PEPS object, note: be sure to initialize the Lattice object before calling the constructor
 */
template<typename T>
PEPS<T>::PEPS() : vector< TArray<T,5> >(Lx * Ly) { }

/**
 * construct constructs a standard PEPS object, note: be sure to initialize the Lattice object before calling the constructor
 * @param D_in cutoff virtual dimension
 */
template<typename T>
PEPS<T>::PEPS(int D_in) : vector< TArray<T,5> >(Lx * Ly) {

   D = D_in;
   
   //corners first

   //r == 0 : c == 0
   (*this)[ 0 ].resize(1,D,d,1,D);

   //r == 0 : c == L - 1
   (*this)[ Lx-1 ].resize(D,D,d,1,1);

   //r == L - 1 : c == 0
   (*this)[ (Ly-1)*Lx ].resize(1,1,d,D,D);

   //r == L - 1 : c == L - 1
   (*this)[ (Ly-1)*Lx + Lx-1 ].resize(D,1,d,D,1);

   //sides:

   //r == 0
   for(int c = 1;c < Lx - 1;++c)
      (*this)[ c ].resize(D,D,d,1,D);

   //r == Ly - 1
   for(int c = 1;c < Lx - 1;++c)
      (*this)[ (Ly-1)*Lx + c ].resize(D,1,d,D,D);

   //c == 0
   for(int r = 1;r < Ly - 1;++r)
      (*this)[ r*Lx ].resize(1,D,d,D,D);

   //c == Lx - 1
   for(int r = 1;r < Ly - 1;++r)
      (*this)[ r*Lx + Lx - 1 ].resize(D,D,d,D,1);

   //the rest is full
   for(int r = 1;r < Ly - 1;++r)
      for(int c = 1;c < Lx - 1;++c)
         (*this)[ r*Lx + c ].resize(D,D,d,D,D);

   //now initialize with random numbers
   for(int r = 0;r < Ly;++r)
      for(int c = 0;c < Lx;++c){

         (*this)[ r*Lx + c ].generate(rgen<T>);

         Normalize((*this)[ r*Lx + c ]);
         Scal((T)D,(*this)[ r*Lx + c ]);

      }

}

/**
 * copy constructor
 */
template<typename T>
PEPS<T>::PEPS(const PEPS<T> &peps_copy) : vector< TArray<T,5> >(peps_copy) {

   D = peps_copy.gD();

}

/**
 * empty destructor
 */
template<typename T>
PEPS<T>::~PEPS(){ }

/**
 * access to the individual tensors: const version
 * @param r row index
 * @param c col index
 * @return the tensor on site (r,c)
 */
template<typename T>
const TArray<T,5> &PEPS<T>::operator()(int r,int c) const {

   return (*this)[r*Lx + c];

}

/**
 * access to the individual tensors: const version
 * @param r row index
 * @param c col index
 * @return the tensor on site (r,c)
 */
template<typename T>
TArray<T,5> &PEPS<T>::operator()(int r,int c) {

   return (*this)[r*Lx + c];

}

/**
 * @return the cutoff virtual dimension
 */
template<typename T>
int PEPS<T>::gD() const {

   return D;

}

/**
 * @param D_in value to the D to
 */
template<typename T>
void PEPS<T>::sD(int D_in) {

   this->D = D_in;

}

/**
 * rescale all the tensors, set largest value to one
 */
template<>
void PEPS<double>::rescale_tensors(){

   for(int row = 0;row < Ly;++row)
      for(int col = 0;col < Lx;++col)
         (*this)(row,col).rescale();

}

/**
 * @param row the row index...
 * rescale all the tensors, set largest value on row 'row' to one
 */
template<>
void PEPS<double>::rescale_tensors(int row){

   for(int col = 0;col < Lx;++col)
      (*this)(row,col).rescale();

}

/**
 * @param D_in value to the D to
 */
template<typename T>
void PEPS<T>::fill_Random(){

   //now initialize with random numbers
   for(int r = 0;r < Ly;++r)
      for(int c = 0;c < Lx;++c){

         (*this)[ r*Lx + c ].generate(rgen<T>);

         Normalize((*this)[ r*Lx + c ]);
         Scal((T)D,(*this)[ r*Lx + c ]);

      }

}

/**
 * initialize the peps a completely spin polarized state
 * @param D_in input D
 * @param option up or down spin
 * @param noise level of noise to add
 */
template<>
void PEPS<double>::initialize_ising(int D_in,int option,double noise) {

   D = D_in;

   //bottom row, first site
   (*this)[0].resize(1,D,d,1,D);
   (*this)[0] = 0.0;

   (*this)[0](0,0,option,0,0) = 1.0;

   //bottom row, middle sites
   for(int col = 1;col < Lx - 1;++col){

      (*this)[col].resize(D,D,d,1,D);
      (*this)[col] = 0.0;

      (*this)[col](0,0,option,0,0) = 1.0;

   }

   //bottom row, last site
   (*this)[Lx-1].resize(D,D,d,1,1);
   (*this)[Lx-1] = 0.0;

   (*this)[Lx-1](0,0,option,0,0) = 1.0;

   //middle rows
   for(int row = 1;row < Ly - 1;++row){

      //leftmost middle site: col == 0
      (*this)[row*Lx].resize(1,D,d,D,D);
      (*this)[row*Lx] = 0.0;

      (*this)[row*Lx](0,0,option,0,0) = 1.0;

      //middle sites on row 'row'
      for(int col = 1;col < Lx - 1;++col){

         (*this)[row*Lx + col].resize(D,D,d,D,D);
         (*this)[row*Lx + col] = 0.0;

         (*this)[row*Lx + col](0,0,option,0,0) = 1.0;

      }

      //rightmost site on row 'row'
      (*this)[row*Lx + Lx - 1].resize(D,D,d,D,1);
      (*this)[row*Lx + Lx - 1] = 0.0;

      (*this)[row*Lx + Lx - 1](0,0,option,0,0) = 1.0;

   }

   //top row

   //leftmost site
   (*this)[(Ly - 1)*Lx].resize(1,1,d,D,D);
   (*this)[(Ly - 1)*Lx] = 0.0;

   (*this)[(Ly - 1)*Lx](0,0,option,0,0) = 1.0;

   //top row, middle sites
   for(int col = 1;col < Lx - 1;++col){

      (*this)[(Ly - 1)*Lx + col].resize(D,1,d,D,D);
      (*this)[(Ly - 1)*Lx + col] = 0.0;

      (*this)[(Ly - 1)*Lx + col](0,0,option,0,0) = 1.0;

   }

   //top row rightmost site
   (*this)[(Ly - 1)*Lx + Lx - 1].resize(D,1,d,D,1);
   (*this)[(Ly - 1)*Lx + Lx - 1] = 0.0;

   (*this)[(Ly - 1)*Lx + Lx - 1](0,0,option,0,0) = 1.0;

   //add some noise
   for(int row = 0;row < Lx;++row)
      for(int col = 0;col < Ly;++col){

         IVector<5> dim = (*this)[row*Lx + col].shape();

         for(int i = 0;i < dim[0];++i)
            for(int j = 0;j < dim[1];++j)
               for(int k = 0;k < dim[2];++k)
                  for(int l = 0;l < dim[3];++l)
                     for(int m = 0;m < dim[4];++m)
                        (*this)[row*Lx + col](i,j,k,l,m) += noise * rgen<double>();

      }

}

/**
 * initialize the peps to the direct sum of two antiferromagnetic D=1 structures
 * @param f jastrow factor
 */
template<>
void PEPS<double>::initialize_jastrow(double f) {

   enum {i,j,k,l,m,n,o,p,q,r,s};

   D = 2;

   //bottom row, first site
   (*this)[0].resize(1,D,d,1,D);
   (*this)[0] = 0.0;

   (*this)[0](0,0,0,0,0) = 1.0;
   (*this)[0](0,1,1,0,1) = 1.0;

   //bottom row, middle sites
   for(int col = 1;col < Lx - 1;++col){

      (*this)[col].resize(D,D,d,1,D);
      (*this)[col] = 0.0;

      (*this)[col](0,0,0,0,0) = f;
      (*this)[col](1,0,0,0,0) = 1.0;
      (*this)[col](0,1,1,0,1) = 1.0;
      (*this)[col](1,1,1,0,1) = f;

   }

   //bottom row, last site
   (*this)[Lx-1].resize(D,D,d,1,1);
   (*this)[Lx-1] = 0.0;

   (*this)[Lx-1](0,0,0,0,0) = f;
   (*this)[Lx-1](1,0,0,0,0) = 1.0;

   (*this)[Lx-1](0,1,1,0,0) = 1.0;
   (*this)[Lx-1](1,1,1,0,0) = f;

   //middle sites
   for(int row = 1;row < Ly - 1;++row){

      //leftmost middle site: col == 0
      (*this)[row*Lx].resize(1,D,d,D,D);
      (*this)[row*Lx] = 0.0;

      (*this)[row*Lx](0,0,0,0,0) = f;
      (*this)[row*Lx](0,0,0,1,0) = 1.0;
      (*this)[row*Lx](0,1,1,0,1) = 1.0;
      (*this)[row*Lx](0,1,1,1,1) = f;

      //middle sites on row 'row'
      for(int col = 1;col < Lx - 1;++col){

         (*this)[row*Lx + col].resize(D,D,d,D,D);
         (*this)[row*Lx + col] = 0.0;

         (*this)[row*Lx + col](0,0,0,0,0) = f*f;
         (*this)[row*Lx + col](0,0,0,1,0) = f;
         (*this)[row*Lx + col](1,0,0,0,0) = f;
         (*this)[row*Lx + col](1,0,0,1,0) = 1.0;

         (*this)[row*Lx + col](0,1,1,0,1) = 1.0;
         (*this)[row*Lx + col](0,1,1,1,1) = f;
         (*this)[row*Lx + col](1,1,1,0,1) = f;
         (*this)[row*Lx + col](1,1,1,1,1) = f*f;

      }

      //rightmost site on row 'row'
      (*this)[row*Lx + Lx - 1].resize(D,D,d,D,1);
      (*this)[row*Lx + Lx - 1] = 0.0;

      (*this)[row*Lx + Lx - 1](0,0,0,0,0) = f*f;
      (*this)[row*Lx + Lx - 1](0,0,0,1,0) = f;
      (*this)[row*Lx + Lx - 1](1,0,0,0,0) = f;
      (*this)[row*Lx + Lx - 1](1,0,0,1,0) = 1.0;

      (*this)[row*Lx + Lx - 1](0,1,1,0,0) = 1.0;
      (*this)[row*Lx + Lx - 1](1,1,1,0,0) = f;
      (*this)[row*Lx + Lx - 1](0,1,1,1,0) = f;
      (*this)[row*Lx + Lx - 1](1,1,1,1,0) = f*f;

   }

   //top row
   //leftmost site
   (*this)[(Ly - 1)*Lx].resize(1,1,d,D,D);
   (*this)[(Ly - 1)*Lx] = 0.0;

   (*this)[(Ly - 1)*Lx](0,0,0,0,0) = f;
   (*this)[(Ly - 1)*Lx](0,0,0,1,0) = 1.0;
   (*this)[(Ly - 1)*Lx](0,0,1,0,1) = 1.0;
   (*this)[(Ly - 1)*Lx](0,0,1,1,1) = f;

   //top row, middle sites
   for(int col = 1;col < Lx - 1;++col){

      (*this)[(Ly - 1)*Lx + col].resize(D,1,d,D,D);
      (*this)[(Ly - 1)*Lx + col] = 0.0;

      (*this)[(Ly - 1)*Lx + col](0,0,0,0,0) = f*f;
      (*this)[(Ly - 1)*Lx + col](0,0,0,1,0) = f;
      (*this)[(Ly - 1)*Lx + col](1,0,0,0,0) = f;
      (*this)[(Ly - 1)*Lx + col](1,0,0,1,0) = 1.0;

      (*this)[(Ly - 1)*Lx + col](0,0,1,0,1) = 1.0;
      (*this)[(Ly - 1)*Lx + col](0,0,1,1,1) = f;
      (*this)[(Ly - 1)*Lx + col](1,0,1,0,1) = f;
      (*this)[(Ly - 1)*Lx + col](1,0,1,1,1) = f*f;

   }

   //top row rightmost site
   (*this)[(Ly - 1)*Lx + Lx - 1].resize(D,1,d,D,1);
   (*this)[(Ly - 1)*Lx + Lx - 1] = 0.0;

   (*this)[(Ly - 1)*Lx + Lx - 1](0,0,0,0,0) = f*f;
   (*this)[(Ly - 1)*Lx + Lx - 1](0,0,0,1,0) = f;
   (*this)[(Ly - 1)*Lx + Lx - 1](1,0,0,0,0) = f;
   (*this)[(Ly - 1)*Lx + Lx - 1](1,0,0,1,0) = 1.0;

   (*this)[(Ly - 1)*Lx + Lx - 1](0,0,1,0,0) = 1.0;
   (*this)[(Ly - 1)*Lx + Lx - 1](0,0,1,1,0) = f;
   (*this)[(Ly - 1)*Lx + Lx - 1](1,0,1,0,0) = f;
   (*this)[(Ly - 1)*Lx + Lx - 1](1,0,1,1,0) = f*f;

}

/**
 * increase the bond dimension by one
 * @param D_in bond dimension to grow to
 * @param noise level of noise added to the initial state
 */
template<>
void PEPS<double>::grow_bond_dimension(int D_in,double noise) {

   D = D_in;

   DArray<5> tmp;

   //bottom row, first site
   tmp.resize(1,D,d,1,D);

   tmp.generate(rgen<double>);
   Scal(noise,tmp);

   for(int i = 0;i < (*this)[0].shape(0);++i)
      for(int j = 0;j < (*this)[0].shape(1);++j)
         for(int k = 0;k < (*this)[0].shape(2);++k)
            for(int l = 0;l < (*this)[0].shape(3);++l)
               for(int m = 0;m < (*this)[0].shape(4);++m)
                  tmp(i,j,k,l,m) += (*this)[0](i,j,k,l,m);

   (*this)[0] = std::move(tmp);

   //bottom row, middle sites
   for(int col = 1;col < Lx - 1;++col){

      tmp.clear();
      tmp.resize(D,D,d,1,D);

      tmp.generate(rgen<double>);
      Scal(noise,tmp);

      for(int i = 0;i < (*this)[col].shape(0);++i)
         for(int j = 0;j < (*this)[col].shape(1);++j)
            for(int k = 0;k < (*this)[col].shape(2);++k)
               for(int l = 0;l < (*this)[col].shape(3);++l)
                  for(int m = 0;m < (*this)[col].shape(4);++m)
                     tmp(i,j,k,l,m) += (*this)[col](i,j,k,l,m);

      (*this)[col] = std::move(tmp);

   }

   //bottom row, last site
   tmp.clear();
   tmp.resize(D,D,d,1,1);

   tmp.generate(rgen<double>);
   Scal(noise,tmp);

   for(int i = 0;i < (*this)[Lx-1].shape(0);++i)
      for(int j = 0;j < (*this)[Lx-1].shape(1);++j)
         for(int k = 0;k < (*this)[Lx-1].shape(2);++k)
            for(int l = 0;l < (*this)[Lx-1].shape(3);++l)
               for(int m = 0;m < (*this)[Lx-1].shape(4);++m)
                  tmp(i,j,k,l,m) += (*this)[Lx - 1](i,j,k,l,m);

   (*this)[Lx-1] = std::move(tmp);

   //middle sites
   for(int row = 1;row < Ly - 1;++row){

      //leftmost middle site: col == 0
      tmp.clear();
      tmp.resize(1,D,d,D,D);

      tmp.generate(rgen<double>);
      Scal(noise,tmp);

      for(int i = 0;i < (*this)[row*Lx].shape(0);++i)
         for(int j = 0;j < (*this)[row*Lx].shape(1);++j)
            for(int k = 0;k < (*this)[row*Lx].shape(2);++k)
               for(int l = 0;l < (*this)[row*Lx].shape(3);++l)
                  for(int m = 0;m < (*this)[row*Lx].shape(4);++m)
                     tmp(i,j,k,l,m) += (*this)[row*Lx](i,j,k,l,m);

      (*this)[row*Lx] = std::move(tmp);

      //middle sites on row 'row'
      for(int col = 1;col < Lx - 1;++col){

         tmp.clear();
         tmp.resize(D,D,d,D,D);

         tmp.generate(rgen<double>);
         Scal(noise,tmp);

         for(int i = 0;i < (*this)[row*Lx + col].shape(0);++i)
            for(int j = 0;j < (*this)[row*Lx + col].shape(1);++j)
               for(int k = 0;k < (*this)[row*Lx + col].shape(2);++k)
                  for(int l = 0;l < (*this)[row*Lx + col].shape(3);++l)
                     for(int m = 0;m < (*this)[row*Lx + col].shape(4);++m)
                        tmp(i,j,k,l,m) += (*this)[row*Lx + col](i,j,k,l,m);

         (*this)[row*Lx + col] = std::move(tmp);

      }

      //rightmost site on row 'row'
      tmp.clear();
      tmp.resize(D,D,d,D,1);

      tmp.generate(rgen<double>);
      Scal(noise,tmp);

      for(int i = 0;i < (*this)[row*Lx + Lx - 1].shape(0);++i)
         for(int j = 0;j < (*this)[row*Lx + Lx - 1].shape(1);++j)
            for(int k = 0;k < (*this)[row*Lx + Lx - 1].shape(2);++k)
               for(int l = 0;l < (*this)[row*Lx + Lx - 1].shape(3);++l)
                  for(int m = 0;m < (*this)[row*Lx + Lx - 1].shape(4);++m)
                     tmp(i,j,k,l,m) += (*this)[row*Lx + Lx - 1](i,j,k,l,m);

      (*this)[row*Lx + Lx - 1] = std::move(tmp);

   }

   //top row
   //leftmost site
   tmp.clear();
   tmp.resize(1,1,d,D,D);

   tmp.generate(rgen<double>);
   Scal(noise,tmp);

   for(int i = 0;i < (*this)[(Ly - 1)*Lx].shape(0);++i)
      for(int j = 0;j < (*this)[(Ly - 1)*Lx].shape(1);++j)
         for(int k = 0;k < (*this)[(Ly - 1)*Lx].shape(2);++k)
            for(int l = 0;l < (*this)[(Ly - 1)*Lx].shape(3);++l)
               for(int m = 0;m < (*this)[(Ly - 1)*Lx].shape(4);++m)
                  tmp(i,j,k,l,m) += (*this)[(Ly - 1)*Lx](i,j,k,l,m);

   (*this)[(Ly - 1)*Lx] = std::move(tmp);

   //top row, middle sites
   for(int col = 1;col < Lx - 1;++col){

      tmp.clear();
      tmp.resize(D,1,d,D,D);

      tmp.generate(rgen<double>);
      Scal(noise,tmp);

      for(int i = 0;i < (*this)[(Ly - 1)*Lx + col].shape(0);++i)
         for(int j = 0;j < (*this)[(Ly - 1)*Lx + col].shape(1);++j)
            for(int k = 0;k < (*this)[(Ly - 1)*Lx + col].shape(2);++k)
               for(int l = 0;l < (*this)[(Ly - 1)*Lx + col].shape(3);++l)
                  for(int m = 0;m < (*this)[(Ly - 1)*Lx + col].shape(4);++m)
                     tmp(i,j,k,l,m) += (*this)[(Ly - 1)*Lx + col](i,j,k,l,m);

      (*this)[(Ly - 1)*Lx + col] = std::move(tmp);

   }

   //top row rightmost site
   tmp.clear();
   tmp.resize(D,1,d,D,1);

   tmp.generate(rgen<double>);
   Scal(noise,tmp);

   for(int i = 0;i < (*this)[(Ly - 1)*Lx + Lx - 1].shape(0);++i)
      for(int j = 0;j < (*this)[(Ly - 1)*Lx + Lx - 1].shape(1);++j)
         for(int k = 0;k < (*this)[(Ly - 1)*Lx + Lx - 1].shape(2);++k)
            for(int l = 0;l < (*this)[(Ly - 1)*Lx + Lx - 1].shape(3);++l)
               for(int m = 0;m < (*this)[(Ly - 1)*Lx + Lx - 1].shape(4);++m)
                  tmp(i,j,k,l,m) += (*this)[(Ly - 1)*Lx + Lx - 1](i,j,k,l,m);

   (*this)[(Ly - 1)*Lx + Lx - 1] = std::move(tmp);

}

/**
 * @param peps_i peps to take the overlap with
 * @param init boolean, if true the top-bottom environment has already been calculated and we just need to do the MPO-MPO contraction
 * @return the inner product of two PEPS <psi1|psi2> 
 */
template<>
double PEPS<double>::dot(PEPS<double> &peps_i,bool init) const {

   int half = Ly/2;

   int b_stop = half+1;
   int t_stop = half;

   if(Ly/2 >= Ly - 3){

      if(Ly == 4){

         b_stop = 1;
         t_stop = 0;

      }
      else if( Ly == 6 ){

         b_stop = 2;
         t_stop = 1;

      }

   }

   if(!init){

      //construct bottom environment until half
      env.gb(0).fill('b',peps_i);

      for(int i = 1;i <= b_stop;++i)
         env.add_layer('b',i,peps_i);

      env.gt(Ly - 3).fill('t',peps_i);

      for(int i = Ly - 4;i >= t_stop;--i)
         env.add_layer('t',i,peps_i);

   }

   return env.gb(b_stop).dot(env.gt(t_stop));

}

/** 
 * normalize the peps approximately, using a contraction with auxiliary dimension
 * @param init if true the environment has already been initialized and does not need to be calculated
 */
template<>
void PEPS<double>::normalize(bool init){

   double val = sqrt(this->dot(*this,init));
   val = pow(val,1.0/(double)this->size());

   //now initialize with random numbers
   for(int r = 0;r < Ly;++r)
      for(int c = 0;c < Lx;++c)
         Scal(1.0/val,(*this)[ r*Lx + c ]);

}

/**
 * scale the peps with a number
 * @param val scalar to be multiplied with the peps
 */
template<typename T>
void PEPS<T>::scal(T val){

   val = pow(val,(T)1.0/(T)this->size());

   //now initialize with random numbers
   for(int r = 0;r < Ly;++r)
      for(int c = 0;c < Lx;++c)
         Scal(val,(*this)[ r*Lx + c ]);

}

/**
 * @param mpx will be written to file
 * @param filename name of the file
 * save the MPX object to a file in binary format.
 */
template<typename T>
void PEPS<T>::save(const char *filename){

   for(int row = 0;row < Ly;++row)
      for(int col = 0;col < Lx;++col){

         char name[200];

         sprintf(name,"%s/site_(%d,%d).peps",filename,row,col);

         std::ofstream fout(name);
         fout.precision(16);

         int Da = (*this)(row,col).shape(0);
         int Db = (*this)(row,col).shape(1);
         int Dc = (*this)(row,col).shape(2);
         int Dd = (*this)(row,col).shape(3);
         int De = (*this)(row,col).shape(4);

         fout << Da << "\t" << Db << "\t" << Dc << "\t" << Dd << "\t" << De << endl;

         for(int a_ = 0;a_ < Da;++a_)
            for(int b_ = 0;b_ < Db;++b_)
               for(int c_ = 0;c_ < Dc;++c_)
                  for(int d_ = 0;d_ < Dd;++d_)
                     for(int e_ = 0;e_ < De;++e_)
                        fout << std::scientific << a_ << "\t" << b_ << "\t" << c_ << "\t" << d_ << "\t" << e_ << "\t" << (*this)(row,col)(a_,b_,c_,d_,e_) << endl;

      }

}

/**
 * @param mpx will be constructed from file
 * @param filename name of the file
 * load the MPX object from a file in binary format.
 */
template<typename T>
void PEPS<T>::load(const char *filename){

   for(int row = 0;row < Ly;++row)
      for(int col = 0;col < Lx;++col){

         char name[200];

         sprintf(name,"%s/site_(%d,%d).peps",filename,row,col);

         std::ifstream fin(name);

         int Da,Db,Dc,Dd,De;

         fin >> Da >> Db >> Dc >> Dd >> De;

         (*this)(row,col).resize(Da,Db,Dc,Dd,De);

         for(int a_ = 0;a_ < Da;++a_)
            for(int b_ = 0;b_ < Db;++b_)
               for(int c_ = 0;c_ < Dc;++c_)
                  for(int d_ = 0;d_ < Dd;++d_)
                     for(int e_ = 0;e_ < De;++e_)
                        fin >> a_ >> b_ >> c_ >> d_ >> e_ >> (*this)(row,col)(a_,b_,c_,d_,e_);

      }

}

/**
 * evaluate the expectation value of the energy for the nn-Heisenberg model
 * beware, the environments have to be constructed beforehand!
 */
template<>
double PEPS<double>::energy(){

   // ---- || evaluate the energy in an MPO/MPS manner, first from bottom to top, then left to right || ----
   int delta = ham.gdelta();

   // #################################################################
   // ### ---- from bottom to top: contract in mps/mpo fashion ---- ### 
   // #################################################################

   // -- (1) -- || bottom row: similar to overlap calculation

   //first construct the right renormalized operators
   vector< DArray<5> > R(Lx);
   contractions::init_ro('b',*this,R); 

   //left going unity
   DArray<5> L(1,1,1,1,1);
   L = 1.0;

   //left going operators: Li
   std::vector< DArray<5> > Li_u( delta ); //upper site with extra operator
   std::vector< DArray<5> > Li_d( delta ); //lower site with extra operator

   //some storage stuff:
   DArray<4> tmp4;
   DArray<4> tmp4bis;

   DArray<5> tmp5;
   DArray<5> tmp5bis;
   DArray<5> perm5;

   DArray<6> tmp6;
   DArray<6> tmp6bis;
   DArray<6> perm6;

   DArray<7> tmp7;
   DArray<7> tmp7bis;
   DArray<7> perm7;

   DArray<8> tmp8;
   DArray<8> tmp8bis;
   DArray<8> perm8;

   DArray<9> tmp9;
   DArray<9> tmp9bis;
   DArray<9> perm9;

   int M,N,K;

   enum {j,k,l,m,n,o};

   //peps contracted with a local operator
   DArray<5> peps_op;

   //energy will be stored here
   double val = 0.0;

   //connect top to upper peps: (1,0)
   Gemm(CblasTrans,CblasNoTrans,1.0,env.gt(0)[0],(*this)(1,0),0.0,tmp5);

   Permute(tmp5,shape(1,4,3,0,2),tmp5bis);

   //construct left hamiltonian operators, first upper sites
   for(int i = 0;i < delta;++i){

      //contract peps with the operator
      peps_op.clear();
      Contract(1.0,ham.gL(i),shape(j,k),(*this)(1,0),shape(l,m,k,n,o),0.0,peps_op,shape(l,m,j,n,o));

      //paste on the upper site
      M = tmp5bis.shape(0) * tmp5bis.shape(1) * tmp5bis.shape(2);
      N = peps_op.shape(3) * peps_op.shape(4);
      K = tmp5bis.shape(3) * tmp5bis.shape(4);

      tmp6.resize( shape(tmp5bis.shape(0),tmp5bis.shape(1),tmp5bis.shape(2),1,D,D) );
      blas::gemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, M, N, K, 1.0,tmp5bis.data(),K,peps_op.data(),N,0.0,tmp6.data(),N);

      tmp6bis.clear();
      Permute(tmp6,shape(0,1,5,3,2,4),tmp6bis);

      //contract with bottom environment for Li_u,left going operators
      M = tmp6bis.shape(0) * tmp6bis.shape(1) * tmp6bis.shape(2);
      N = env.gb(0)[0].shape(3);
      K = env.gb(0)[0].shape(0) * env.gb(0)[0].shape(1) * env.gb(0)[0].shape(2);

      Li_u[i].resize( shape(tmp6bis.shape(0),tmp6bis.shape(1),tmp6bis.shape(2),D,D) );
      blas::gemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, M, N, K, 1.0,tmp6bis.data(),K,env.gb(0)[0].data(),N,0.0,Li_u[i].data(),N);

      //vertical energy contribution: first peps(0,0) with operator to tmp6bis
      peps_op.clear();
      Contract(1.0,ham.gR(i),shape(j,k),(*this)(0,0),shape(l,m,k,n,o),0.0,peps_op,shape(l,m,j,n,o));

      M = tmp6bis.shape(0) * tmp6bis.shape(1) * tmp6bis.shape(2) * tmp6bis.shape(4);
      N = (*this)(0,0).shape(2) * (*this)(0,0).shape(4);
      K = (*this)(0,0).shape(1);

      tmp6.resize( shape(tmp6bis.shape(0),tmp6bis.shape(1),tmp6bis.shape(2),tmp6bis.shape(4),d,D ) );
      blas::gemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, M, N, K, 1.0,tmp6bis.data(),K,peps_op.data(),N,0.0,tmp6.data(),N);

      tmp6bis.clear();
      Permute(tmp6,shape(0,1,2,5,3,4),tmp6bis);

      //then add on regular peps
      M = tmp6bis.shape(0) * tmp6bis.shape(1) * tmp6bis.shape(2) * tmp6bis.shape(3);
      N = (*this)(0,0).shape(4);
      K = (*this)(0,0).shape(1) * (*this)(0,0).shape(2);

      tmp5.resize( shape(tmp6bis.shape(0),tmp6bis.shape(1),tmp6bis.shape(2),tmp6bis.shape(3),D) );
      blas::gemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, M, N, K, 1.0,tmp6bis.data(),K,(*this)(0,0).data(),N,0.0,tmp5.data(),N);

      //now permute
      Permute(tmp5,shape(0,1,2,4,3),Li_d[i]);

      //and contract with R[0]
      val += ham.gcoef_n(i) * Dot(Li_d[i],R[0]);
      cout << 0 << "\t" << i << "\t" << val << endl;

   }

   //Now do the bottom left operators: paste on the second 'regular' peps
   M = tmp5bis.shape(0) * tmp5bis.shape(1) * tmp5bis.shape(2);
   N = (*this)(1,0).shape(3) * (*this)(1,0).shape(4);
   K = tmp5bis.shape(3) * tmp5bis.shape(4);

   tmp5.resize( shape(tmp5bis.shape(0),tmp5bis.shape(1),tmp5bis.shape(2),D,D) );
   blas::gemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, M, N, K, 1.0,tmp5bis.data(),K,(*this)(1,0).data(),N,0.0,tmp5.data(),N);

   tmp5bis.clear();
   Permute(tmp5,shape(0,1,4,2,3),tmp5bis);

   //add on the different lower opeators
   for(int i = 0;i < delta;++i){

      peps_op.clear();
      Contract(1.0,ham.gL(i),shape(j,k),(*this)(0,0),shape(l,m,k,n,o),0.0,peps_op,shape(l,m,j,n,o));

      M = tmp5bis.shape(0) * tmp5bis.shape(1) * tmp5bis.shape(2) * tmp5bis.shape(3);
      N = (*this)(0,0).shape(2) * (*this)(0,0).shape(4);
      K = (*this)(0,0).shape(1);

      tmp6.resize( shape(tmp5bis.shape(0),tmp5bis.shape(1),tmp5bis.shape(2),tmp5bis.shape(3),d,D ) );
      blas::gemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, M, N, K, 1.0,tmp5bis.data(),K,peps_op.data(),N,0.0,tmp6.data(),N);

      tmp6bis.clear();
      Permute(tmp6,shape(0,1,2,5,3,4),tmp6bis);

      //then add on regular peps
      M = tmp6bis.shape(0) * tmp6bis.shape(1) * tmp6bis.shape(2) * tmp6bis.shape(3);
      N = (*this)(0,0).shape(4);
      K = (*this)(0,0).shape(1) * (*this)(0,0).shape(2);

      tmp5.resize( shape(tmp6bis.shape(0),tmp6bis.shape(1),tmp6bis.shape(2),tmp6bis.shape(3),D) );
      blas::gemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, M, N, K, 1.0,tmp6bis.data(),K,(*this)(0,0).data(),N,0.0,tmp5.data(),N);

      //now permute to get correct left going down
      Permute(tmp5,shape(0,1,2,4,3),Li_d[i]);

   }

   //and finally construct left unity (in R[0]), contract with bottom environment
   M = tmp5bis.shape(0) * tmp5bis.shape(1) * tmp5bis.shape(2);
   N = env.gb(0)[0].shape(3);
   K = env.gb(0)[0].shape(0) * env.gb(0)[0].shape(1) * env.gb(0)[0].shape(2);

   blas::gemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, M, N, K, 1.0,tmp5bis.data(),K,env.gb(0)[0].data(),N,0.0,tmp5.data(),N);

   //now for the middle terms
   for(int col = 1;col < Lx - 1;++col){

      //close down the left renormalized operators from the previous site

      //start with Left Up
      for(int i = 0;i < delta;++i){

         //connect top to left
         tmp7.clear();
         Gemm(CblasTrans,CblasNoTrans,1.0,env.gt(0)[col],Li_u[i],0.0,tmp7);

         tmp7bis.clear();
         Permute(tmp7,shape(1,2,4,5,6,3,0),tmp7bis);

         //add peps
         tmp8.clear();
         Gemm(CblasNoTrans,CblasNoTrans,1.0,tmp7bis,(*this)(1,col),0.0,tmp8);

         tmp8bis.clear();
         Permute(tmp8,shape(1,3,4,6,7,2,0,5),tmp8bis);

         //and yet another peps
         tmp7.clear();
         Gemm(CblasNoTrans,CblasNoTrans,1.0,tmp8bis,(*this)(1,col),0.0,tmp7);

         tmp7bis.clear();
         Permute(tmp7,shape(0,4,6,2,5,1,3),tmp7bis);

         //paste bottom regular peps on
         tmp8.clear();
         Gemm(CblasNoTrans,CblasNoTrans,1.0,tmp7bis,(*this)(0,col),0.0,tmp8);

         tmp8bis.clear();
         Permute(tmp8,shape(0,1,2,7,3,4,5,6),tmp8bis);

         //and construct the left up-down operator
         peps_op.clear();
         Contract(1.0,ham.gR(i),shape(j,k),(*this)(0,col),shape(l,m,k,n,o),0.0,peps_op,shape(l,m,j,n,o));

         Li_u[i].clear();
         Gemm(CblasNoTrans,CblasNoTrans,1.0,tmp8bis,peps_op,0.0,Li_u[i]);

         val += ham.gcoef_nn(i) * Dot(Li_u[i],R[col]);

      }

      //Left down - close down with a diagonal and horizontal term!
      for(int i = 0;i < delta;++i){

         //first add on top to Left renormalized operator
         tmp7.clear();
         Gemm(CblasTrans,CblasNoTrans,1.0,env.gt(0)[col],Li_d[i],0.0,tmp7);

         tmp7bis.clear();
         Permute(tmp7,shape(1,2,4,5,6,3,0),tmp7bis);

         //add peps
         tmp8.clear();
         Gemm(CblasNoTrans,CblasNoTrans,1.0,tmp7bis,(*this)(1,col),0.0,tmp8);

         tmp8bis.clear();
         Permute(tmp8,shape(1,3,4,6,7,2,0,5),tmp8bis);

         //1) do the diagonal link first
         peps_op.clear();
         Contract(1.0,ham.gR(i),shape(j,k),(*this)(1,col),shape(l,m,k,n,o),0.0,peps_op,shape(l,m,j,n,o));

         tmp7.clear();
         Gemm(CblasNoTrans,CblasNoTrans,1.0,tmp8bis,peps_op,0.0,tmp7);

         tmp7bis.clear();
         Permute(tmp7,shape(0,4,6,1,2,3,5),tmp7bis);

         //now contract with bottom environment:
         M = tmp7bis.shape(0) * tmp7bis.shape(1) * tmp7bis.shape(2);
         N = env.gb(0)[col].shape(3);
         K = env.gb(0)[col].shape(0) * env.gb(0)[col].shape(1) * env.gb(0)[col].shape(2);

         Li_d[i].resize( shape(tmp7bis.shape(0),tmp7bis.shape(1),tmp7bis.shape(2),D,D) );
         blas::gemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, M, N, K, 1.0,tmp7bis.data(),K,env.gb(0)[col].data(),N,0.0,Li_d[i].data(),N);

         val += ham.gcoef_nn(i) * Dot(Li_d[i],R[col]);

         //2) then do the horizontal gate:

         //add regular peps
         tmp7.clear();
         Gemm(CblasNoTrans,CblasNoTrans,1.0,tmp8bis,(*this)(1,col),0.0,tmp7);

         tmp7bis.clear();
         Permute(tmp7,shape(0,4,6,2,5,1,3),tmp7bis);

         //paste bottom regular peps on
         tmp8.clear();
         Gemm(CblasNoTrans,CblasNoTrans,1.0,tmp7bis,(*this)(0,col),0.0,tmp8);

         tmp8bis.clear();
         Permute(tmp8,shape(0,1,2,7,3,4,5,6),tmp8bis);

         //construct the left down-down operator (horizontal gate)
         peps_op.clear();
         Contract(1.0,ham.gR(i),shape(j,k),(*this)(0,col),shape(l,m,k,n,o),0.0,peps_op,shape(l,m,j,n,o));

         Gemm(CblasNoTrans,CblasNoTrans,1.0,tmp8bis,peps_op,0.0,Li_d[i]);

         val += ham.gcoef_n(i) * Dot(Li_d[i],R[col]);
         cout << col << "\t" << i << "\t" << val << endl;

      }

      //now construct left going operators and evaluate the vertical term

      //first add top to left unity R[col - 1]:
      tmp7.clear();
      Gemm(CblasTrans,CblasNoTrans,1.0,env.gt(0)[col],R[col - 1],0.0,tmp7);

      tmp7bis.clear();
      Permute(tmp7,shape(1,2,4,5,6,3,0),tmp7bis);

      //add peps
      tmp8.clear();
      Gemm(CblasNoTrans,CblasNoTrans,1.0,tmp7bis,(*this)(1,col),0.0,tmp8);

      perm8.clear();
      Permute(tmp8,shape(1,3,4,6,7,2,0,5),perm8);

      //add upper operator for Lu_i and vertical contribution
      for(int i = 0;i < delta;++i){

         //add left operator on top site
         peps_op.clear();
         Contract(1.0,ham.gL(i),shape(j,k),(*this)(1,col),shape(l,m,k,n,o),0.0,peps_op,shape(l,m,j,n,o));

         tmp7.clear();
         Gemm(CblasNoTrans,CblasNoTrans,1.0,perm8,peps_op,0.0,tmp7);

         tmp7bis.clear();
         Permute(tmp7,shape(0,4,6,2,5,1,3),tmp7bis);

         //paste bottom regular peps on
         tmp8.clear();
         Gemm(CblasNoTrans,CblasNoTrans,1.0,tmp7bis,(*this)(0,col),0.0,tmp8);

         tmp8bis.clear();
         Permute(tmp8,shape(0,1,2,7,3,4,5,6),tmp8bis);

         //1) first make the vertical energy contribution
         peps_op.clear();
         Contract(1.0,ham.gR(i),shape(j,k),(*this)(0,col),shape(l,m,k,n,o),0.0,peps_op,shape(l,m,j,n,o));

         Gemm(CblasNoTrans,CblasNoTrans,1.0,tmp8bis,peps_op,0.0,Li_u[i]);

         val += ham.gcoef_n(i) * Dot(Li_u[i],R[col]);
         cout << col << "\t" << i << "\t" << val << endl;

         //2) the construct the 'real' Lu_i
         Gemm(CblasNoTrans,CblasNoTrans,1.0,tmp8bis,(*this)(0,col),0.0,Li_u[i]);

      }

      //add regular peps to perm8
      tmp7.clear();
      Gemm(CblasNoTrans,CblasNoTrans,1.0,perm8,(*this)(1,col),0.0,tmp7);

      tmp7bis.clear();
      Permute(tmp7,shape(0,4,6,2,5,1,3),tmp7bis);

      //paste bottom regular peps on
      tmp8.clear();
      Gemm(CblasNoTrans,CblasNoTrans,1.0,tmp7bis,(*this)(0,col),0.0,tmp8);

      tmp8bis.clear();
      Permute(tmp8,shape(0,1,2,7,3,4,5,6),tmp8bis);

      //then construct the lower operator Ld_i
      for(int i = 0;i < delta;++i){

         //construct lower left operator
         peps_op.clear();
         Contract(1.0,ham.gL(i),shape(j,k),(*this)(0,col),shape(l,m,k,n,o),0.0,peps_op,shape(l,m,j,n,o));

         Gemm(CblasNoTrans,CblasNoTrans,1.0,tmp8bis,peps_op,0.0,Li_d[i]);

      }

      //finally make new unity by adding one more regular peps
      Gemm(CblasNoTrans,CblasNoTrans,1.0,tmp8bis,(*this)(0,col),0.0,R[col]);

   }

   //last site of bottom row: close down the left operators

   //start with Left Up
   for(int i = 0;i < delta;++i){

      //connect top to left
      tmp7.clear();
      Gemm(CblasTrans,CblasNoTrans,1.0,env.gt(0)[Lx - 1],Li_u[i],0.0,tmp7);

      tmp7bis.clear();
      Permute(tmp7,shape(1,2,4,5,6,3,0),tmp7bis);

      //add peps
      tmp8.clear();
      Gemm(CblasNoTrans,CblasNoTrans,1.0,tmp7bis,(*this)(1,Lx-1),0.0,tmp8);

      tmp8bis.clear();
      Permute(tmp8,shape(1,3,4,6,7,2,0,5),tmp8bis);

      //and yet another peps
      tmp7.clear();
      Gemm(CblasNoTrans,CblasNoTrans,1.0,tmp8bis,(*this)(1,Lx-1),0.0,tmp7);

      tmp7bis.clear();
      Permute(tmp7,shape(0,4,6,2,5,1,3),tmp7bis);

      //paste bottom regular peps on
      tmp8.clear();
      Gemm(CblasNoTrans,CblasNoTrans,1.0,tmp7bis,(*this)(0,Lx-1),0.0,tmp8);

      tmp8bis.clear();
      Permute(tmp8,shape(0,1,2,7,3,4,5,6),tmp8bis);

      //and construct the left up-down operator
      peps_op.clear();
      Contract(1.0,ham.gR(i),shape(j,k),(*this)(0,Lx-1),shape(l,m,k,n,o),0.0,peps_op,shape(l,m,j,n,o));

      val += ham.gcoef_nn(i) * blas::dot(tmp8bis.size(),tmp8bis.data(),1,peps_op.data(),1);

   }

   //Left down - close down with a diagonal and horizontal term!
   for(int i = 0;i < delta;++i){

      //first add on top to Left renormalized operator
      tmp7.clear();
      Gemm(CblasTrans,CblasNoTrans,1.0,env.gt(0)[Lx-1],Li_d[i],0.0,tmp7);

      tmp7bis.clear();
      Permute(tmp7,shape(1,2,4,5,6,3,0),tmp7bis);

      //add peps
      tmp8.clear();
      Gemm(CblasNoTrans,CblasNoTrans,1.0,tmp7bis,(*this)(1,Lx-1),0.0,tmp8);

      tmp8bis.clear();
      Permute(tmp8,shape(1,3,4,6,7,2,0,5),tmp8bis);

      //1) do the diagonal link first
      peps_op.clear();
      Contract(1.0,ham.gR(i),shape(j,k),(*this)(1,Lx-1),shape(l,m,k,n,o),0.0,peps_op,shape(l,m,j,n,o));

      tmp7.clear();
      Gemm(CblasNoTrans,CblasNoTrans,1.0,tmp8bis,peps_op,0.0,tmp7);

      tmp7bis.clear();
      Permute(tmp7,shape(0,4,6,1,2,3,5),tmp7bis);

      //inner product with bottom environment leads to energy contribution:
      val += ham.gcoef_nn(i) * blas::dot(tmp7bis.size(),tmp7bis.data(),1,env.gb(0)[Lx-1].data(),1);

      //2) then do the horizontal gate:

      //add regular peps
      tmp7.clear();
      Gemm(CblasNoTrans,CblasNoTrans,1.0,tmp8bis,(*this)(1,Lx-1),0.0,tmp7);

      tmp7bis.clear();
      Permute(tmp7,shape(0,4,6,2,5,1,3),tmp7bis);

      //paste bottom regular peps on
      tmp8.clear();
      Gemm(CblasNoTrans,CblasNoTrans,1.0,tmp7bis,(*this)(0,Lx-1),0.0,tmp8);

      tmp8bis.clear();
      Permute(tmp8,shape(0,1,2,7,3,4,5,6),tmp8bis);

      //construct the left down-down operator (horizonatal gate)
      peps_op.clear();
      Contract(1.0,ham.gR(i),shape(j,k),(*this)(0,Lx-1),shape(l,m,k,n,o),0.0,peps_op,shape(l,m,j,n,o));

      val += ham.gcoef_n(i) * blas::dot(tmp8bis.size(),tmp8bis.data(),1,peps_op.data(),1);

   }

   //finally the last vertical term:
   tmp7.clear();
   Gemm(CblasTrans,CblasNoTrans,1.0,env.gt(0)[Lx-1],R[Lx - 2],0.0,tmp7);

   tmp7bis.clear();
   Permute(tmp7,shape(1,2,4,5,6,3,0),tmp7bis);

   //add peps
   tmp8.clear();
   Gemm(CblasNoTrans,CblasNoTrans,1.0,tmp7bis,(*this)(1,Lx-1),0.0,tmp8);

   perm8.clear();
   Permute(tmp8,shape(1,3,4,6,7,2,0,5),perm8);

   //add upper operator for Lu_i and vertical contribution
   for(int i = 0;i < delta;++i){

      //add left operator on top site
      peps_op.clear();
      Contract(1.0,ham.gL(i),shape(j,k),(*this)(1,Lx-1),shape(l,m,k,n,o),0.0,peps_op,shape(l,m,j,n,o));

      tmp7.clear();
      Gemm(CblasNoTrans,CblasNoTrans,1.0,perm8,peps_op,0.0,tmp7);

      tmp7bis.clear();
      Permute(tmp7,shape(0,4,6,2,5,1,3),tmp7bis);

      //paste bottom regular peps on
      tmp8.clear();
      Gemm(CblasNoTrans,CblasNoTrans,1.0,tmp7bis,(*this)(0,Lx-1),0.0,tmp8);

      tmp8bis.clear();
      Permute(tmp8,shape(0,1,2,7,3,4,5,6),tmp8bis);

      //make the vertical energy contribution
      peps_op.clear();
      Contract(1.0,ham.gR(i),shape(j,k),(*this)(0,Lx-1),shape(l,m,k,n,o),0.0,peps_op,shape(l,m,j,n,o));

      val += ham.gcoef_n(i) * blas::dot(tmp8bis.size(),tmp8bis.data(),1,peps_op.data(),1);

   }
/*
   // -- (2) -- now move from bottom to top calculating everything like an MPO/MPS expectation value

   //Right renormalized operators
   vector< DArray<6> > RO(Lx);

   //array of delta left upper and lower renormalized operators needed
   std::vector< DArray<6> > LOi_u(delta);
   std::vector< DArray<6> > LOi_d(delta);

   for(int row = 1;row < Ly - 2;++row){

      //first create right renormalized operator
      contractions::init_ro(row,*this,RO);

      // --- now move from left to right to get the expecation value of the interactions ---

      //start with the first site, evaluate the vertical term and the left going terms

      //paste upper row regular peps on top environment
      tmp5.clear();
      Gemm(CblasTrans,CblasNoTrans,1.0,env.gt(row)[0],(*this)(row+1,0),0.0,tmp5);

      perm5.clear();
      Permute(tmp5,shape(1,3,4,0,2),perm5);

      //construct Left-Up operator and vertical energy contribution
      for(int i = 0;i < delta;++i){

         //paste operator on upper row peps
         peps_op.clear();
         Contract(1.0,ham.gL(i),shape(j,k),(*this)(row+1,0),shape(l,m,k,n,o),0.0,peps_op,shape(l,m,j,n,o));

         //and add onto perm5
         M = perm5.shape(0) * perm5.shape(1) * perm5.shape(2);
         N = peps_op.shape(3) * peps_op.shape(4);
         K = peps_op.shape(1) * peps_op.shape(2);

         tmp5.resize( shape(env.gt(row)[0].shape(3),D,D,D,D) );
         blas::gemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, M, N, K, 1.0,perm5.data(),K,peps_op.data(),N,0.0,tmp5.data(),N);

         tmp5bis.clear();
         Permute(tmp5,shape(0,2,4,3,1),tmp5bis);

         //paste lower row regular peps on
         M = tmp5bis.shape(0) * tmp5bis.shape(1) * tmp5bis.shape(2) * tmp5bis.shape(3);
         N = (*this)(row,0).shape(2) * (*this)(row,0).shape(3) * (*this)(row,0).shape(4);
         K = tmp5bis.shape(4);

         tmp7.resize( shape(env.gt(row)[0].shape(3),D,D,D,d,D,D) );
         blas::gemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, M, N, K, 1.0,tmp5bis.data(),K,(*this)(row,0).data(),N,0.0,tmp7.data(),N);

         perm7.clear();
         Permute(tmp7,shape(0,1,2,5,6,3,4),perm7);

         //first evaluate vertical energy contribution
         peps_op.clear();
         Contract(1.0,ham.gR(i),shape(j,k),(*this)(row,0),shape(l,m,k,n,o),0.0,peps_op,shape(l,m,j,n,o));

         //add on lower site peps with operator
         M = perm7.shape(0) * perm7.shape(1) * perm7.shape(2) * perm7.shape(3) * perm7.shape(4);
         N = (*this)(row,0).shape(3) * (*this)(row,0).shape(4);
         K = perm7.shape(5) * perm7.shape(6);

         tmp7.resize( shape(env.gt(row)[0].shape(3),D,D,D,D,D,D) );
         blas::gemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, M, N, K, 1.0,perm7.data(),K,peps_op.data(),N,0.0,tmp7.data(),N);

         tmp7bis.clear();
         Permute(tmp7,shape(0,1,2,4,6,3,5),tmp7bis);

         M = tmp7bis.shape(0) * tmp7bis.shape(1) * tmp7bis.shape(2) * tmp7bis.shape(3) * tmp7bis.shape(4);
         N = env.gb(row-1)[0].shape(3);
         K = tmp7bis.shape(5) * tmp7bis.shape(6);

         LOi_u[i].resize( shape(env.gt(row)[0].shape(3),D,D,D,D,env.gb(row-1)[0].shape(3)) );
         blas::gemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, M, N, K, 1.0,tmp7bis.data(),K,env.gb(row-1)[0].data(),N,0.0,LOi_u[i].data(),N);

         //vertical energy
         val += ham.gcoef_n(i) * Dot(LOi_u[i],RO[0]);

         //then construct leftgoing upper operator

         //add on lower site regular peps
         M = perm7.shape(0) * perm7.shape(1) * perm7.shape(2) * perm7.shape(3) * perm7.shape(4);
         N = (*this)(row,0).shape(3) * (*this)(row,0).shape(4);
         K = perm7.shape(5) * perm7.shape(6);

         tmp7.resize( shape(env.gt(row)[0].shape(3),D,D,D,D,D,D) );
         blas::gemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, M, N, K, 1.0,perm7.data(),K,(*this)(row,0).data(),N,0.0,tmp7.data(),N);

         tmp7bis.clear();
         Permute(tmp7,shape(0,1,2,4,6,3,5),tmp7bis);

         M = tmp7bis.shape(0) * tmp7bis.shape(1) * tmp7bis.shape(2) * tmp7bis.shape(3) * tmp7bis.shape(4);
         N = env.gb(row-1)[0].shape(3);
         K = tmp7bis.shape(5) * tmp7bis.shape(6);

         blas::gemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, M, N, K, 1.0,tmp7bis.data(),K,env.gb(row-1)[0].data(),N,0.0,LOi_u[i].data(),N);

      }

      //next regular upper peps onto perm5 to construct lower left operator and unity
      M = perm5.shape(0) * perm5.shape(1) * perm5.shape(2);
      N = (*this)(row+1,0).shape(3) * (*this)(row+1,0).shape(4);
      K = (*this)(row+1,0).shape(1) * (*this)(row+1,0).shape(2);

      tmp5.resize( shape(env.gt(row)[0].shape(3),D,D,D,D) );
      blas::gemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, M, N, K, 1.0,perm5.data(),K,(*this)(row+1,0).data(),N,0.0,tmp5.data(),N);

      tmp5bis.clear();
      Permute(tmp5,shape(0,2,4,3,1),tmp5bis);

      //paste lower row regular peps on
      M = tmp5bis.shape(0) * tmp5bis.shape(1) * tmp5bis.shape(2) * tmp5bis.shape(3);
      N = (*this)(row,0).shape(2) * (*this)(row,0).shape(3) * (*this)(row,0).shape(4);
      K = tmp5bis.shape(4);

      tmp7.resize( shape(env.gt(row)[0].shape(3),D,D,D,d,D,D) );
      blas::gemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, M, N, K, 1.0,tmp5bis.data(),K,(*this)(row,0).data(),N,0.0,tmp7.data(),N);

      perm7.clear();
      Permute(tmp7,shape(0,1,2,5,6,3,4),perm7);

      //then construct the Left-Down operator
      for(int i = 0;i < delta;++i){

         //add on lower site peps with operator
         peps_op.clear();
         Contract(1.0,ham.gL(i),shape(j,k),(*this)(row,0),shape(l,m,k,n,o),0.0,peps_op,shape(l,m,j,n,o));

         M = perm7.shape(0) * perm7.shape(1) * perm7.shape(2) * perm7.shape(3) * perm7.shape(4);
         N = (*this)(row,0).shape(3) * (*this)(row,0).shape(4);
         K = perm7.shape(5) * perm7.shape(6);

         tmp7.resize( shape(env.gt(row)[0].shape(3),D,D,D,D,D,D) );
         blas::gemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, M, N, K, 1.0,perm7.data(),K,peps_op.data(),N,0.0,tmp7.data(),N);

         tmp7bis.clear();
         Permute(tmp7,shape(0,1,2,4,6,3,5),tmp7bis);

         M = tmp7bis.shape(0) * tmp7bis.shape(1) * tmp7bis.shape(2) * tmp7bis.shape(3) * tmp7bis.shape(4);
         N = env.gb(row-1)[0].shape(3);
         K = tmp7bis.shape(5) * tmp7bis.shape(6);

         LOi_d[i].resize( shape(env.gt(row)[0].shape(3),D,D,D,D,env.gb(row-1)[0].shape(3)) );
         blas::gemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, M, N, K, 1.0,tmp7bis.data(),K,env.gb(row-1)[0].data(),N,0.0,LOi_d[i].data(),N);

      }

      //and finally the left unit operator
      M = perm7.shape(0) * perm7.shape(1) * perm7.shape(2) * perm7.shape(3) * perm7.shape(4);
      N = (*this)(row,0).shape(3) * (*this)(row,0).shape(4);
      K = perm7.shape(5) * perm7.shape(6);

      tmp7.resize( shape(env.gt(row)[0].shape(3),D,D,D,D,D,D) );
      blas::gemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, M, N, K, 1.0,perm7.data(),K,(*this)(row,0).data(),N,0.0,tmp7.data(),N);

      tmp7bis.clear();
      Permute(tmp7,shape(0,1,2,4,6,3,5),tmp7bis);

      M = tmp7bis.shape(0) * tmp7bis.shape(1) * tmp7bis.shape(2) * tmp7bis.shape(3) * tmp7bis.shape(4);
      N = env.gb(row-1)[0].shape(3);
      K = tmp7bis.shape(5) * tmp7bis.shape(6);

      RO[0].resize( shape(env.gt(row)[0].shape(3),D,D,D,D,env.gb(row-1)[0].shape(3)) );
      blas::gemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, M, N, K, 1.0,tmp7bis.data(),K,env.gb(row-1)[0].data(),N,0.0,RO[0].data(),N);

      // --- now for the middle sites, close down the operators on the left and construct new ones --- 
      for(int col = 1;col < Lx - 1;++col){

         //close down the left up and down operators with two diagonal and one horizontal contribution
         tmp8.clear();
         Gemm(CblasNoTrans,CblasTrans,1.0,env.gb(row-1)[col],RO[col],0.0,tmp8);

         tmp8bis.clear();
         Permute(tmp8,shape(2,7,0,1,3,4,5,6),tmp8bis);

         //add regular peps on lower site
         tmp9.clear();
         Gemm(CblasNoTrans,CblasNoTrans,1.0,(*this)(row,col),tmp8bis,0.0,tmp9);

         perm9.clear();
         Permute(tmp9,shape(2,4,8,0,1,3,5,6,7),perm9);

         //first close down horizontal, Lu-Rd diagonal and vertical
         for(int i = 0;i < delta;++i){

            //--- HORIZONTAL AND LU-RD DIAGONAL ---

            //act with operator on peps on lower site: hermitian adjoint because of change in direction when going from right to left
            peps_op.clear();
            Contract(1.0,ham.gR(i),shape(k,j),(*this)(row,col),shape(l,m,k,n,o),0.0,peps_op,shape(l,m,j,n,o));

            tmp8.clear();
            Gemm(CblasNoTrans,CblasNoTrans,1.0,peps_op,perm9,0.0,tmp8);

            tmp8bis.clear();
            Permute(tmp8,shape(3,7,1,6,0,2,4,5),tmp8bis);

            //add regular peps on upper site
            tmp9.clear();
            Gemm(CblasNoTrans,CblasNoTrans,1.0,(*this)(row+1,col),tmp8bis,0.0,tmp9);

            tmp9bis.clear();
            Permute(tmp9,shape(2,3,4,0,1,5,6,7,8),tmp9bis);

            //yet another regular on upper site
            tmp8.clear();
            Gemm(CblasNoTrans,CblasNoTrans,1.0,(*this)(row+1,col),tmp9bis,0.0,tmp8);

            tmp8bis.clear();
            Permute(tmp8,shape(1,3,7,0,2,4,5,6),tmp8bis);

            //finally top environment for closure
            tmp6.clear();
            Gemm(CblasNoTrans,CblasNoTrans,1.0,env.gt(row)[col],tmp8bis,0.0,tmp6);

            //now close down horizontal
            val += ham.gcoef_n(i) * Dot(tmp6,LOi_d[i]);

            //and left-up right-down diagonal
            val += ham.gcoef_nn(i) * Dot(tmp6,LOi_u[i]);

            //--- VERTICAL GATE ---

            //act with operator on peps on upper site, for vertical energy contribution: hermitian adjoint
            peps_op.clear();
            Contract(1.0,ham.gL(i),shape(k,j),(*this)(row+1,col),shape(l,m,k,n,o),0.0,peps_op,shape(l,m,j,n,o));

            //add it on the intermediary with already an operator on the lower site (tmp9bis)
            Gemm(CblasNoTrans,CblasNoTrans,1.0,peps_op,tmp9bis,0.0,tmp8);

            Permute(tmp8,shape(1,3,7,0,2,4,5,6),tmp8bis);

            //finally top environment for closure
            Gemm(CblasNoTrans,CblasNoTrans,1.0,env.gt(row)[col],tmp8bis,0.0,tmp6);

            //close down vertical with left unity
            val += ham.gcoef_n(i) * Dot(tmp6,RO[col-1]);

         }

         //now close down the left-down right-up diagonal
         for(int i = 0;i < delta;++i){

            //add regular lower-site peps to intermediate perm9
            tmp8.clear();
            Gemm(CblasNoTrans,CblasNoTrans,1.0,(*this)(row,col),perm9,0.0,tmp8);

            tmp8bis.clear();
            Permute(tmp8,shape(3,7,1,6,0,2,4,5),tmp8bis);

            //add regular peps on upper site
            tmp9.clear();
            Gemm(CblasNoTrans,CblasNoTrans,1.0,(*this)(row+1,col),tmp8bis,0.0,tmp9);

            tmp9bis.clear();
            Permute(tmp9,shape(2,3,4,0,1,5,6,7,8),tmp9bis);

            //act with operator on peps on upper site: hermitian adjoint
            peps_op.clear();
            Contract(1.0,ham.gR(i),shape(k,j),(*this)(row+1,col),shape(l,m,k,n,o),0.0,peps_op,shape(l,m,j,n,o));

            //add it on the intermediary
            tmp8.clear();
            Gemm(CblasNoTrans,CblasNoTrans,1.0,peps_op,tmp9bis,0.0,tmp8);

            tmp8bis.clear();
            Permute(tmp8,shape(1,3,7,0,2,4,5,6),tmp8bis);

            //finally top environment for closure
            Gemm(CblasNoTrans,CblasNoTrans,1.0,env.gt(row)[col],tmp8bis,0.0,tmp6);

            //close down diagonal Left-down right-up
            val += ham.gcoef_nn(i) * Dot(tmp6,LOi_d[i]);

         }

         //construct new left going operators

         //first add top to left unit
         tmp8.clear();
         Gemm(CblasTrans,CblasNoTrans,1.0,env.gt(row)[col],RO[col - 1],0.0,tmp8);

         tmp8bis.clear();
         Permute(tmp8,shape(1,2,4,5,6,7,3,0),tmp8bis);

         //add regular upper peps to left
         tmp9.clear();
         Gemm(CblasNoTrans,CblasNoTrans,1.0,tmp8bis,(*this)(row+1,col),0.0,tmp9);

         perm9.clear();
         Permute(tmp9,shape(1,3,4,5,7,8,2,0,6),perm9);

         //construct Left Up operator first
         for(int i = 0;i < delta;++i){

            //add operator to upper peps
            peps_op.clear();
            Contract(1.0,ham.gL(i),shape(j,k),(*this)(row+1,col),shape(l,m,k,n,o),0.0,peps_op,shape(l,m,j,n,o));

            tmp8.clear();
            Gemm(CblasNoTrans,CblasNoTrans,1.0,perm9,peps_op,0.0,tmp8);

            tmp8bis.clear();
            Permute(tmp8,shape(0,5,7,2,3,6,1,4),tmp8bis);

            //add lower regular peps on
            tmp9.clear();
            Gemm(CblasNoTrans,CblasNoTrans,1.0,tmp8bis,(*this)(row,col),0.0,tmp9);

            tmp9bis.clear();
            Permute(tmp9,shape(0,1,2,4,7,8,3,5,6),tmp9bis);

            //and another lower regular peps
            tmp8.clear();
            Gemm(CblasNoTrans,CblasNoTrans,1.0,tmp9bis,(*this)(row,col),0.0,tmp8);

            tmp8bis.clear();
            Permute(tmp8,shape(0,1,2,5,7,3,4,6),tmp8bis);

            //finally construct Lu_i by contracting with bottom environment
            LOi_u[i].clear();
            Gemm(CblasNoTrans,CblasNoTrans,1.0,tmp8bis,env.gb(row-1)[col],0.0,LOi_u[i]);

         }

         //add on regular upper peps to intermediate perm9
         tmp8.clear();
         Gemm(CblasNoTrans,CblasNoTrans,1.0,perm9,(*this)(row+1,col),0.0,tmp8);

         tmp8bis.clear();
         Permute(tmp8,shape(0,5,7,2,3,6,1,4),tmp8bis);

         //add lower regular peps on
         tmp9.clear();
         Gemm(CblasNoTrans,CblasNoTrans,1.0,tmp8bis,(*this)(row,col),0.0,tmp9);

         tmp9bis.clear();
         Permute(tmp9,shape(0,1,2,4,7,8,3,5,6),tmp9bis);

         //construct lower-left going operator
         for(int i = 0;i < delta;++i){

            //add operator to lower peps
            peps_op.clear();
            Contract(1.0,ham.gL(i),shape(j,k),(*this)(row,col),shape(l,m,k,n,o),0.0,peps_op,shape(l,m,j,n,o));

            //paste it on intermediate tmp9bis
            tmp8.clear();
            Gemm(CblasNoTrans,CblasNoTrans,1.0,tmp9bis,peps_op,0.0,tmp8);

            tmp8bis.clear();
            Permute(tmp8,shape(0,1,2,5,7,3,4,6),tmp8bis);

            //finally construct Ld_i by contracting with bottom environment
            LOi_d[i].clear();
            Gemm(CblasNoTrans,CblasNoTrans,1.0,tmp8bis,env.gb(row-1)[col],0.0,LOi_d[i]);

         }

         //finally left unit operator: add on another regular lower peps
         Gemm(CblasNoTrans,CblasNoTrans,1.0,tmp9bis,(*this)(row,col),0.0,tmp8);

         Permute(tmp8,shape(0,1,2,5,7,3,4,6),tmp8bis);

         //finally construct Left Unit by contracting with bottom environment
         Gemm(CblasNoTrans,CblasNoTrans,1.0,tmp8bis,env.gb(row-1)[col],0.0,RO[col]);

      }

      //close down incoming operators on last site + vertical gate:

      //attach peps Lx - 1 to bottom
      tmp5.clear();
      Gemm(CblasNoTrans,CblasTrans,1.0,(*this)(row,Lx-1),env.gb(row-1)[Lx-1],0.0,tmp5);

      tmp5bis.clear();
      Permute(tmp5,shape(2,4,0,3,1),tmp5bis);

      //close down Left Up for diagonal, Left Down for horizontal and evaluate vertical
      for(int i = 0;i < delta;++i){

         //add operator on lower peps: hermitian adjoint
         peps_op.clear();
         Contract(1.0,ham.gR(i),shape(k,j),(*this)(row,Lx-1),shape(l,m,k,n,o),0.0,peps_op,shape(l,m,j,n,o));

         M = (*this)(row,Lx-1).shape(0) * (*this)(row,Lx-1).shape(1);
         N = tmp5bis.shape(2) * tmp5bis.shape(3) * tmp5bis.shape(4);
         K = (*this)(row,Lx-1).shape(2) * (*this)(row,Lx-1).shape(3);

         tmp5.resize(shape(D,D,D,env.gb(row-1)[Lx-1].shape(0),D));
         blas::gemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, M, N, K, 1.0,peps_op.data(),K,tmp5bis.data(),N,0.0,tmp5.data(),N);

         //add regular upper peps
         M = (*this)(row+1,Lx-1).shape(0) * (*this)(row+1,Lx-1).shape(1) * (*this)(row+1,Lx-1).shape(2);
         N = env.gb(row-1)[Lx-1].shape(0) * D * D * D;
         K = (*this)(row+1,Lx-1).shape(3);

         tmp7.resize( shape(D,D,d,D,D,D,env.gb(row-1)[Lx-1].shape(0)) );
         blas::gemm(CblasRowMajor, CblasNoTrans, CblasTrans, M, N, K, 1.0,(*this)(row+1,Lx-1).data(),K,tmp5.data(),K,0.0,tmp7.data(),N);

         perm7.clear();
         Permute(tmp7,shape(2,4,0,1,3,5,6),perm7);

         // --- 1) add operator for Vertical contribution
         peps_op.clear();
         Contract(1.0,ham.gL(i),shape(k,j),(*this)(row+1,Lx-1),shape(l,m,k,n,o),0.0,peps_op,shape(l,m,j,n,o));

         M = (*this)(row+1,Lx-1).shape(0) * (*this)(row+1,Lx-1).shape(1);
         N = env.gb(row-1)[Lx-1].shape(0) * D * D * D * D;
         K = (*this)(row+1,Lx-1).shape(2) * (*this)(row+1,Lx-1).shape(3);

         tmp7.resize( shape(D,D,D,D,D,D,env.gb(row-1)[Lx-1].shape(0)) );
         blas::gemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, M, N, K, 1.0,peps_op.data(),K,perm7.data(),N,0.0,tmp7.data(),N);

         tmp7bis.clear();
         Permute(tmp7,shape(1,3,0,2,4,5,6),tmp7bis);

         //finally top environment
         M = env.gt(row)[Lx-1].shape(0);
         N = env.gb(row-1)[Lx-1].shape(0) * D * D * D * D;
         K = env.gt(row)[Lx-1].shape(1) * env.gt(row)[Lx-1].shape(2);

         tmp6.resize(env.gt(row)[Lx-1].shape(0),D,D,D,D,env.gb(row-1)[Lx-1].shape(0));
         blas::gemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, M, N, K, 1.0,env.gt(row)[Lx-1].data(),K,tmp7bis.data(),N,0.0,tmp6.data(),N);

         //vertical contribution
         val += ham.gcoef_n(i) * Dot(tmp6,RO[Lx-2]);

         // --- 2) add regular peps' for Left Up and Horizontal energy contributions
         M = (*this)(row+1,Lx-1).shape(0) * (*this)(row+1,Lx-1).shape(1);
         N = env.gb(row-1)[Lx-1].shape(0) * D * D * D * D;
         K = (*this)(row+1,Lx-1).shape(2) * (*this)(row+1,Lx-1).shape(3);

         tmp7.resize( shape(D,D,D,D,D,D,env.gb(row-1)[Lx-1].shape(0)) );
         blas::gemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, M, N, K, 1.0,(*this)(row+1,Lx-1).data(),K,perm7.data(),N,0.0,tmp7.data(),N);

         tmp7bis.clear();
         Permute(tmp7,shape(1,3,0,2,4,5,6),tmp7bis);

         //finally top environment
         M = env.gt(row)[Lx-1].shape(0);
         N = env.gb(row-1)[Lx-1].shape(0) * D * D * D * D;
         K = env.gt(row)[Lx-1].shape(1) * env.gt(row)[Lx-1].shape(2);

         tmp6.resize(env.gt(row)[Lx-1].shape(0),D,D,D,D,env.gb(row-1)[Lx-1].shape(0));
         blas::gemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, M, N, K, 1.0,env.gt(row)[Lx-1].data(),K,tmp7bis.data(),N,0.0,tmp6.data(),N);

         //horizontal energy contribution
         val += ham.gcoef_n(i) * Dot(tmp6,LOi_d[i]);

         //diagonal energy contribution: Lu \ Rd
         val += ham.gcoef_nn(i) * Dot(tmp6,LOi_u[i]);

      }

      //finally close down Left Down for diagonal Ld / Ru
      for(int i = 0;i < delta;++i){

         //add regular lower peps
         M = (*this)(row,Lx-1).shape(0) * (*this)(row,Lx-1).shape(1);
         N = tmp5bis.shape(2) * tmp5bis.shape(3) * tmp5bis.shape(4);
         K = (*this)(row,Lx-1).shape(2) * (*this)(row,Lx-1).shape(3);

         tmp5.resize(shape(D,D,D,env.gb(row-1)[Lx-1].shape(0),D));
         blas::gemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, M, N, K, 1.0,(*this)(row,Lx-1).data(),K,tmp5bis.data(),N,0.0,tmp5.data(),N);

         //add regular upper peps
         M = (*this)(row+1,Lx-1).shape(0) * (*this)(row+1,Lx-1).shape(1) * (*this)(row+1,Lx-1).shape(2);
         N = env.gb(row-1)[Lx-1].shape(0) * D * D * D;
         K = (*this)(row+1,Lx-1).shape(3);

         tmp7.resize( shape(D,D,d,D,D,D,env.gb(row-1)[Lx-1].shape(0)) );
         blas::gemm(CblasRowMajor, CblasNoTrans, CblasTrans, M, N, K, 1.0,(*this)(row+1,Lx-1).data(),K,tmp5.data(),K,0.0,tmp7.data(),N);

         Permute(tmp7,shape(2,4,0,1,3,5,6),tmp7bis);

         //add right operator to upper peps: hermitian adjoint
         peps_op.clear();
         Contract(1.0,ham.gR(i),shape(k,j),(*this)(row+1,Lx-1),shape(l,m,k,n,o),0.0,peps_op,shape(l,m,j,n,o));

         M = (*this)(row+1,Lx-1).shape(0) * (*this)(row+1,Lx-1).shape(1);
         N = env.gb(row-1)[Lx-1].shape(0) * D * D * D * D;
         K = (*this)(row+1,Lx-1).shape(2) * (*this)(row+1,Lx-1).shape(3);

         tmp7.resize( shape(D,D,D,D,D,D,env.gb(row-1)[Lx-1].shape(0)) );
         blas::gemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, M, N, K, 1.0,peps_op.data(),K,tmp7bis.data(),N,0.0,tmp7.data(),N);

         tmp7bis.clear();
         Permute(tmp7,shape(1,3,0,2,4,5,6),tmp7bis);

         //finally top environment
         M = env.gt(row)[Lx-1].shape(0);
         N = env.gb(row-1)[Lx-1].shape(0) * D * D * D * D;
         K = env.gt(row)[Lx-1].shape(1) * env.gt(row)[Lx-1].shape(2);

         tmp6.resize(env.gt(row)[Lx-1].shape(0),D,D,D,D,env.gb(row-1)[Lx-1].shape(0));
         blas::gemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, M, N, K, 1.0,env.gt(row)[Lx-1].data(),K,tmp7bis.data(),N,0.0,tmp6.data(),N);

         //diagonal energy contribution: Ld / Ru
         val += ham.gcoef_nn(i) * Dot(tmp6,LOi_d[i]);

      }

   }

   // -- (3) -- || top row = Ly-1: again similar to overlap calculation

   //first construct the right renormalized operators
   contractions::init_ro('t',*this,R);

   //construct the left going operators and the vertical energy contribution
   
   //first construct Left Up operator and vertcial energy
   for(int i = 0;i < delta;++i){

      peps_op.clear();
      Contract(1.0,ham.gL(i),shape(j,k),(*this)(Ly-1,0),shape(l,m,k,n,o),0.0,peps_op,shape(l,m,j,n,o));

      //connect upper regular and operator peps together
      Gemm(CblasTrans,CblasNoTrans,1.0,(*this)(Ly-1,0),peps_op,0.0,tmp4);

      tmp4bis.clear();
      Permute(tmp4,shape(1,3,2,0),tmp4bis);

      //paste on lower regular peps
      M = tmp4bis.shape(0) * tmp4bis.shape(1) * tmp4bis.shape(2);
      N = (*this)(Ly-2,0).shape(2) * (*this)(Ly-2,0).shape(3) * (*this)(Ly-2,0).shape(4);
      K = tmp4bis.shape(3);

      tmp6.resize( shape(D,D,D,d,D,D) );
      blas::gemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, M, N, K, 1.0,tmp4bis.data(),K,(*this)(Ly-2,0).data(),N,0.0,tmp6.data(),N);

      perm6.clear();
      Permute(tmp6,shape(0,1,4,5,2,3),perm6);

      //first do vertical energy contribution: operator on lower site
      peps_op.clear();
      Contract(1.0,ham.gR(i),shape(j,k),(*this)(Ly-2,0),shape(l,m,k,n,o),0.0,peps_op,shape(l,m,j,n,o));

      M = perm6.shape(0) * perm6.shape(1) * perm6.shape(2) * perm6.shape(3);
      N = (*this)(Ly-2,0).shape(3) * (*this)(Ly-2,0).shape(4);
      K = perm6.shape(4) * perm6.shape(5);

      tmp6.resize( shape(D,D,D,D,D,D) );
      blas::gemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, M, N, K, 1.0,perm6.data(),K,peps_op.data(),N,0.0,tmp6.data(),N);

      tmp6bis.clear();
      Permute(tmp6,shape(0,1,3,5,2,4),tmp6bis);

      //finally contract with bottom environment
      M = tmp6bis.shape(0) * tmp6bis.shape(1) * tmp6bis.shape(2) * tmp6bis.shape(3);
      N = env.gb(Ly-3)[0].shape(3);
      K = tmp6bis.shape(4) * tmp6bis.shape(5);

      Li_u[i].resize( shape(D,D,D,D,env.gb(Ly-3)[0].shape(3)) );
      blas::gemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, M, N, K, 1.0,tmp6bis.data(),K,env.gb(Ly-3)[0].data(),N,0.0,Li_u[i].data(),N);

      //vertical energy:
      val += ham.gcoef_n(i) * Dot(Li_u[i],R[0]);

      //now construct the left up going operator: add regular lower peps to intermediate perm6
      M = perm6.shape(0) * perm6.shape(1) * perm6.shape(2) * perm6.shape(3);
      N = (*this)(Ly-2,0).shape(3) * (*this)(Ly-2,0).shape(4);
      K = perm6.shape(4) * perm6.shape(5);

      tmp6.resize( shape(D,D,D,D,D,D) );
      blas::gemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, M, N, K, 1.0,perm6.data(),K,(*this)(Ly-2,0).data(),N,0.0,tmp6.data(),N);

      tmp6bis.clear();
      Permute(tmp6,shape(0,1,3,5,2,4),tmp6bis);

      //finally contract with bottom environment
      M = tmp6bis.shape(0) * tmp6bis.shape(1) * tmp6bis.shape(2) * tmp6bis.shape(3);
      N = env.gb(Ly-3)[0].shape(3);
      K = tmp6bis.shape(4) * tmp6bis.shape(5);

      blas::gemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, M, N, K, 1.0,tmp6bis.data(),K,env.gb(Ly-3)[0].data(),N,0.0,Li_u[i].data(),N);

   }

   //then construct lower left operator and unity
   Gemm(CblasTrans,CblasNoTrans,1.0,(*this)(Ly-1,0),(*this)(Ly-1,0),0.0,tmp4);

   tmp4bis.clear();
   Permute(tmp4,shape(1,3,2,0),tmp4bis);

   //paste on lower regular peps
   M = tmp4bis.shape(0) * tmp4bis.shape(1) * tmp4bis.shape(2);
   N = (*this)(Ly-2,0).shape(2) * (*this)(Ly-2,0).shape(3) * (*this)(Ly-2,0).shape(4);
   K = tmp4bis.shape(3);

   tmp6.resize( shape(D,D,D,d,D,D) );
   blas::gemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, M, N, K, 1.0,tmp4bis.data(),K,(*this)(Ly-2,0).data(),N,0.0,tmp6.data(),N);

   perm6.clear();
   Permute(tmp6,shape(0,1,4,5,2,3),perm6);

   //construct lower left operator
   for(int i = 0;i < delta;++i){

      //first do vertical energy contribution: operator on lower site
      peps_op.clear();
      Contract(1.0,ham.gL(i),shape(j,k),(*this)(Ly-2,0),shape(l,m,k,n,o),0.0,peps_op,shape(l,m,j,n,o));

      //attach to perm6
      M = perm6.shape(0) * perm6.shape(1) * perm6.shape(2) * perm6.shape(3);
      N = (*this)(Ly-2,0).shape(3) * (*this)(Ly-2,0).shape(4);
      K = perm6.shape(4) * perm6.shape(5);

      tmp6.resize( shape(D,D,D,D,D,D) );
      blas::gemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, M, N, K, 1.0,perm6.data(),K,peps_op.data(),N,0.0,tmp6.data(),N);

      tmp6bis.clear();
      Permute(tmp6,shape(0,1,3,5,2,4),tmp6bis);

      //finally contract with bottom environment to get Left Down
      M = tmp6bis.shape(0) * tmp6bis.shape(1) * tmp6bis.shape(2) * tmp6bis.shape(3);
      N = env.gb(Ly-3)[0].shape(3);
      K = tmp6bis.shape(4) * tmp6bis.shape(5);

      Li_d[i].resize( shape(D,D,D,D,env.gb(Ly-3)[0].shape(3)) );
      blas::gemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, M, N, K, 1.0,tmp6bis.data(),K,env.gb(Ly-3)[0].data(),N,0.0,Li_d[i].data(),N);

   }

   //and finally unity,  add regular lower peps to intermediate perm6
   M = perm6.shape(0) * perm6.shape(1) * perm6.shape(2) * perm6.shape(3);
   N = (*this)(Ly-2,0).shape(3) * (*this)(Ly-2,0).shape(4);
   K = perm6.shape(4) * perm6.shape(5);

   tmp6.resize( shape(D,D,D,D,D,D) );
   blas::gemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, M, N, K, 1.0,perm6.data(),K,(*this)(Ly-2,0).data(),N,0.0,tmp6.data(),N);

   tmp6bis.clear();
   Permute(tmp6,shape(0,1,3,5,2,4),tmp6bis);

   //finally contract with bottom environment
   M = tmp6bis.shape(0) * tmp6bis.shape(1) * tmp6bis.shape(2) * tmp6bis.shape(3);
   N = env.gb(Ly-3)[0].shape(3);
   K = tmp6bis.shape(4) * tmp6bis.shape(5);

   blas::gemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, M, N, K, 1.0,tmp6bis.data(),K,env.gb(Ly-3)[0].data(),N,0.0,R[0].data(),N);

   //middle of the chain:
   for(int col = 1;col < Lx-1;++col){

      //close down the left renormalized terms from the previous site
      tmp7.clear();
      Gemm(CblasNoTrans,CblasTrans,1.0,env.gb(Ly-3)[col],R[col],0.0,tmp7);

      tmp7bis.clear();
      Permute(tmp7,shape(2,6,0,1,3,4,5),tmp7bis);

      tmp8.clear();
      Gemm(CblasNoTrans,CblasNoTrans,1.0,(*this)(Ly-2,col),tmp7bis,0.0,tmp8);

      perm8.clear();
      Permute(tmp8,shape(2,4,7,0,1,3,5,6),perm8);

      //close down left up LU \ RD and horizontal, and make vertical contribution
      for(int i = 0;i < delta;++i){

         //add operator to lower site peps: hermitian adjoint because of change of direction
         peps_op.clear();
         Contract(1.0,ham.gR(i),shape(k,j),(*this)(Ly-2,col),shape(l,m,k,n,o),0.0,peps_op,shape(l,m,j,n,o));

         tmp7.clear();
         Gemm(CblasNoTrans,CblasNoTrans,1.0,peps_op,perm8,0.0,tmp7);

         tmp7bis.clear();
         Permute(tmp7,shape(3,6,0,2,4,1,5),tmp7bis);

         //add regular upper peps to intermediate
         tmp8.clear();
         Gemm(CblasNoTrans,CblasNoTrans,1.0,(*this)(Ly-1,col),tmp7bis,0.0,tmp8);

         tmp8bis.clear();
         Permute(tmp8,shape(1,2,6,7,0,3,4,5),tmp8bis);

         //add operator on upper peps for vertical contribution: hermitian adjoint
         peps_op.clear();
         Contract(1.0,ham.gL(i),shape(k,j),(*this)(Ly-1,col),shape(l,m,k,n,o),0.0,peps_op,shape(l,m,j,n,o));

         tmp5.clear();
         Gemm(CblasNoTrans,CblasNoTrans,1.0,peps_op,tmp8bis,0.0,tmp5);

         val += ham.gcoef_n(i) * Dot(tmp5,R[col-1]);

         //add regular upper peps to intermediate tmp8bis for closure of Left Up and Horizontal
         Gemm(CblasNoTrans,CblasNoTrans,1.0,(*this)(Ly-1,col),tmp8bis,0.0,tmp5);

         //close down with Left Down for horizontal contribution
         val += ham.gcoef_n(i) * Dot(tmp5,Li_d[i]);

         //and Left Up \ Right Down diagonal:
         val += ham.gcoef_nn(i) * Dot(tmp5,Li_u[i]);

      }

      //close down the Left-Down diagonal and upper horizontal:
      for(int i = 0;i < delta;++i){

         //add regular lower site peps to intermediate perm8
         Gemm(CblasNoTrans,CblasNoTrans,1.0,(*this)(Ly-2,col),perm8,0.0,tmp7);

         Permute(tmp7,shape(3,6,0,2,4,1,5),tmp7bis);

         //add regular upper peps to intermediate
         Gemm(CblasNoTrans,CblasNoTrans,1.0,(*this)(Ly-1,col),tmp7bis,0.0,tmp8);

         Permute(tmp8,shape(1,2,6,7,0,3,4,5),tmp8bis);

         //add operator on upper peps for diagonal energy evaluation: hermitian adjoint
         Contract(1.0,ham.gR(i),shape(k,j),(*this)(Ly-1,col),shape(l,m,k,n,o),0.0,peps_op,shape(l,m,j,n,o));

         Gemm(CblasNoTrans,CblasNoTrans,1.0,peps_op,tmp8bis,0.0,tmp5);

         //Left Down / Right Up diagonal:
         val += ham.gcoef_nn(i) * Dot(tmp5,Li_d[i]);

         //upper horizontal with Left-Up
         val += ham.gcoef_n(i) * Dot(tmp5,Li_u[i]);

      }

      //construct left renormalized operators for next site

      //add regular upper peps to Left Unity
      tmp8.clear();
      Gemm(CblasTrans,CblasNoTrans,1.0,R[col-1],(*this)(Ly-1,col),0.0,tmp8);

      perm8.clear();
      Permute(tmp8,shape(1,2,3,6,7,0,4,5),perm8);

      //first construct left up:
      for(int i = 0;i < delta;++i){

         //add operator to upper peps
         peps_op.clear();
         Contract(1.0,ham.gL(i),shape(j,k),(*this)(Ly-1,col),shape(l,m,k,n,o),0.0,peps_op,shape(l,m,j,n,o));

         tmp7.clear();
         Gemm(CblasNoTrans,CblasNoTrans,1.0,perm8,peps_op,0.0,tmp7);

         tmp7bis.clear();
         Permute(tmp7,shape(1,2,4,5,6,0,3),tmp7bis);

         //now add regular lower peps
         tmp8.clear();
         Gemm(CblasNoTrans,CblasNoTrans,1.0,tmp7bis,(*this)(Ly-2,col),0.0,tmp8);

         tmp8bis.clear();
         Permute(tmp8,shape(1,2,4,6,7,0,3,5),tmp8bis);

         //and another
         tmp7.clear();
         Gemm(CblasNoTrans,CblasNoTrans,1.0,tmp8bis,(*this)(Ly-2,col),0.0,tmp7);

         tmp7bis.clear();
         Permute(tmp7,shape(1,2,4,6,0,3,5),tmp7bis);

         //finally construct Left-Up operator by contracting with bottom environment
         Li_u[i].clear();
         Gemm(CblasNoTrans,CblasNoTrans,1.0,tmp7bis,env.gb(Ly-3)[col],0.0,Li_u[i]);

      }

      //then construct the Left-Down operator and unity

      //first construct left up:

      //add regular upper peps to intermediate perm8
      tmp7.clear();
      Gemm(CblasNoTrans,CblasNoTrans,1.0,perm8,(*this)(Ly-1,col),0.0,tmp7);

      tmp7bis.clear();
      Permute(tmp7,shape(1,2,4,5,6,0,3),tmp7bis);

      //now add regular lower peps
      tmp8.clear();
      Gemm(CblasNoTrans,CblasNoTrans,1.0,tmp7bis,(*this)(Ly-2,col),0.0,tmp8);

      tmp8bis.clear();
      Permute(tmp8,shape(1,2,4,6,7,0,3,5),tmp8bis);

      //construct Left Down operator
      for(int i = 0;i < delta;++i){

         //add operator to lower peps
         peps_op.clear();
         Contract(1.0,ham.gL(i),shape(j,k),(*this)(Ly-2,col),shape(l,m,k,n,o),0.0,peps_op,shape(l,m,j,n,o));

         //and add it to tmp8bis
         tmp7.clear();
         Gemm(CblasNoTrans,CblasNoTrans,1.0,tmp8bis,peps_op,0.0,tmp7);

         tmp7bis.clear();
         Permute(tmp7,shape(1,2,4,6,0,3,5),tmp7bis);

         //construct Left-Down operator by contracting with bottom environment
         Li_d[i].clear();
         Gemm(CblasNoTrans,CblasNoTrans,1.0,tmp7bis,env.gb(Ly-3)[col],0.0,Li_d[i]);

      }

      //and lastly make the unit operator

      //add lower regular peps to intemediate tmp8bis
      tmp7.clear();
      Gemm(CblasNoTrans,CblasNoTrans,1.0,tmp8bis,(*this)(Ly-2,col),0.0,tmp7);

      tmp7bis.clear();
      Permute(tmp7,shape(1,2,4,6,0,3,5),tmp7bis);

      //construct Left-Down operator by contracting with bottom environment
      Gemm(CblasNoTrans,CblasNoTrans,1.0,tmp7bis,env.gb(Ly-3)[col],0.0,R[col]);

   }

   //finally close down on last top sites

   //attach bottom environment to lower peps
   tmp5.clear();
   Gemm(CblasNoTrans,CblasTrans,1.0,(*this)(Ly-2,Lx-1),env.gb(Ly-3)[Lx-1],0.0,tmp5);

   perm5.clear();
   Permute(tmp5,shape(2,4,0,1,3),perm5);

   //first close down Left-Up diagonal, horizontal and construct vertical
   for(int i = 0;i < delta;++i){

      //add operator on lower peps and contract with tmp5: (notice that the hermitian adjoint of the operator is applied, since the order of operations is reversed)
      peps_op.clear();
      Contract(1.0,ham.gR(i),shape(k,j),(*this)(Ly-2,Lx-1),shape(l,m,k,n,o),0.0,peps_op,shape(l,m,j,n,o));
      
      M = (*this)(Ly-2,Lx-1).shape(0) * (*this)(Ly-2,Lx-1).shape(1);
      N = perm5.shape(2) * perm5.shape(3) * perm5.shape(4);
      K = (*this)(Ly-2,Lx-1).shape(2) * (*this)(Ly-2,Lx-1).shape(3);

      tmp5.resize(shape( D,D,D,D,env.gb(Ly-3)[Lx-1].shape(0) ) );
      blas::gemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, M, N, K, 1.0,peps_op.data(),K,perm5.data(),N,0.0,tmp5.data(),N);

      tmp5bis.clear();
      Permute(tmp5,shape(3,0,1,2,4),tmp5bis);

      //add regular upper peps to block
      M = (*this)(Ly-1,Lx-1).shape(0) * (*this)(Ly-1,Lx-1).shape(2);
      N = tmp5bis.shape(1) * tmp5bis.shape(2) * tmp5bis.shape(3) * tmp5bis.shape(4);
      K = (*this)(Ly-1,Lx-1).shape(3);

      tmp6.resize(shape( D,d,D,D,D,env.gb(Ly-3)[Lx-1].shape(0) ) );
      blas::gemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, M, N, K, 1.0,(*this)(Ly-1,Lx-1).data(),K,tmp5bis.data(),N,0.0,tmp6.data(),N);

      tmp6bis.clear();
      Permute(tmp6,shape(1,3,0,2,4,5),tmp6bis);

      //now for the vertical energy contribution, add operator on upper peps: again hermitian adjoint
      peps_op.clear();
      Contract(1.0,ham.gL(i),shape(k,j),(*this)(Ly-1,Lx-1),shape(l,m,k,n,o),0.0,peps_op,shape(l,m,j,n,o));

      M = (*this)(Ly-1,Lx-1).shape(0);
      N = tmp6bis.shape(2) * tmp6bis.shape(3) * tmp6bis.shape(4) * tmp6bis.shape(5);
      K = tmp6bis.shape(0) * tmp6bis.shape(1);

      tmp5.resize(shape( D,D,D,D,env.gb(Ly-3)[Lx-1].shape(0) ) );
      blas::gemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, M, N, K, 1.0,peps_op.data(),K,tmp6bis.data(),N,0.0,tmp5.data(),N);

      val += ham.gcoef_n(i) * Dot(tmp5,R[Lx-2]);

      //for horizontal and Left-Up closure, add regular upper peps to intermediate
      M = (*this)(Ly-1,Lx-1).shape(0);
      N = tmp6bis.shape(2) * tmp6bis.shape(3) * tmp6bis.shape(4) * tmp6bis.shape(5);
      K = tmp6bis.shape(0) * tmp6bis.shape(1);

      tmp5.resize(shape( D,D,D,D,env.gb(Ly-3)[Lx-1].shape(0) ) );
      blas::gemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, M, N, K, 1.0,(*this)(Ly-1,Lx-1).data(),K,tmp6bis.data(),N,0.0,tmp5.data(),N);

      //horizontal with Left-Down
      val += ham.gcoef_n(i) * Dot(tmp5,Li_d[i]);

      //diagonal with Left-Up : LU \ RD
      val += ham.gcoef_nn(i) * Dot(tmp5,Li_u[i]);

   }

   //now close diagonal with left down, and horizontal on top row:
   for(int i = 0;i < delta;++i){

      //add regular lower peps to intermediate perm5
      M = (*this)(Ly-2,Lx-1).shape(0) * (*this)(Ly-2,Lx-1).shape(1);
      N = perm5.shape(2) * perm5.shape(3) * perm5.shape(4);
      K = (*this)(Ly-2,Lx-1).shape(2) * (*this)(Ly-2,Lx-1).shape(3);

      tmp5.resize(shape( D,D,D,D,env.gb(Ly-3)[Lx-1].shape(0) ) );
      blas::gemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, M, N, K, 1.0,(*this)(Ly-2,Lx-1).data(),K,perm5.data(),N,0.0,tmp5.data(),N);

      tmp5bis.clear();
      Permute(tmp5,shape(3,0,1,2,4),tmp5bis);

      //add regular upper peps to block
      M = (*this)(Ly-1,Lx-1).shape(0) * (*this)(Ly-1,Lx-1).shape(2);
      N = tmp5bis.shape(1) * tmp5bis.shape(2) * tmp5bis.shape(3) * tmp5bis.shape(4);
      K = (*this)(Ly-1,Lx-1).shape(3);

      tmp6.resize(shape( D,d,D,D,D,env.gb(Ly-3)[Lx-1].shape(0) ) );
      blas::gemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, M, N, K, 1.0,(*this)(Ly-1,Lx-1).data(),K,tmp5bis.data(),N,0.0,tmp6.data(),N);

      tmp6bis.clear();
      Permute(tmp6,shape(1,3,0,2,4,5),tmp6bis);

      //add operator on upper peps: again hermitian adjoint
      peps_op.clear();
      Contract(1.0,ham.gR(i),shape(k,j),(*this)(Ly-1,Lx-1),shape(l,m,k,n,o),0.0,peps_op,shape(l,m,j,n,o));

      M = (*this)(Ly-1,Lx-1).shape(0);
      N = tmp6bis.shape(2) * tmp6bis.shape(3) * tmp6bis.shape(4) * tmp6bis.shape(5);
      K = tmp6bis.shape(0) * tmp6bis.shape(1);

      tmp5.resize(shape( D,D,D,D,env.gb(Ly-3)[Lx-1].shape(0) ) );
      blas::gemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, M, N, K, 1.0,peps_op.data(),K,tmp6bis.data(),N,0.0,tmp5.data(),N);

      //upper horizontal contribution
      val += ham.gcoef_n(i) * Dot(tmp5,Li_u[i]);

      //diagonal with Left-Down : LD / RU
      val += ham.gcoef_nn(i) * Dot(tmp5,Li_d[i]);

   }
*/
   return val;

}

/**
 * 'canonicalize' a single peps row
 * @param dir Left or Right canonicalization
 * @param norm if true: normalize, else not
 */
template<typename T>
void PEPS<T>::canonicalize(int row,const BTAS_SIDE &dir,bool norm){

   if(dir == Left){//QR

      TArray<T,2> R;
      TArray<T,5> tmp;

      for(int i = 0;i < global::Lx - 1;++i){

         R.clear();

         //do QR
         Geqrf((*this)(row,i),R);

         //paste to next matrix
         tmp.clear();

         Contract((T)1.0,R,shape(1),(*this)(row,i + 1),shape(0),(T)0.0,tmp);

         (*this)(row,i + 1) = std::move(tmp);

      }

      if(norm){

         T nrm = sqrt( Dotc( (*this)(row,Lx-1),(*this)(row,Lx-1) ) );
         Scal(1.0/nrm,(*this)(row,Lx-1));

      }

   }
   else{//LQ

      TArray<T,2> L;
      TArray<T,5> tmp;

      for(int i = global::Lx - 1;i > 0;--i){

         L.clear();

         //do QR
         Gelqf(L,(*this)(row,i));

         //paste to previous matrix
         tmp.clear();

         Contract((T)1.0,(*this)(row,i - 1),shape(4),L,shape(0),(T)0.0,tmp);

         (*this)(row,i - 1) = std::move(tmp);

      }

      if(norm){

         T nrm = sqrt(Dotc((*this)(row,0),(*this)(row,0)));
         Scal(1.0/nrm,(*this)(row,0));

      }

   }
   
}
 
//forward declarations for types to be used!
template PEPS<double>::PEPS();
template PEPS< complex<double> >::PEPS();

template PEPS<double>::PEPS(int);
template PEPS< complex<double> >::PEPS(int);

template PEPS<double>::PEPS(const PEPS<double> &);
template PEPS< complex<double> >::PEPS(const PEPS< complex<double> > &);

template PEPS<double>::~PEPS();
template PEPS< complex<double> >::~PEPS();

template TArray<double,5> &PEPS<double>::operator()(int r,int c);
template TArray<complex<double>,5> &PEPS< complex<double> >::operator()(int r,int c);

template const TArray<double,5> &PEPS<double>::operator()(int r,int c) const;
template const TArray<complex<double>,5> &PEPS< complex<double> >::operator()(int r,int c) const;

template int PEPS<double>::gD() const;
template int PEPS< complex<double> >::gD() const;

template void PEPS<double>::sD(int);
template void PEPS< complex<double> >::sD(int);

template void PEPS<double>::scal(double val);
template void PEPS< complex<double> >::scal(complex<double> val);

template void PEPS<double>::load(const char *filename);
template void PEPS< complex<double> >::load(const char *filename);

template void PEPS<double>::fill_Random();
template void PEPS< complex<double> >::fill_Random();

template void PEPS<double>::save(const char *filename);
template void PEPS< complex<double> >::save(const char *filename);

template void PEPS<double>::canonicalize(int row,const BTAS_SIDE &dir,bool norm);
template void PEPS< complex<double> >::canonicalize(int row,const BTAS_SIDE &dir,bool norm);
