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
 * rescale all the tensors, set largest value to num
 * @param num number to rescale to
 */
template<>
void PEPS<double>::rescale_tensors(double num){

   for(int row = 0;row < Ly;++row)
      for(int col = 0;col < Lx;++col)
         (*this)(row,col).rescale(num);

}

/**
 * @param row the row index...
 * @param num number to rescale to
 * rescale all the tensors, set largest value on row 'row' to one
 */
template<>
void PEPS<double>::rescale_tensors(int row,double num){

   for(int col = 0;col < Lx;++col)
      (*this)(row,col).rescale(num);

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

   enum {j,k,l,m,n,o};

   //peps contracted with a local operator
   DArray<5> peps_op;

   //energy will be stored here
   double val = 0.0;

   //loop over the columns
   for(int col = 0;col < Lx-1;++col){

      // (A) construct left going operators and evaluate the vertical term

      //first add top to left unity L:
      tmp7.clear();
      Contract(1.0,L,shape(0),env.gt(0)[col],shape(0),0.0,tmp7);

      //add peps
      tmp8.clear();
      Contract(1.0,tmp7,shape(0,4),(*this)(1,col),shape(0,1),0.0,tmp8);

      //add upper operator for Lu_i and vertical contribution
      for(int i = 0;i < delta;++i){

         //add left operator on top site
         peps_op.clear();
         Contract(1.0,ham.gL(i),shape(j,k),(*this)(1,col),shape(l,m,k,n,o),0.0,peps_op,shape(l,m,j,n,o));

         tmp7.clear();
         Contract(1.0,tmp8,shape(0,3,5),peps_op,shape(0,1,2),0.0,tmp7);

         //paste bottom regular peps on
         tmp8bis.clear();
         Contract(1.0,tmp7,shape(0,3),(*this)(0,col),shape(0,1),0.0,tmp8bis);

         //1) first make the vertical energy contribution
         peps_op.clear();
         Contract(1.0,ham.gR(i),shape(j,k),(*this)(0,col),shape(l,m,k,n,o),0.0,peps_op,shape(l,m,j,n,o));

         Li_u[i].clear();
         Contract(1.0,tmp8bis,shape(0,3,5,6),peps_op,shape(0,1,2,3),0.0,Li_u[i]);

         val += ham.gcoef_n(i) * Dot(Li_u[i],R[col]);

         //2) the construct the 'real' Lu_i
         Contract(1.0,tmp8bis,shape(0,3,5,6),(*this)(0,col),shape(0,1,2,3),0.0,Li_u[i]);

      }

      //add regular peps to tmp8
      tmp7.clear();
      Contract(1.0,tmp8,shape(0,3,5),(*this)(1,col),shape(0,1,2),0.0,tmp7);
      
      //paste bottom regular peps on
      tmp8.clear();
      Contract(1.0,tmp7,shape(0,3),(*this)(0,col),shape(0,1),0.0,tmp8);

      //then construct the lower operator Ld_i
      for(int i = 0;i < delta;++i){

         //construct lower left operator
         peps_op.clear();
         Contract(1.0,ham.gL(i),shape(j,k),(*this)(0,col),shape(l,m,k,n,o),0.0,peps_op,shape(l,m,j,n,o));

         Li_d[i].clear();
         Contract(1.0,tmp8,shape(0,3,5,6),peps_op,shape(0,1,2,3),0.0,Li_d[i]);

      }

      //finally make new unity by adding one more regular peps
      L.clear();
      Contract(1.0,tmp8,shape(0,3,5,6),(*this)(0,col),shape(0,1,2,3),0.0,L);

      // (B) then close down the left renormalized operators for lurd,ldru and horizontal terms

      //start with Left Up
      for(int i = 0;i < delta;++i){

         //first add top to left unity L:
         tmp7.clear();
         Contract(1.0,Li_u[i],shape(0),env.gt(0)[col+1],shape(0),0.0,tmp7);

         //add peps
         tmp8.clear();
         Contract(1.0,tmp7,shape(0,4),(*this)(1,col+1),shape(0,1),0.0,tmp8);

         tmp7.clear();
         Contract(1.0,tmp8,shape(0,3,5),(*this)(1,col+1),shape(0,1,2),0.0,tmp7);

         //paste bottom regular peps on
         tmp8.clear();
         Contract(1.0,tmp7,shape(0,3),(*this)(0,col+1),shape(0,1),0.0,tmp8);

         //and construct the left up-down operator
         peps_op.clear();
         Contract(1.0,ham.gR(i),shape(j,k),(*this)(0,col+1),shape(l,m,k,n,o),0.0,peps_op,shape(l,m,j,n,o));

         Li_u[i].clear();
         Contract(1.0,tmp8,shape(0,3,5,6),peps_op,shape(0,1,2,3),0.0,Li_u[i]);

         val += ham.gcoef_nn(i) * Dot(Li_u[i],R[col+1]);

      }

      //Left down - close down with a diagonal and horizontal term!
      for(int i = 0;i < delta;++i){

         //first add on top to Left renormalized operator
         tmp7.clear();
         Contract(1.0,Li_d[i],shape(0),env.gt(0)[col+1],shape(0),0.0,tmp7);

         //add peps
         tmp8.clear();
         Contract(1.0,tmp7,shape(0,4),(*this)(1,col+1),shape(0,1),0.0,tmp8);

         //1) do the diagonal link first
         peps_op.clear();
         Contract(1.0,ham.gR(i),shape(j,k),(*this)(1,col+1),shape(l,m,k,n,o),0.0,peps_op,shape(l,m,j,n,o));

         tmp7.clear();
         Contract(1.0,tmp8,shape(0,3,5),peps_op,shape(0,1,2),0.0,tmp7);

         //paste bottom regular peps on
         tmp8bis.clear();
         Contract(1.0,tmp7,shape(0,3),(*this)(0,col+1),shape(0,1),0.0,tmp8bis);

         //and again
         tmp5.clear();
         Contract(1.0,tmp8bis,shape(0,3,5,6),(*this)(0,col+1),shape(0,1,2,3),0.0,tmp5);

         //energy:
         val += ham.gcoef_nn(i) * Dot(tmp5,R[col+1]);

         //2) then do the horizontal gate:

         //add regular peps
         tmp7.clear();
         Contract(1.0,tmp8,shape(0,3,5),(*this)(1,col+1),shape(0,1,2),0.0,tmp7);

         //paste bottom regular peps on
         tmp8.clear();
         Contract(1.0,tmp7,shape(0,3),(*this)(0,col+1),shape(0,1),0.0,tmp8);

         //construct the left down-down operator (horizontal gate)
         peps_op.clear();
         Contract(1.0,ham.gR(i),shape(j,k),(*this)(0,col+1),shape(l,m,k,n,o),0.0,peps_op,shape(l,m,j,n,o));

         tmp5.clear();
         Contract(1.0,tmp8,shape(0,3,5,6),peps_op,shape(0,1,2,3),0.0,tmp5);

         val += ham.gcoef_n(i) * Dot(tmp5,R[col+1]);

      }

   }

   //finally the last vertical term:
   tmp7.clear();
   Contract(1.0,L,shape(0),env.gt(0)[Lx-1],shape(0),0.0,tmp7);

   //add peps
   tmp8.clear();
   Contract(1.0,tmp7,shape(0,4),(*this)(1,Lx-1),shape(0,1),0.0,tmp8);

   //add upper operator
   for(int i = 0;i < delta;++i){

      //add left operator on top site
      peps_op.clear();
      Contract(1.0,ham.gL(i),shape(j,k),(*this)(1,Lx-1),shape(l,m,k,n,o),0.0,peps_op,shape(l,m,j,n,o));

      tmp7.clear();
      Contract(1.0,tmp8,shape(0,3,5),peps_op,shape(0,1,2),0.0,tmp7);

      //paste bottom regular peps on
      tmp8bis.clear();
      Contract(1.0,tmp7,shape(0,3),(*this)(0,Lx-1),shape(0,1),0.0,tmp8bis);

      //1) first make the vertical energy contribution
      peps_op.clear();
      Contract(1.0,ham.gR(i),shape(j,k),(*this)(0,Lx-1),shape(l,m,k,n,o),0.0,peps_op,shape(l,m,j,n,o));

      Li_u[i].clear();
      Contract(1.0,tmp8bis,shape(0,3,5,6),peps_op,shape(0,1,2,3),0.0,Li_u[i]);

      val += ham.gcoef_n(i) * Dot(Li_u[i],R[Lx-1]);

   }

   // -- (2) -- now move from bottom to top calculating everything like an MPO/MPS expectation value

   //Right renormalized operators
   vector< DArray<6> > RO(Lx);

   //left unity
   DArray<6> LO;

   //array of delta left upper and lower renormalized operators needed
   std::vector< DArray<6> > LOi_u(delta);
   std::vector< DArray<6> > LOi_d(delta);

   for(int row = 1;row < Ly - 2;++row){

      //first create right renormalized operator
      contractions::init_ro(row,*this,RO);

      LO.resize(shape(1,1,1,1,1,1));
      LO = 1.0;

      // --- move from left to right to get the expecation value of the interactions ---
      for(int col = 0;col < Lx - 1;++col){
         
         // (A) construct left going operators and evaluate vertical term

         //first add top to left unit
         tmp8.clear();
         Contract(1.0,LO,shape(0),env.gt(row)[col],shape(0),0.0,tmp8);

         //add regular upper peps to left
         tmp9.clear();
         Contract(1.0,tmp8,shape(0,5),(*this)(row+1,col),shape(0,1),0.0,tmp9);

         //construct Left Up operator first
         for(int i = 0;i < delta;++i){

            //add operator to upper peps
            peps_op.clear();
            Contract(1.0,ham.gL(i),shape(j,k),(*this)(row+1,col),shape(l,m,k,n,o),0.0,peps_op,shape(l,m,j,n,o));

            tmp8.clear();
            Contract(1.0,tmp9,shape(0,4,6),peps_op,shape(0,1,2),0.0,tmp8);

            //add lower regular peps on
            tmp9bis.clear();
            Contract(1.0,tmp8,shape(0,4),(*this)(row,col),shape(0,1),0.0,tmp9bis);

            //add operator to lower peps for vertical gate
            peps_op.clear();
            Contract(1.0,ham.gR(i),shape(j,k),(*this)(row,col),shape(l,m,k,n,o),0.0,peps_op,shape(l,m,j,n,o));

            //and add on tmp9bis
            tmp8.clear();
            Contract(1.0,tmp9bis,shape(0,4,6),peps_op,shape(0,1,2),0.0,tmp8);

            LOi_u[i].clear();
            Contract(1.0,tmp8,shape(0,4,6),env.gb(row-1)[col],shape(0,1,2),0.0,LOi_u[i]);

            //add vertical energy contribution
            val += ham.gcoef_n(i) * Dot(LOi_u[i],RO[col]);

            //then construct Lu_i by adding regular peps to tmp9bis
            Contract(1.0,tmp9bis,shape(0,4,6),(*this)(row,col),shape(0,1,2),0.0,tmp8);

            //and add bottom environment
            Contract(1.0,tmp8,shape(0,4,6),env.gb(row-1)[col],shape(0,1,2),0.0,LOi_u[i]);

         }

         //add on regular upper peps to intermediate perm9
         tmp8.clear();
         Contract(1.0,tmp9,shape(0,4,6),(*this)(row+1,col),shape(0,1,2),0.0,tmp8);

         //add lower regular peps on
         tmp9.clear();
         Contract(1.0,tmp8,shape(0,4),(*this)(row,col),shape(0,1),0.0,tmp9);

         //construct lower-left going operator
         for(int i = 0;i < delta;++i){

            //add operator to lower peps
            peps_op.clear();
            Contract(1.0,ham.gL(i),shape(j,k),(*this)(row,col),shape(l,m,k,n,o),0.0,peps_op,shape(l,m,j,n,o));

            //and add on tmp9
            tmp8.clear();
            Contract(1.0,tmp9,shape(0,4,6),peps_op,shape(0,1,2),0.0,tmp8);

            //and add bottom environment for LOi_d[i] construction
            LOi_d[i].clear();
            Contract(1.0,tmp8,shape(0,4,6),env.gb(row-1)[col],shape(0,1,2),0.0,LOi_d[i]);

         }

         //finally left unit operator: add on another regular lower peps
         tmp8.clear();
         Contract(1.0,tmp9,shape(0,4,6),(*this)(row,col),shape(0,1,2),0.0,tmp8);

         LO.clear();
         Contract(1.0,tmp8,shape(0,4,6),env.gb(row-1)[col],shape(0,1,2),0.0,LO);

         // (B) close down the left up and down operators with two diagonal and one horizontal contribution
         
         //start with Left Up
         for(int i = 0;i < delta;++i){

            //first add top to left-renormalized operator
            tmp8.clear();
            Contract(1.0,LOi_u[i],shape(0),env.gt(row)[col+1],shape(0),0.0,tmp8);

            //add regular upper peps to left
            tmp9.clear();
            Contract(1.0,tmp8,shape(0,5),(*this)(row+1,col+1),shape(0,1),0.0,tmp9);

            //and another
            tmp8.clear();
            Contract(1.0,tmp9,shape(0,4,6),(*this)(row+1,col+1),shape(0,1,2),0.0,tmp8);

            //add lower regular peps on
            tmp9.clear();
            Contract(1.0,tmp8,shape(0,4),(*this)(row,col+1),shape(0,1),0.0,tmp9);

            //add operator to lower peps for vertical gate
            peps_op.clear();
            Contract(1.0,ham.gR(i),shape(j,k),(*this)(row,col+1),shape(l,m,k,n,o),0.0,peps_op,shape(l,m,j,n,o));

            //and add on tmp9
            tmp8.clear();
            Contract(1.0,tmp9,shape(0,4,6),peps_op,shape(0,1,2),0.0,tmp8);

            //add bottom environment
            tmp6.clear();
            Contract(1.0,tmp8,shape(0,4,6),env.gb(row-1)[col+1],shape(0,1,2),0.0,tmp6);

            //diagonal-lurd energy contribution
            val += ham.gcoef_nn(i) * Dot(tmp6,RO[col+1]);

         }

         //Left down - close down with a diagonal and horizontal term!
         for(int i = 0;i < delta;++i){

            //first add top to left-renormalized operator
            tmp8.clear();
            Contract(1.0,LOi_d[i],shape(0),env.gt(row)[col+1],shape(0),0.0,tmp8);

            //add regular upper peps to left
            tmp9.clear();
            Contract(1.0,tmp8,shape(0,5),(*this)(row+1,col+1),shape(0,1),0.0,tmp9);

            //1) do the diagonal link first
            peps_op.clear();
            Contract(1.0,ham.gR(i),shape(j,k),(*this)(row+1,col+1),shape(l,m,k,n,o),0.0,peps_op,shape(l,m,j,n,o));

            tmp8.clear();
            Contract(1.0,tmp9,shape(0,4,6),peps_op,shape(0,1,2),0.0,tmp8);

            //add lower regular peps on
            tmp9bis.clear();
            Contract(1.0,tmp8,shape(0,4),(*this)(row,col+1),shape(0,1),0.0,tmp9bis);

            //and another
            tmp8.clear();
            Contract(1.0,tmp9bis,shape(0,4,6),(*this)(row,col+1),shape(0,1,2),0.0,tmp8);

            //and add bottom environment
            tmp6.clear();
            Contract(1.0,tmp8,shape(0,4,6),env.gb(row-1)[col+1],shape(0,1,2),0.0,tmp6);

            //diagnoal-ldru energy:
            val += ham.gcoef_nn(i) * Dot(tmp6,RO[col+1]);

            //2) then do the horizontal gate:

            //add regular top peps to intermediate tmp9
            tmp8.clear();
            Contract(1.0,tmp9,shape(0,4,6),(*this)(row+1,col+1),shape(0,1,2),0.0,tmp8);

            //add lower regular peps on
            tmp9bis.clear();
            Contract(1.0,tmp8,shape(0,4),(*this)(row,col+1),shape(0,1),0.0,tmp9bis);

            //construct the left down-down operator (horizontal gate)
            peps_op.clear();
            Contract(1.0,ham.gR(i),shape(j,k),(*this)(row,col+1),shape(l,m,k,n,o),0.0,peps_op,shape(l,m,j,n,o));

            //add interaction term
            tmp8.clear();
            Contract(1.0,tmp9bis,shape(0,4,6),peps_op,shape(0,1,2),0.0,tmp8);

            //and add bottom environment
            tmp6.clear();
            Contract(1.0,tmp8,shape(0,4,6),env.gb(row-1)[col+1],shape(0,1,2),0.0,tmp6);

            //horizontal energy contribution
            val += ham.gcoef_n(i) * Dot(tmp6,RO[col+1]);

         }

      }

      //last vertical gate

      //first add top to left unit
      tmp8.clear();
      Contract(1.0,LO,shape(0),env.gt(row)[Lx-1],shape(0),0.0,tmp8);

      //add regular upper peps to left
      tmp9.clear();
      Contract(1.0,tmp8,shape(0,5),(*this)(row+1,Lx-1),shape(0,1),0.0,tmp9);

      //construct Left Up operator first
      for(int i = 0;i < delta;++i){

         //add operator to upper peps
         peps_op.clear();
         Contract(1.0,ham.gL(i),shape(j,k),(*this)(row+1,Lx-1),shape(l,m,k,n,o),0.0,peps_op,shape(l,m,j,n,o));

         tmp8.clear();
         Contract(1.0,tmp9,shape(0,4,6),peps_op,shape(0,1,2),0.0,tmp8);

         //add lower regular peps on
         tmp9bis.clear();
         Contract(1.0,tmp8,shape(0,4),(*this)(row,Lx-1),shape(0,1),0.0,tmp9bis);

         //add operator to lower peps for vertical gate
         peps_op.clear();
         Contract(1.0,ham.gR(i),shape(j,k),(*this)(row,Lx-1),shape(l,m,k,n,o),0.0,peps_op,shape(l,m,j,n,o));

         //and add on tmp9bis
         tmp8.clear();
         Contract(1.0,tmp9bis,shape(0,4,6),peps_op,shape(0,1,2),0.0,tmp8);

         LOi_u[i].clear();
         Contract(1.0,tmp8,shape(0,4,6),env.gb(row-1)[Lx-1],shape(0,1,2),0.0,LOi_u[i]);

         //add vertical energy contribution
         val += ham.gcoef_n(i) * Dot(LOi_u[i],RO[Lx-1]);

      }

   }

   // -- (3) -- || top row = Ly-1: again similar to overlap calculation

   //first construct the right renormalized operators
   contractions::init_ro('t',*this,R);

   L.resize(shape(1,1,1,1,1));
   L = 1.0;

   //construct the left going operators and the vertical energy contribution

   //middle of the chain:
   for(int col = 0;col < Lx-1;++col){

      // (A) construct left renormalized operators and vertical energy contribution
      
      //add top peps to left
      tmp8.clear();
      Contract(1.0,L,shape(0),(*this)(Ly-1,col),shape(0),0.0,tmp8);

      //construct Left Up operator first
      for(int i = 0;i < delta;++i){

         //add operator to upper peps
         peps_op.clear();
         Contract(1.0,ham.gL(i),shape(j,k),(*this)(Ly-1,col),shape(l,m,k,n,o),0.0,peps_op,shape(l,m,j,n,o));

         //add to tmp8
         tmp7.clear();
         Contract(1.0,tmp8,shape(0,4,5),peps_op,shape(0,1,2),0.0,tmp7);

         //add regular peps
         tmp8bis.clear();
         Contract(1.0,tmp7,shape(0,3),(*this)(Ly-2,col),shape(0,1),0.0,tmp8bis);

         //right operator to lower peps
         peps_op.clear();
         Contract(1.0,ham.gR(i),shape(j,k),(*this)(Ly-2,col),shape(l,m,k,n,o),0.0,peps_op,shape(l,m,j,n,o));

         //add to tmp8bis
         tmp7.clear();
         Contract(1.0,tmp8bis,shape(0,3,5),peps_op,shape(0,1,2),0.0,tmp7);

         //finally contract with bottom environment
         tmp5.clear();
         Contract(1.0,tmp7,shape(0,3,5),env.gb(Ly-3)[col],shape(0,1,2),0.0,tmp5);

         //vertical energy
         val += ham.gcoef_n(i) * Dot(tmp5,R[col]);

         //construct Lu

         //add regular peps to tmp8bis
         tmp7.clear();
         Contract(1.0,tmp8bis,shape(0,3,5),(*this)(Ly-2,col),shape(0,1,2),0.0,tmp7);

         //finally contract with bottom environment
         Li_u[i].clear();
         Contract(1.0,tmp7,shape(0,3,5),env.gb(Ly-3)[col],shape(0,1,2),0.0,Li_u[i]);

      }

      //add regular upper peps to tmp8
      tmp7.clear();
      Contract(1.0,tmp8,shape(0,4,5),(*this)(Ly-1,col),shape(0,1,2),0.0,tmp7);

      //and add regular lower peps
      tmp8.clear();
      Contract(1.0,tmp7,shape(0,3),(*this)(Ly-2,col),shape(0,1),0.0,tmp8);

      //then construct Left Down operator
      for(int i = 0;i < delta;++i){

         //left operator to lower peps
         peps_op.clear();
         Contract(1.0,ham.gL(i),shape(j,k),(*this)(Ly-2,col),shape(l,m,k,n,o),0.0,peps_op,shape(l,m,j,n,o));

         //add to tmp8
         tmp7.clear();
         Contract(1.0,tmp8,shape(0,3,5),peps_op,shape(0,1,2),0.0,tmp7);

         //finally contract with bottom environment
         Li_d[i].clear();
         Contract(1.0,tmp7,shape(0,3,5),env.gb(Ly-3)[col],shape(0,1,2),0.0,Li_d[i]);

      }

      //finally make left unity

      //add lower regular peps to tmp8
      tmp7.clear();
      Contract(1.0,tmp8,shape(0,3,5),(*this)(Ly-2,col),shape(0,1,2),0.0,tmp7);

      //finally contract with bottom environment
      L.clear();
      Contract(1.0,tmp7,shape(0,3,5),env.gb(Ly-3)[col],shape(0,1,2),0.0,L);

      // (B) close down the left operators
      
      //start with left-up
      
      //add top peps to left
      for(int i = 0;i < delta;++i){

         // -- top row horizontal -- 
         tmp8.clear();
         Contract(1.0,Li_u[i],shape(0),(*this)(Ly-1,col+1),shape(0),0.0,tmp8);

         //add operator to upper peps
         peps_op.clear();
         Contract(1.0,ham.gR(i),shape(j,k),(*this)(Ly-1,col+1),shape(l,m,k,n,o),0.0,peps_op,shape(l,m,j,n,o));

         //add to tmp8
         tmp7.clear();
         Contract(1.0,tmp8,shape(0,4,5),peps_op,shape(0,1,2),0.0,tmp7);

         //add regular peps
         tmp8bis.clear();
         Contract(1.0,tmp7,shape(0,3),(*this)(Ly-2,col+1),shape(0,1),0.0,tmp8bis);

         //add another regular lower peps to tmp8bis
         tmp7.clear();
         Contract(1.0,tmp8bis,shape(0,3,5),(*this)(Ly-2,col+1),shape(0,1,2),0.0,tmp7);

         //finally contract with bottom environment
         tmp5.clear();
         Contract(1.0,tmp7,shape(0,3,5),env.gb(Ly-3)[col+1],shape(0,1,2),0.0,tmp5);

         //top row horizontal energy contribution
         val += ham.gcoef_n(i) * Dot(tmp5,R[col+1]);

         //-- lurd diagonal --

         //add regular upper peps to tmp8
         tmp7.clear();
         Contract(1.0,tmp8,shape(0,4,5),(*this)(Ly-1,col+1),shape(0,1,2),0.0,tmp7);

         //add regular lower peps
         tmp8.clear();
         Contract(1.0,tmp7,shape(0,3),(*this)(Ly-2,col+1),shape(0,1),0.0,tmp8);

         //add operator to lower peps
         peps_op.clear();
         Contract(1.0,ham.gR(i),shape(j,k),(*this)(Ly-2,col+1),shape(l,m,k,n,o),0.0,peps_op,shape(l,m,j,n,o));

         //add to tmp8bis
         tmp7.clear();
         Contract(1.0,tmp8,shape(0,3,5),peps_op,shape(0,1,2),0.0,tmp7);

         //finally contract with bottom environment
         tmp5.clear();
         Contract(1.0,tmp7,shape(0,3,5),env.gb(Ly-3)[col+1],shape(0,1,2),0.0,tmp5);

         //lurd diagonal energy contribution
         val += ham.gcoef_nn(i) * Dot(tmp5,R[col+1]);

      }
 
      //the nleft down
      
      //add top peps to left
      for(int i = 0;i < delta;++i){

         //-- ldru diagonal --
         tmp8.clear();
         Contract(1.0,Li_d[i],shape(0),(*this)(Ly-1,col+1),shape(0),0.0,tmp8);

         //add operator to upper peps
         peps_op.clear();
         Contract(1.0,ham.gR(i),shape(j,k),(*this)(Ly-1,col+1),shape(l,m,k,n,o),0.0,peps_op,shape(l,m,j,n,o));

         //add to tmp8
         tmp7.clear();
         Contract(1.0,tmp8,shape(0,4,5),peps_op,shape(0,1,2),0.0,tmp7);

         //add regular peps
         tmp8bis.clear();
         Contract(1.0,tmp7,shape(0,3),(*this)(Ly-2,col+1),shape(0,1),0.0,tmp8bis);

         //add another regular lower peps to tmp8bis
         tmp7.clear();
         Contract(1.0,tmp8bis,shape(0,3,5),(*this)(Ly-2,col+1),shape(0,1,2),0.0,tmp7);

         //finally contract with bottom environment
         tmp5.clear();
         Contract(1.0,tmp7,shape(0,3,5),env.gb(Ly-3)[col+1],shape(0,1,2),0.0,tmp5);

         //ldru-diagonal energy contribution
         val += ham.gcoef_nn(i) * Dot(tmp5,R[col+1]);

         //bottom row horizontal

         //add regular upper peps to tmp8
         tmp7.clear();
         Contract(1.0,tmp8,shape(0,4,5),(*this)(Ly-1,col+1),shape(0,1,2),0.0,tmp7);

         //add regular lower peps
         tmp8.clear();
         Contract(1.0,tmp7,shape(0,3),(*this)(Ly-2,col+1),shape(0,1),0.0,tmp8);

         //add operator to lower peps
         peps_op.clear();
         Contract(1.0,ham.gR(i),shape(j,k),(*this)(Ly-2,col+1),shape(l,m,k,n,o),0.0,peps_op,shape(l,m,j,n,o));

         //add to tmp8bis
         tmp7.clear();
         Contract(1.0,tmp8,shape(0,3,5),peps_op,shape(0,1,2),0.0,tmp7);

         //finally contract with bottom environment
         tmp5.clear();
         Contract(1.0,tmp7,shape(0,3,5),env.gb(Ly-3)[col+1],shape(0,1,2),0.0,tmp5);

         //bottom row horizontal energy contribution
         val += ham.gcoef_n(i) * Dot(tmp5,R[col+1]);

      }

   }

   //last row vertical update

   //add top peps to left
   tmp8.clear();
   Contract(1.0,L,shape(0),(*this)(Ly-1,Lx-1),shape(0),0.0,tmp8);

   //construct Left Up operator first
   for(int i = 0;i < delta;++i){

      //add operator to upper peps
      peps_op.clear();
      Contract(1.0,ham.gL(i),shape(j,k),(*this)(Ly-1,Lx-1),shape(l,m,k,n,o),0.0,peps_op,shape(l,m,j,n,o));

      //add to tmp8
      tmp7.clear();
      Contract(1.0,tmp8,shape(0,4,5),peps_op,shape(0,1,2),0.0,tmp7);

      //add regular peps
      tmp8bis.clear();
      Contract(1.0,tmp7,shape(0,3),(*this)(Ly-2,Lx-1),shape(0,1),0.0,tmp8bis);

      //right operator to lower peps
      peps_op.clear();
      Contract(1.0,ham.gR(i),shape(j,k),(*this)(Ly-2,Lx-1),shape(l,m,k,n,o),0.0,peps_op,shape(l,m,j,n,o));

      //add to tmp8bis
      tmp7.clear();
      Contract(1.0,tmp8bis,shape(0,3,5),peps_op,shape(0,1,2),0.0,tmp7);

      //finally contract with bottom environment
      tmp5.clear();
      Contract(1.0,tmp7,shape(0,3,5),env.gb(Ly-3)[Lx-1],shape(0,1,2),0.0,tmp5);

      //vertical energy
      val += ham.gcoef_n(i) * Dot(tmp5,R[Lx-1]);

   }

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
