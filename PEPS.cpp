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
double PEPS<double>::dot(const PEPS<double> &peps_i,bool init) const {

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
                        fout << a_ << "\t" << b_ << "\t" << c_ << "\t" << d_ << "\t" << e_ << "\t" << (*this)(row,col)(a_,b_,c_,d_,e_) << endl;

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
   vector< DArray<5> > R(Lx - 1);

   contractions::init_ro('b',*this,R); 

   //left going operators: Li
   std::vector< DArray<5> > Li_u( delta ); //upper site with extra operator
   std::vector< DArray<5> > Li_d( delta ); //lower site with extra operator


   //some storage stuff:
   DArray<5> tmp5;
   DArray<5> tmp5bis;

   DArray<6> tmp6;
   DArray<6> tmp6bis;

   DArray<7> tmp7;
   DArray<7> tmp7bis;

   DArray<8> tmp8;
   DArray<8> tmp8bis;
   DArray<8> perm8;

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
      val += ham.gcoef(i) * Dot(Li_d[i],R[0]);

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

   blas::gemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, M, N, K, 1.0,tmp5bis.data(),K,env.gb(0)[0].data(),N,0.0,R[0].data(),N);

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

         val += ham.gcoef(i) * global::J2 * Dot(Li_u[i],R[col]);

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

         val += ham.gcoef(i) * global::J2 * Dot(Li_d[i],R[col]);

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

         //construct the left down-down operator (horizonatal gate)
         peps_op.clear();
         Contract(1.0,ham.gR(i),shape(j,k),(*this)(0,col),shape(l,m,k,n,o),0.0,peps_op,shape(l,m,j,n,o));

         Gemm(CblasNoTrans,CblasNoTrans,1.0,tmp8bis,peps_op,0.0,Li_d[i]);

         val += ham.gcoef(i) * Dot(Li_d[i],R[col]);

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

         val += ham.gcoef(i) * Dot(Li_u[i],R[col]);

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

      val += ham.gcoef(i) * global::J2 * blas::dot(tmp8bis.size(),tmp8bis.data(),1,peps_op.data(),1);

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
      val += ham.gcoef(i) * global::J2 * blas::dot(tmp7bis.size(),tmp7bis.data(),1,env.gb(0)[Lx-1].data(),1);

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

      val += ham.gcoef(i) * blas::dot(tmp8bis.size(),tmp8bis.data(),1,peps_op.data(),1);

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

      val += blas::dot(tmp8bis.size(),tmp8bis.data(),1,peps_op.data(),1);

   }

   // -- (2) -- now move from bottom to top calculating everything like an MPO/MPS expectation value

   //Right renormalized operators
   vector< DArray<6> > RO(Lx - 1);

   //array of delta left renormalized operators needed
   std::vector< DArray<6> > LOi(delta);

   //for(int row = 1;row < Ly - 2;++row){
   int row = 1;

   //first create right renormalized operator
   contractions::init_ro(row,*this,RO);

   // --- now move from left to right to get the expecation value of the interactions ---
   
   //start with the first site, evaluate the vertical term and the left going terms

   //paste regular peps on top environment
   tmp5.clear();
   Gemm(CblasTrans,CblasNoTrans,1.0,env.gt(row)[0],(*this)(row,0),0.0,tmp5);

   tmp5bis.clear();
   Permute(tmp5,shape(1,3,4,0,2),tmp5bis);

   //construct  Left-Up operator and vertical energy contribution
   for(int i = 0;i < delta;++i){

      peps_op.clear();
      Contract(1.0,ham.gL(i),shape(j,k),(*this)(row,0),shape(l,m,k,n,o),0.0,peps_op,shape(l,m,j,n,o));

      M = 
         N
         K

      // add peps_i's to intermediate
      tmp8.clear();
      Contract(1.0,tmp7,shape(4,1),peps_op,shape(0,2),0.0,tmp8);

      tmp8bis.clear();
      Contract(1.0,tmp8,shape(3,6),env.gb(row-1)[0],shape(1,2),0.0,tmp8bis);

      //move to a DArray<3> object: order (top-env,(*this)-row,bottom-env)
      LOi[i] = tmp8bis.reshape_clear(shape(env.gt(row)[0].shape(3),(*this)(row,0).shape(4),(*this)(row,0).shape(4),env.gb(row-1)[0].shape(3)));

   }

   // 4) 1 -- finally construct left renormalized operator with unity
   Contract(1.0,tmp7,shape(1,4),(*this)(row,0),shape(1,2),0.0,tmp8);
   Contract(1.0,tmp8,shape(3,6),env.gb(row-1)[0],shape(1,2),0.0,tmp8bis);

   //move to a DArray<3> object: order (top-env,(*this)-row,bottom-env)
   RO[0] = tmp8bis.reshape_clear(shape(env.gt(row)[0].shape(3),(*this)(row,0).shape(4),(*this)(row,0).shape(4),env.gb(row-1)[0].shape(3)));

   // --- now for the middle sites, close down the operators on the left and construct new ones --- 
   for(int col = 1;col < Lx - 1;++col){

      //add top and one peps to the right side
      tmp6.clear();
      Contract(1.0,env.gt(row)[col],shape(3),RO[col],shape(0),0.0,tmp6);

      tmp7.clear();
      Contract(1.0,tmp6,shape(1,3),(*this)(row,col),shape(1,4),0.0,tmp7);

      for(int i = 0;i < delta;++i){

         //contract the peps with right operators: peps_i's
         peps_op.clear();
         Contract(1.0,ham.gR(i),shape(1),(*this)(row,col),shape(2),0.0,peps_op);

         //close down the LOi's with peps_i's
         tmp6.clear();
         Contract(1.0,tmp7,shape(1,2,5),peps_op,shape(2,4,0),0.0,tmp6);

         tmp6bis.clear();
         Permute(tmp6,shape(0,2,4,3,5,1),tmp6bis);

         RO[col].clear();
         Gemm(CblasNoTrans,CblasTrans,1.0,tmp6bis,env.gb(row - 1)[col],0.0,RO[col]);

         //expectation value:
         val += ham.gcoef(i) * Dot(LOi[i],RO[col]);

      }

      //add local term
      if(ham.gis_local()){

         //contract the peps with right operators: peps_i's
         peps_op.clear();
         Contract(1.0,ham.gB(),shape(1),(*this)(row,col),shape(2),0.0,peps_op);

         //close down the left unit with peps_op
         tmp6.clear();
         Contract(1.0,tmp7,shape(1,2,5),peps_op,shape(2,4,0),0.0,tmp6);

         tmp6bis.clear();
         Permute(tmp6,shape(0,2,4,3,5,1),tmp6bis);

         Gemm(CblasNoTrans,CblasTrans,1.0,tmp6bis,env.gb(row - 1)[col],0.0,RO[col]);

         //expectation value:
         val += Dot(RO[col],RO[col-1]);

      }

      // now construct the new left going renormalized operators

      //first attach top to left unity
      tmp6.clear();
      Contract(1.0,env.gt(row)[col],shape(0),RO[col-1],shape(0),0.0,tmp6);

      //add peps to it, intermediary
      tmp7.clear();
      Contract(1.0,tmp6,shape(3,0),(*this)(row,col),shape(0,1),0.0,tmp7);

      //paste left operators to intermediary to make new left operators
      for(int i = 0;i < delta;++i){

         Contract(1.0,ham.gL(i),shape(1),(*this)(row,col),shape(2),0.0,peps_op);

         tmp6.clear();
         Contract(1.0,tmp7,shape(4,2,0),peps_op,shape(0,1,2),0.0,tmp6);

         tmp6bis.clear();
         Permute(tmp6,shape(0,3,5,1,2,4),tmp6bis);

         LOi[i].clear();
         Gemm(CblasNoTrans,CblasNoTrans,1.0,tmp6bis,env.gb(row - 1)[col],0.0,LOi[i]);

      }

      //finally construct new left unity
      Contract(1.0,tmp7,shape(2,0,4),(*this)(row,col),shape(0,1,2),0.0,tmp6);
      Permute(tmp6,shape(0,3,5,1,2,4),tmp6bis);

      RO[col].clear();
      Gemm(CblasNoTrans,CblasNoTrans,1.0,tmp6bis,env.gb(row - 1)[col],0.0,RO[col]);

   }

   //last site on the right: close down on the incomings
   for(int i = 0;i < delta;++i){

      peps_op.clear();
      Contract(1.0,ham.gR(i),shape(1),(*this)(row,Lx - 1),shape(2),0.0,peps_op);

      tmp6.clear();
      Contract(1.0,env.gt(row)[Lx - 1],shape(0),LOi[i],shape(0),0.0,tmp6);

      //add peps to it
      tmp7.clear();
      Contract(1.0,tmp6,shape(3,0),(*this)(row,Lx - 1),shape(0,1),0.0,tmp7);

      tmp6.clear();
      Contract(1.0,tmp7,shape(4,2,0),peps_op,shape(0,1,2),0.0,tmp6);

      val += ham.gcoef(i) * blas::dot(tmp6.size(), tmp6.data(), 1, env.gb(row-1)[Lx-1].data(), 1);

   }

   if(ham.gis_local()){

      peps_op.clear();
      Contract(1.0,ham.gB(),shape(1),(*this)(row,Lx - 1),shape(2),0.0,peps_op);

      tmp6.clear();
      Contract(1.0,env.gt(row)[Lx - 1],shape(0),RO[Lx - 2],shape(0),0.0,tmp6);

      //add peps to it
      tmp7.clear();
      Contract(1.0,tmp6,shape(3,0),(*this)(row,Lx - 1),shape(0,1),0.0,tmp7);

      tmp6.clear();
      Contract(1.0,tmp7,shape(4,2,0),peps_op,shape(0,1,2),0.0,tmp6);

      val += blas::dot(tmp6.size(), tmp6.data(), 1, env.gb(row-1)[Lx-1].data(), 1);

   }

   //   }

   // -- (3) -- || top row = Ly-1: again similar to overlap calculation

   //first construct the right renormalized operators
   contractions::init_ro(ham.gis_local(),'t',*this,R);

   //construct the left operator with two open physical bonds
   tmp5.clear();
   Contract(1.0,(*this)(Ly-1,0),shape(0,3),env.gb(Ly-2)[0],shape(0,2),0.0,tmp5);

   tmp5bis.clear();
   Permute(tmp5,shape(1,0,3,2,4),tmp5bis);

   //add local energy term
   if(ham.gis_local()){

      peps_op.clear();
      Contract(1.0,ham.gB(),shape(1),(*this)(Ly-1,0),shape(2),0.0,peps_op);

      tmp4.clear();
      Contract(1.0,peps_op,shape(0,2,3),tmp5bis,shape(0,1,2),0.0,tmp4);

      val += blas::dot(tmp4.size(), tmp4.data(), 1, R[0].data(), 1);

   }

   //multiply intermediate with Left hamiltonian operators
   for(int i = 0;i < delta;++i){

      peps_op.clear();
      Contract(1.0,ham.gL(i),shape(1),(*this)(Ly-1,0),shape(2),0.0,peps_op);

      tmp4.clear();
      Contract(1.0,peps_op,shape(0,2,3),tmp5bis,shape(0,1,2),0.0,tmp4);

      Li[i] = tmp4.reshape_clear(shape(tmp4.shape(1),tmp4.shape(2),tmp4.shape(3)));

   }

   //and finally unity
   Contract(1.0,(*this)(Ly-1,0),shape(2,1,3),tmp5bis,shape(0,1,2),0.0,tmp4);
   R[0] = tmp4.reshape_clear(shape(tmp4.shape(1),tmp4.shape(2),tmp4.shape(3)));

   //middle of the chain:
   for(int col = 1;col < Lx-1;++col){

      //close down the left renormalized terms from the previous site

      //construct the intermediate contraction (paste top right)
      tmp5.clear();
      Contract(1.0,env.gb(Ly-2)[col],shape(3),R[col],shape(2),0.0,tmp5);

      //and upper peps to tmp5
      tmp6.clear();
      Contract(1.0,(*this)(Ly-1,col),shape(3,4),tmp5,shape(2,4),0.0,tmp6);

      //paste right hamiltonian operators to intermediary an contract with left renormalized operator
      for(int i = 0;i < delta;++i){

         peps_op.clear();
         Contract(1.0,ham.gR(i),shape(1),(*this)(Ly - 1,col),shape(2),0.0,peps_op);

         tmp5.clear();
         Contract(1.0,peps_op,shape(0,3,4),tmp6,shape(2,4,5),0.0,tmp5);

         //contract with left hamiltonian operator
         val += ham.gcoef(i) * blas::dot(Li[i].size(), Li[i].data(), 1, tmp5.data(), 1);

      }

      //add local hamiltonian term
      if(ham.gis_local()){

         peps_op.clear();
         Contract(1.0,ham.gB(),shape(1),(*this)(Ly - 1,col),shape(2),0.0,peps_op);

         tmp5.clear();
         Contract(1.0,peps_op,shape(0,3,4),tmp6,shape(2,4,5),0.0,tmp5);

         //contract with left hamiltonian operator
         val += blas::dot(R[col-1].size(), R[col-1].data(), 1, tmp5.data(), 1);

      }

      //construct left renormalized operators for next site: first construct intermediary
      tmp5.clear();
      Contract(1.0,env.gb(Ly-2)[col],shape(0),R[col-1],shape(2),0.0,tmp5);

      tmp6.clear();
      Contract(1.0,(*this)(Ly-1,col),shape(3,0),tmp5,shape(1,4),0.0,tmp6);

      //now paste on left hamiltonian operators 
      for(int i = 0;i < delta;++i){

         Contract(1.0,ham.gL(i),shape(1),(*this)(Ly - 1,col),shape(2),0.0,peps_op);

         tmp5.clear();
         Contract(1.0,peps_op,shape(0,3,1),tmp6,shape(1,3,5),0.0,tmp5);

         Li[i] = tmp5.reshape_clear(shape((*this)(Ly-1,col).shape(4),(*this)(Ly-1,col).shape(4),env.gb(Ly-2)[col].shape(3)));

      }

      //finally construct new unity on the left
      tmp5.clear();
      Contract(1.0,(*this)(Ly-1,col),shape(2,3,0),tmp6,shape(1,3,5),0.0,tmp5);

      R[col] = tmp5.reshape_clear(shape((*this)(Ly-1,col).shape(4),(*this)(Ly-1,col).shape(4),env.gb(Ly-2)[col].shape(3)));

   }

   //finally close down on last top site

   //construct intermediate
   tmp7.clear();
   Contract(1.0,(*this)(Ly-1,Lx-1),shape(3),env.gb(Ly-2)[Lx - 1],shape(2),0.0,tmp7);

   // close down Li[i] right hamiltonian operators
   for(int i = 0;i < delta;++i){

      peps_op.clear();
      Contract(1.0,ham.gR(i),shape(1),(*this)(Ly - 1,Lx - 1),shape(2),0.0,peps_op);

      tmp8.clear();
      Contract(1.0,peps_op,shape(0,3),tmp7,shape(2,5),0.0,tmp8);

      val += ham.gcoef(i) * blas::dot(Li[i].size(),Li[i].data(),1,tmp8.data(),1);

   }

   //close down local term
   if(ham.gis_local()){

      peps_op.clear();
      Contract(1.0,ham.gB(),shape(1),(*this)(Ly - 1,Lx - 1),shape(2),0.0,peps_op);

      tmp8.clear();
      Contract(1.0,peps_op,shape(0,3),tmp7,shape(2,5),0.0,tmp8);

      val += blas::dot(R[Lx - 2].size(),R[Lx - 2].data(),1,tmp8.data(),1);

   }

   // #################################################################
   // ###   ---- from right left : contract in mps/mpo fashion ---- ### 
   // #################################################################

   // -- (1) -- || right column: similar to overlap calculation

   //first construct the right renormalized operators
   R.resize(Ly - 1);

   contractions::init_ro(false,'r',*this,R);

   tmp5.clear();
   Contract(1.0,env.gl(Lx - 2)[0],shape(0,1),(*this)(0,Lx - 1),shape(3,0),0.0,tmp5);

   tmp5bis.clear();
   Permute(tmp5,shape(1,2,3,0,4),tmp5bis);

   //multiply intermediate with Left hamiltonian operators
   for(int i = 0;i < delta;++i){

      peps_op.clear();
      Contract(1.0,ham.gL(i),shape(1),(*this)(0,Lx - 1),shape(2),0.0,peps_op);

      tmp4.clear();
      Contract(1.0,tmp5bis,shape(2,3,4),peps_op,shape(0,1,4),0.0,tmp4);

      Li[i] = tmp4.reshape_clear(shape(tmp4.shape(0),tmp4.shape(1),tmp4.shape(2)));

   }

   //left unit
   tmp4.clear();
   Contract(1.0,tmp5bis,shape(2,3,4),(*this)(0,Lx-1),shape(2,0,4),0.0,tmp4);

   R[0] = tmp4.reshape_clear(shape(tmp4.shape(0),tmp4.shape(1),tmp4.shape(2)));

   //now for the middle terms
   for(int row = 1;row < Ly - 1;++row){

      //first close down the interaction terms from the previous site

      //construct the intermediate contraction (paste top right)
      tmp5.clear();
      Contract(1.0,env.gl(Lx - 2)[row],shape(3),R[row],shape(0),0.0,tmp5);

      //and upper peps to tmp5
      tmp6.clear();
      Contract(1.0,tmp5,shape(1,3),(*this)(row,Lx-1),shape(0,1),0.0,tmp6);

      //paste right hamiltonian operator to intermediary and constact with left renormalized operator
      for(int i = 0;i < delta;++i){

         peps_op.clear();
         Contract(1.0,ham.gR(i),shape(1),(*this)(row,Lx-1),shape(2),0.0,peps_op);

         tmp5.clear();
         Contract(1.0,tmp6,shape(3,1,2),peps_op,shape(0,1,2),0.0,tmp5);

         //contract with left S+
         val += ham.gcoef(i) * blas::dot(Li[i].size(),Li[i].data(),1,tmp5.data(),1);

      }

      //construct left renormalized operators for next site: first construct intermediary
      tmp5.clear();
      Contract(1.0,R[row-1],shape(0),env.gl(Lx - 2)[row],shape(0),0.0,tmp5);

      tmp6.clear();
      Contract(1.0,tmp5,shape(0,2),(*this)(row,Lx - 1),shape(3,0),0.0,tmp6);

      // then contract with left Hamiltonian operator
      for(int i = 0;i < delta;++i){

         Contract(1.0,ham.gL(i),shape(1),(*this)(row,Lx-1),shape(2),0.0,peps_op);

         tmp5.clear();
         Contract(1.0,tmp6,shape(0,1,4),peps_op,shape(3,1,0),0.0,tmp5);

         Li[i] = tmp5.reshape_clear(shape(env.gl(Lx - 2)[row].shape(3),(*this)(row,Lx - 1).shape(1),(*this)(row,Lx - 1).shape(1)));

      }

      // finally construct new unity on the left
      Contract(1.0,tmp6,shape(0,1,4),(*this)(row,Lx-1),shape(3,0,2),0.0,tmp5);
      R[row] = tmp5.reshape_clear(shape(env.gl(Lx - 2)[row].shape(3),(*this)(row,Lx - 1).shape(1),(*this)(row,Lx - 1).shape(1)));

   }

   //last site: close down the left operators with right hamiltonian operators

   //construct intermediate first
   tmp7.clear();
   Contract(1.0,env.gl(Ly-2)[Lx - 1],shape(1),(*this)(Ly-1,Lx-1),shape(0),0.0,tmp7);

   //contract with right operator for energy
   for(int i = 0;i < delta;++i){

      peps_op.clear();
      Contract(1.0,ham.gR(i),shape(1),(*this)(Ly-1,Lx-1),shape(2),0.0,peps_op);

      tmp8.clear();
      Contract(1.0,tmp7,shape(4,1),peps_op,shape(0,1),0.0,tmp8);

      val += ham.gcoef(i) * blas::dot(Li[i].size(), Li[i].data(), 1, tmp8.data(), 1);

   }

   // -- (2) -- now move from right to left calculating everything like an MPO/MPS expectation value

   //Right renormalized operators
   RO.resize(Ly - 1);

   //loop over the columns
   for(int col = Lx - 2;col > 0;--col){

      //first create right renormalized operator
      contractions::init_ro(false,'V',col,*this,RO);

      // --- now move from left to right to get the expecation value of the interactions ---

      // --- First construct the left going operators for the first site -----

      //paste left environment on
      tmp7.clear();
      Contract(1.0,env.gl(col - 1)[0],shape(1),(*this)(0,col),shape(0),0.0,tmp7);

      // add left hamiltonian operators to intermediate
      for(int i = 0;i < delta;++i){

         peps_op.clear();
         Contract(1.0,ham.gL(i),shape(1),(*this)(0,col),shape(2),0.0,peps_op);

         tmp8.clear();
         Contract(1.0,tmp7,shape(4,1),peps_op,shape(0,1),0.0,tmp8);

         tmp8bis.clear();
         Contract(1.0,tmp8,shape(4,7),env.gr(col)[0],shape(1,2),0.0,tmp8bis);

         //move to a DArray<3> object: order (top-env,(*this)-row,bottom-env)
         LOi[i] = tmp8bis.reshape_clear(shape(env.gl(col - 1)[0].shape(3),(*this)(0,col).shape(1),(*this)(0,col).shape(1),env.gr(col)[0].shape(3)));

      }

      // construct left renormalized operator with unity
      Contract(1.0,tmp7,shape(1,4),(*this)(0,col),shape(0,2),0.0,tmp8);
      Contract(1.0,tmp8,shape(4,7),env.gr(col)[0],shape(1,2),0.0,tmp8bis);

      //move to a DArray<3> object: order (top-env,(*this)-row,bottom-env)
      RO[0] = tmp8bis.reshape_clear(shape(env.gl(col - 1)[0].shape(3),(*this)(0,col).shape(1),(*this)(0,col).shape(1),env.gr(col)[0].shape(3)));

      // --- now for the middle sites, close down the operators on the left and construct new ones --- 
      for(int row = 1;row < Ly - 1;++row){

         //first add 'left' and one peps to the right side to construct intermediate
         tmp6.clear();
         Contract(1.0,env.gl(col - 1)[row],shape(3),RO[row],shape(0),0.0,tmp6);

         tmp7.clear();
         Contract(1.0,tmp6,shape(1,3),(*this)(row,col),shape(0,1),0.0,tmp7);

         // close down LOp with right hamiltonian operators pasted on intermediate for energy contribution
         for(int i = 0;i < delta;++i){

            peps_op.clear();
            Contract(1.0,ham.gR(i),shape(1),(*this)(row,col),shape(2),0.0,peps_op);

            tmp6.clear();
            Contract(1.0,tmp7,shape(4,1,2),peps_op,shape(0,1,2),0.0,tmp6);

            tmp6bis.clear();
            Permute(tmp6,shape(0,2,4,3,5,1),tmp6bis);

            RO[row].clear();
            Gemm(CblasNoTrans,CblasTrans,1.0,tmp6bis,env.gr(col)[row],0.0,RO[row]);

            //expectation value:
            val += ham.gcoef(i) * Dot(LOi[i],RO[row]);

         }

         // now construct the new left going renormalized operators

         //first attach top to left unity, make intermediary
         tmp6.clear();
         Contract(1.0,env.gl(col - 1)[row],shape(0),RO[row-1],shape(0),0.0,tmp6);

         //add peps to it
         tmp7.clear();
         Contract(1.0,tmp6,shape(0,3),(*this)(row,col),shape(0,3),0.0,tmp7);

         // add left hamiltonian operator to intermediary to construct next left renormalized operator
         for(int i = 0;i < delta;++i){

            Contract(1.0,ham.gL(i),shape(1),(*this)(row,col),shape(2),0.0,peps_op);

            tmp6.clear();
            Contract(1.0,tmp7,shape(5,0,2),peps_op,shape(0,1,3),0.0,tmp6);

            tmp6bis.clear();
            Permute(tmp6,shape(0,2,4,1,3,5),tmp6bis);

            LOi[i].clear();
            Gemm(CblasNoTrans,CblasNoTrans,1.0,tmp6bis,env.gr(col)[row],0.0,LOi[i]);

         }

         // finally construct new left unity
         Contract(1.0,tmp7,shape(0,5,2),(*this)(row,col),shape(0,2,3),0.0,tmp6);
         Permute(tmp6,shape(0,2,4,1,3,5),tmp6bis);

         RO[row].clear();
         Gemm(CblasNoTrans,CblasNoTrans,1.0,tmp6bis,env.gr(col)[row],0.0,RO[row]);

      }

      //last site on the right: close down on the incomings
      for(int i = 0;i < delta;++i){

         //add right hamiiltonian operator to peps
         peps_op.clear();
         Contract(1.0,ham.gR(i),shape(1),(*this)(Ly-1,col),shape(2),0.0,peps_op);

         //make intermediary
         tmp6.clear();
         Contract(1.0,env.gl(col - 1)[Ly - 1],shape(0),LOi[i],shape(0),0.0,tmp6);

         tmp7.clear();
         Contract(1.0,tmp6,shape(0,3),(*this)(Ly - 1,col),shape(0,3),0.0,tmp7);

         //contract with right ham operator
         tmp6.clear();
         Contract(1.0,tmp7,shape(0,2,5),peps_op,shape(1,3,0),0.0,tmp6);

         val += ham.gcoef(i) * blas::dot(tmp6.size(), tmp6.data(), 1, env.gr(col)[Ly-1].data(), 1);

      }

   }

   // -- (3) -- || left column = 0: again similar to overlap calculation

   //first construct the right renormalized operators
   contractions::init_ro(false,'l',*this,R);

   //construct the left operator with two open physical bonds
   tmp5.clear();
   Contract(1.0,(*this)(0,0),shape(3,4),env.gr(0)[0],shape(0,2),0.0,tmp5);

   tmp5bis.clear();
   Permute(tmp5,shape(2,0,3,1,4),tmp5bis);

   //contract with left hamiltonian operators to construct left renormalized operators
   for(int i = 0;i < delta;++i){

      peps_op.clear();
      Contract(1.0,ham.gL(i),shape(1),(*this)(0,0),shape(2),0.0,peps_op);

      tmp4.clear();
      Contract(1.0,peps_op,shape(0,1,4),tmp5bis,shape(0,1,2),0.0,tmp4);

      Li[i] = tmp4.reshape_clear( shape( tmp4.shape(0),tmp4.shape(2),tmp4.shape(3) ) );

   }

   //and contract with peps to get left unit operator
   tmp4.clear();
   Contract(1.0,(*this)(0,0),shape(2,0,4),tmp5bis,shape(0,1,2),0.0,tmp4);

   R[0] = tmp4.reshape_clear( shape( tmp4.shape(0),tmp4.shape(2),tmp4.shape(3) ) );

   //middle of the chain:
   for(int row = 1;row < Ly-1;++row){

      //close down the left renormalized operators from the previous site

      //construct the intermediate contraction (paste right)
      tmp5.clear();
      Contract(1.0,env.gr(0)[row],shape(3),R[row],shape(2),0.0,tmp5);

      //and upper peps to tmp5
      tmp6.clear();
      Contract(1.0,(*this)(row,0),shape(4,1),tmp5,shape(2,4),0.0,tmp6);

      //paste right hamiltonian operators to it for energy 
      for(int i = 0;i < delta;++i){

         peps_op.clear();
         Contract(1.0,ham.gR(i),shape(1),(*this)(row,0),shape(2),0.0,peps_op);

         tmp5.clear();
         Contract(1.0,peps_op,shape(0,4,2),tmp6,shape(1,4,5),0.0,tmp5);

         //contract with left S+
         val += ham.gcoef(i) * blas::dot(Li[i].size(), Li[i].data(), 1, tmp5.data(), 1);

      }

      //construct left renormalized operators for next site: first construct intermediary
      tmp5.clear();
      Contract(1.0,env.gr(0)[row],shape(0),R[row-1],shape(2),0.0,tmp5);

      tmp6.clear();
      Contract(1.0,(*this)(row,0),shape(3,4),tmp5,shape(4,1),0.0,tmp6);

      // paste left hamiltonian operator to it
      for(int i = 0;i < delta;++i){

         Contract(1.0,ham.gL(i),shape(1),(*this)(row,0),shape(2),0.0,peps_op);

         tmp5.clear();
         Contract(1.0,peps_op,shape(0,3,4),tmp6,shape(2,5,3),0.0,tmp5);

         Li[i] = tmp5.reshape_clear(shape((*this)(row,0).shape(1),(*this)(row,0).shape(1),env.gr(0)[row].shape(3)));

      }

      // construct new unity on the left
      Contract(1.0,(*this)(row,0),shape(2,3,4),tmp6,shape(2,5,3),0.0,tmp5);

      R[row] = tmp5.reshape_clear(shape((*this)(row,0).shape(1),(*this)(row,0).shape(1),env.gr(0)[row].shape(3)));

   }

   //finally close down on last left site

   //first construct intermediate
   tmp7.clear();
   Contract(1.0,(*this)(Ly-1,0),shape(4),env.gr(0)[Ly - 1],shape(2),0.0,tmp7);

   // paste on right hamiltonian operator to close down left renormalized operator
   for(int i = 0;i < delta;++i){

      peps_op.clear();
      Contract(1.0,ham.gR(i),shape(1),(*this)(Ly-1,0),shape(2),0.0,peps_op);

      tmp8.clear();
      Contract(1.0,peps_op,shape(0,4),tmp7,shape(2,5),0.0,tmp8);

      val += ham.gcoef(i) * blas::dot(Li[i].size(),Li[i].data(),1,tmp8.data(),1);

   }
   */
      return val;

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

template void PEPS<double>::save(const char *filename);
template void PEPS< complex<double> >::save(const char *filename);
