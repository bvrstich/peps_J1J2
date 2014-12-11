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

/** 
 * empty constructor
 */
template<typename T>
MPS<T>::MPS() : vector< TArray<T,3> >() { }

/** 
 * constructor: just sets the length of the vector, nothing is allocates or initialized
 * @param L_in length of the chain
 */
template<typename T>
MPS<T>::MPS(int L_in) : vector< TArray<T,3> >(L_in) { }

/** 
 * standard constructor:
 * @param L_in length of the chain
 * @param d_phys_in physical dimension
 * @param D_in virtual max bond dimension
 * allocates the tensors and fills them randomly
 */
template<typename T>
MPS<T>::MPS(int L_in,int d_phys_in,int D_in) : vector< TArray<T,3> >(L_in) {

   D = D_in;
   d_phys = d_phys_in;

   vector<int> vdim(L_in + 1);

   vdim[0] = 1;

   for(int i = 1;i < L_in;++i){

      int tmp = vdim[i - 1] * d_phys;

      if(tmp < D)
         vdim[i] = tmp;
      else 
         vdim[i] = D;

   }

   vdim[L_in] = 1;

   for(int i = L_in - 1;i > 0;--i){

      int tmp = vdim[i + 1] * d_phys;

      if(tmp < vdim[i])
         vdim[i] = tmp;

   }

   for(int i = 0;i < this->size();++i){

      (*this)[i].resize(vdim[i],d_phys,vdim[i+1]);
      (*this)[i].generate(global::rgen<T>);

   }

}

/**
 * copy constructor
 */
template<typename T>
MPS<T>::MPS(const MPS<T> &mps_copy) : vector< TArray<T,3> >(mps_copy) {

   this->D = mps_copy.gD();
   this->d_phys = mps_copy.gd_phys();

}

/**
 * empty destructor
 */
template<typename T>
MPS<T>::~MPS(){ }

/**
 * @return virtual dimension of the MPS
 */
template<typename T>
int MPS<T>::gD() const {

   return D;

}

/**
 * @return virtual dimension of the MPS
 */
template<typename T>
int MPS<T>::gd_phys() const {

   return d_phys;

}


/**
 * act with an MPO on this MPS, resulting MPS is returned as *this object
 * @param uplo if == 'U' contract with the upper physical index of the MPO, if == 'L', contract with the lower
 * @param mpo the MPO
 */
template<typename T>
void MPS<T>::gemv(char uplo,const MPO<T> &mpo){

   int DO = mpo.gD();

   int L = this->size();

   if(uplo == 'U'){

      //first site
      TArray<T,5> tmp;

      enum {i,j,k,l,m,n,o,p,q};

      //dimensions of the new MPS
      int DL = 1;
      int DR = (*this)[0].shape(2) * mpo[0].shape(3);

      Contract((T)1.0,mpo[0],shape(i,j,k,l),(*this)[0],shape(m,j,n),(T)0.0,tmp,shape(m,i,k,n,l));

      (*this)[0] = tmp.reshape_clear(shape(DL,d_phys,DR));

      //middle sites
      for(int c = 1;c < L - 1;++c){

         DL = DR;
         DR = (*this)[c].shape(2) * mpo[c].shape(3);

         Contract((T)1.0,mpo[c],shape(i,j,k,l),(*this)[c],shape(m,j,n),(T)0.0,tmp,shape(m,i,k,n,l));

         (*this)[c] = tmp.reshape_clear(shape(DL,d_phys,DR));

      }

      DL = DR;
      DR = 1;

      Contract((T)1.0,mpo[L - 1],shape(i,j,k,l),(*this)[L - 1],shape(m,j,n),(T)0.0,tmp,shape(m,i,k,n,l));

      (*this)[L - 1] = tmp.reshape_clear(shape(DL,d_phys,DR));

   }
   else{//L

      //first site
      TArray<T,5> tmp;

      enum {i,j,k,l,m,n,o,p,q};

      //dimensions of the new MPS
      int DL = 1;
      int DR = (*this)[0].shape(2) * mpo[0].shape(3);

      Contract((T)1.0,mpo[0],shape(i,j,k,l),(*this)[0],shape(m,k,n),(T)0.0,tmp,shape(m,i,j,n,l));

      (*this)[0] = tmp.reshape_clear(shape(DL,d_phys,DR));

      //middle sites
      for(int c = 1;c < L - 1;++c){

         DL = DR;
         DR = (*this)[c].shape(2) * mpo[c].shape(3);

         Contract((T)1.0,mpo[c],shape(i,j,k,l),(*this)[c],shape(m,k,n),(T)0.0,tmp,shape(m,i,j,n,l));

         (*this)[c] = tmp.reshape_clear(shape(DL,d_phys,DR));

      }

      DL = DR;
      DR = 1;

      Contract((T)1.0,mpo[L - 1],shape(i,j,k,l),(*this)[L - 1],shape(m,k,n),(T)0.0,tmp,shape(m,i,j,n,l));

      (*this)[L - 1] = tmp.reshape_clear(shape(DL,d_phys,DR));

   }

   //vdim is increased
   D *= DO;

}

/**
 * canonicalize the mps
 * @param dir Left or Right canonicalization
 * @param norm if true: normalize, else not
 */
template<typename T>
void MPS<T>::canonicalize(const BTAS_SIDE &dir,bool norm){

   int length = this->size();

   if(dir == Left){//QR

      TArray<T,2> R;
      TArray<T,3> tmp;

      for(int i = 0;i < length - 1;++i){

         R.clear();

         //do QR
         Geqrf((*this)[i],R);

         //paste to next matrix
         tmp.clear();

         Contract((T)1.0,R,shape(1),(*this)[i + 1],shape(0),(T)0.0,tmp);

         (*this)[i + 1] = std::move(tmp);

      }

      if(norm){

         T nrm = sqrt(Dotc((*this)[length-1],(*this)[length-1]));
         Scal(1.0/nrm,(*this)[length-1]);

      }

   }
   else{//LQ

      TArray<T,2> L;
      TArray<T,3> tmp;

      for(int i = length - 1;i > 0;--i){

         L.clear();

         //do QR
         Gelqf(L,(*this)[i]);

         //paste to previous matrix
         tmp.clear();

         Contract((T)1.0,(*this)[i - 1],shape(2),L,shape(0),(T)0.0,tmp);

         (*this)[i - 1] = std::move(tmp);

      }

      if(norm){

         T nrm = sqrt(Dotc((*this)[0],(*this)[0]));
         Scal(1.0/nrm,(*this)[0]);

      }

   }

}

/**
 * scale the MPS with a constant factor
 * @param alpha scalingfactor
 */
template<>
void MPS<double>::scal(double alpha){

   int sign;

   if(alpha > 0)
      sign = 1;
   else
      sign = -1;

   alpha = pow(fabs(alpha),1.0/(double)this->size());

   Scal(sign * alpha,(*this)[0]);

   for(int i = 1;i < this->size();++i)
      Scal(alpha,(*this)[i]);

}

/**
 * scale the MPS with a constant factor
 * @param alpha scalingfactor
 */
template<>
void MPS< complex<double> >::scal(complex<double> alpha){

   alpha = pow(fabs(alpha),1.0/(complex<double>)this->size());

   Scal(alpha,(*this)[0]);

   for(int i = 1;i < this->size();++i)
      Scal(alpha,(*this)[i]);

}

/**
 * @param bra the bra of the inner product
 * @return the inner product of two MPS's, with *this being the ket
 */
template<typename T>
T MPS<T>::dot(const MPS<T> &bra) const {

   TArray<T,2> E;

   Contract((T)1.0,bra[0],shape(0,1),(*this)[0],shape(0,1),(T)0.0,E);

   TArray<T,3> I;

   for(int i = 1;i < this->size();++i){

      I.clear();

      Contract((T)1.0,bra[i],shape(0),E,shape(0),(T)0.0,I);

      E.clear();

      Contract((T)1.0,I,shape(2,0),(*this)[i],shape(0,1),(T)0.0,E);

   }

   return E(0,0);

}

/**
 * normalize the mps
 * @return the norm before
 */
template<typename T>
T MPS<T>::normalize(){

   T nrm = std::sqrt(this->dot(*this));

   this->scal(1.0/nrm);

   return nrm;

}

template MPS<double>::MPS(int,int,int);
template MPS< complex<double> >::MPS(int,int,int);

template MPS<double>::MPS(int);
template MPS< complex<double> >::MPS(int);

template MPS<double>::MPS();
template MPS< complex<double> >::MPS();

template MPS<double>::MPS(const MPS<double> &);
template MPS< complex<double> >::MPS(const MPS< complex<double> > &);

template MPS<double>::~MPS();
template MPS< complex<double> >::~MPS();

template int MPS<double>::gD() const;
template int MPS< complex<double> >::gD() const;

template void MPS<double>::gemv(char uplo,const MPO<double> &mpo);
template void MPS< complex<double> >::gemv(char uplo,const MPO< complex<double> > &mpo);

template void MPS<double>::canonicalize(const BTAS_SIDE &dir,bool);
template void MPS< complex<double> >::canonicalize(const BTAS_SIDE &dir,bool);

template double MPS<double>::dot(const MPS<double> &bra) const;
template  complex<double>  MPS< complex<double> >::dot(const MPS< complex<double> > &bra) const;

template double MPS<double>::normalize();
template complex<double> MPS< complex<double> >::normalize();
