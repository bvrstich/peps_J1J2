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
template<typename T>
MPO<T>::MPO() : vector< TArray<T,4> >() { }

/** 
 * constructor: just sets the length of the vector, nothing is allocates or initialized
 * @param L_in length of the chain
 */
template<typename T>
MPO<T>::MPO(int L_in) : vector< TArray<T,4> >(L_in) { }

/** 
 * standard constructor:
 * @param L_in length of the chain
 * @param d_phys_in physical dimension
 * @param D_in virtual max bond dimension
 * allocates the tensors with correct dimensions
 */
template<typename T>
MPO<T>::MPO(int L_in,int d_phys_in,int D_in) : vector< TArray<T,4> >(L_in) {

   D = D_in;
   d_phys = d_phys_in;

   vector<int> vdim(L_in + 1);

   vdim[0] = 1;

   for(int i = 1;i < L_in;++i){

      int tmp = vdim[i - 1] * d_phys * d_phys;

      if(tmp < D)
         vdim[i] = tmp;
      else 
         vdim[i] = D;

   }

   vdim[L_in] = 1;

   for(int i = L_in - 1;i > 0;--i){

      int tmp = vdim[i + 1] * d_phys * d_phys;

      if(tmp < vdim[i])
         vdim[i] = tmp;

   }

   for(int i = 0;i < this->size();++i)
      (*this)[i].resize(vdim[i],d_phys,d_phys,vdim[i+1]);

}

/**
 * copy constructor
 */
template<typename T>
MPO<T>::MPO(const MPO<T> &mpo_copy) : vector< TArray<T,4> >(mpo_copy) {

   D = mpo_copy.gD();
   d_phys = mpo_copy.gd_phys();

}

/**
 * empty destructor
 */
template<typename T>
MPO<T>::~MPO(){ }

/**
 * @return virtual dimension of the MPO
 */
template<typename T>
int MPO<T>::gD() const {

   return D;

}

/**
 * @return physical dimension of the MPO
 */
template<typename T>
int MPO<T>::gd_phys() const {

   return d_phys;

}

/**
 * @param bra the bra of the inner product
 * @return the inner product of two MPO's, with *this being the ket
 */
template<typename T>
T MPO<T>::dot(const MPO<T> &bra) const {

   TArray<T,2> E;

   Contract((T)1.0,bra[0],shape(0,1,2),(*this)[0],shape(0,1,2),(T)0.0,E);

   TArray<T,4> I;

   for(int i = 1;i < this->size();++i){

      I.clear();

      Contract((T)1.0,E,shape(0),bra[i],shape(0),(T)0.0,I);

      E.clear();

      Contract((T)1.0,I,shape(0,1,2),(*this)[i],shape(0,1,2),(T)0.0,E);

   }

   return E(0,0);

}

/**
 * @param option 'H'orizontal or 'V'ertical
 * @param rc the row/col index of the peps
 * @param peps the peps 'operator'
 * @return the expectation value of the double peps row within the boundary MPO
 */
template<typename T>
T MPO<T>::expect(const char option,int rc,const PEPS<T> &peps) const {

   TArray<T,4> E;

   if(option == 'H'){

      TArray<T,7> tmp7;
      Contract((T)1.0,(*this)[0],shape(1),peps(rc,0),shape(1),(T)0.0,tmp7);

      TArray<T,8> tmp8;
      Contract((T)1.0,tmp7,shape(1,4),peps(rc,0),shape(1,2),(T)0.0,tmp8);

      TArray<T,8> tmp8bis;
      Contract((T)1.0,tmp8,shape(3,6),(*this)[0],shape(1,2),(T)0.0,tmp8bis);

      E = tmp8bis.reshape_clear(shape(tmp8bis.shape(1),tmp8bis.shape(3),tmp8bis.shape(5),tmp8bis.shape(7)));

      //now for the rest of the rightgoing sweep.
      for(int i = 1;i < Lx;++i){

         TArray<T,6> tmp6;
         Contract((T)1.0,E,shape(0),(*this)[i],shape(0),(T)0.0,tmp6);

         tmp7.clear();
         Contract((T)1.0,tmp6,shape(0,3),peps(rc,i),shape(0,1),(T)0.0,tmp7);

         tmp6.clear();
         Contract((T)1.0,tmp7,shape(0,2,4),peps(rc,i),shape(0,1,2),(T)0.0,tmp6);

         E.clear();
         Contract((T)1.0,tmp6,shape(0,2,4),(*this)[i],shape(0,1,2),(T)0.0,E);

      }

   }
   else{//Vertical

   }

   return E(0,0,0,0);

}


/**
 * Fill the MPO by contracting a the physical dimensions of a peps object
 * @param option 'b'ottom, 't'op, 'l'eft or 'r'ight
 * @param peps input PEPS<double> object
 */
template<>
void MPO<double>::fill(const char option,const PEPS<double> &peps){

   if(option == 'b'){

      enum {i,j,k,l,m,n,o,p,q};

      DArray<8> tmp;

      //share the pointer
      for(int col = 0;col < Lx;col++){

         tmp.share_mem( (*this)[col] );

         Contract(1.0,peps(0,col),shape(i,j,k,l,m),peps(0,col),shape(n,o,k,p,q),0.0,tmp,shape(i,n,j,o,l,p,m,q));

      }

      //canonicalize

   }
   else{

      enum {i,j,k,l,m,n,o,p,q};

      DArray<8> tmp;

      //share the pointer
      for(int col = 0;col < Lx;col++){

         tmp.share_mem( (*this)[col] );

         Contract(1.0,peps(Ly - 1,col),shape(i,j,k,l,m),peps(Ly - 1,col),shape(n,o,k,p,q),0.0,tmp,shape(i,n,j,o,l,p,m,q));

      }

   }

}
/**
 * Fill the MPO with random entries, output is right normalized state
 */
template<>
void MPO<double>::fill_Random() {

   DArray<2> L;

   for(int i = this->size() - 1;i >= 0;--i){

      (*this)[i].generate(rgen<double>);
      Gelqf(L,(*this)[i]);

   }

}

/**
 * scale the MPO with a constant factor
 * @param alpha scalingfactor
 */
template<>
void MPO<double>::scal(double alpha){

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
 * @param rc the row/col index of the peps
 * @param peps the peps 'operator'
 * @return the expectation value of the double peps row within the boundary MPO
 */
template<>
void MPO<double>::normalize() {

   double nrm = sqrt(this->dot(*this));
   this->scal(1.0/nrm);

}

/**
 * canonicalize the mps
 * @param dir Left or Right canonicalization
 * @param norm if true: normalize, else not
 */
template<typename T>
void MPO<T>::canonicalize(const BTAS_SIDE &dir,bool norm){

   int length = this->size();

   if(dir == Left){//QR

      TArray<T,2> R;
      TArray<T,4> tmp;

      for(int i = 0;i < length - 1;++i){

         R.clear();

         //do QR
         Geqrf((*this)[i],R);

         //paste to next matrix
         tmp.clear();

         Contract((T)1.0,R,shape(1),(*this)[i + 1],shape(0),(T)0.0,tmp);

         (*this)[i + 1] = std::move(tmp);

      }

      T nrm = sqrt(Dotc((*this)[length-1],(*this)[length-1]));
      Scal(1.0/nrm,(*this)[length-1]);

      if(!norm){//redistribute over chain

         nrm = pow(nrm,1.0/(double)length);

         for(int i = 0;i < length;++i)
            Scal(nrm,(*this)[i]);

      }

   }
   else{//LQ

      TArray<T,2> L;
      TArray<T,4> tmp;

      for(int i = length - 1;i > 0;--i){

         L.clear();

         //do QR
         Gelqf(L,(*this)[i]);

         //paste to previous matrix
         tmp.clear();

         Contract((T)1.0,(*this)[i - 1],shape(3),L,shape(0),(T)0.0,tmp);

         (*this)[i - 1] = std::move(tmp);

      }

      T nrm = sqrt(Dotc((*this)[0],(*this)[0]));
      Scal(1.0/nrm,(*this)[0]);

      if(!norm){//redistribute over chain

         nrm = pow(nrm,1.0/(double)length);

         for(int i = 0;i < length;++i)
            Scal(nrm,(*this)[i]);

      }

   }

}

template MPO<double>::MPO();
template MPO< complex<double> >::MPO();

template MPO<double>::MPO(int);
template MPO< complex<double> >::MPO(int);

template MPO<double>::MPO(int,int,int);
template MPO< complex<double> >::MPO(int,int,int);

template MPO<double>::MPO(const MPO<double> &);
template MPO< complex<double> >::MPO(const MPO< complex<double> > &);

template MPO<double>::~MPO();
template MPO< complex<double> >::~MPO();

template int MPO<double>::gD() const;
template int MPO< complex<double> >::gD() const;

template double MPO<double>::dot(const MPO<double> &bra) const;
template  complex<double>  MPO< complex<double> >::dot(const MPO< complex<double> > &bra) const;

template double MPO<double>::expect(const char,int,const PEPS<double> &) const;
template  complex<double>  MPO< complex<double> >::expect(const char,int,const PEPS<complex< double> > &) const;

template void MPO<double>::fill_Random();

template void MPO<double>::canonicalize(const BTAS_SIDE &dir,bool norm);
template void MPO< complex<double> >::canonicalize(const BTAS_SIDE &dir,bool norm);
