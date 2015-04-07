#ifndef MPO_H
#define MPO_H

#include <iostream>
#include <fstream>
#include <vector>

#include <btas/common/blas_cxx_interface.h>
#include <btas/common/TVector.h>
#include <btas/DENSE/TArray.h>

using std::ostream;
using std::vector;

using namespace btas;

#include "PEPS.h"

/**
 * @author Brecht Verstichel
 * @date 26-03-2014\n\n
 * This class MPO is a class written for the construction of matrix product states without symmetry.
 * More specifically it will be used for the contraction of PEPS networks. Where the reduction to a MPO-like for is done.
 */
template<typename T>
class MPO : public vector< TArray<T,4> > {

   public:

      MPO();

      MPO(int);

      MPO(int,int,int);

      //copy constructor
      MPO(const MPO &);

      //destructor
      virtual ~MPO();

      T dot(const MPO<T> &bra) const;

      T expect(const char,int,const PEPS<T> &) const;

      int gD() const;

      int gd_phys() const;

      void fill(const char,const PEPS<double> &);

      void fill_Random();

      void scal(T );

      void canonicalize(const BTAS_SIDE &,bool);

      void normalize();

   private:

      //!dimension of the bonds
      int D;

      //!physical dimension
      int d_phys;


};

/**
 * output stream operator overloaded for MPO<T> 
 */
template<typename T>
ostream &operator<<(ostream &output,const MPO<T> &mpo_p){

   for(int s = 0;s < mpo_p.size();++s){

         output << std::endl;
         output << "Tensor on site (" << s << ")\t" << std::endl;
         output << std::endl;
         output << mpo_p[s] << std::endl;

      }

   return output;

}

#endif

/* vim: set ts=3 sw=3 expandtab :*/
