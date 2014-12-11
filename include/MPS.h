#ifndef MPS_H
#define MPS_H

#include <iostream>
#include <fstream>
#include <vector>

#include <btas/common/blas_cxx_interface.h>
#include <btas/common/TVector.h>
#include <btas/DENSE/TArray.h>

using std::ostream;
using std::vector;

using namespace btas;

#include "MPO.h"

/**
 * @author Brecht Verstichel
 * @date 26-03-2014\n\n
 * This class MPS is a class written for the construction of matrix product states without symmetry.
 * More specifically it will be used for the contraction of PEPS networks. Where the reduction to a MPS-like form is done.
 */
template<typename T>
class MPS : public vector< TArray<T,3> > {

   public:

      MPS();

      MPS(int L);

      MPS(int,int,int);

      //copy constructor
      MPS(const MPS &);

      //destructor
      virtual ~MPS();

      int gD() const;

      int gd_phys() const;

      void gemv(char , const MPO<T> &);

      void canonicalize(const BTAS_SIDE &,bool);

      void scal(T );

      T dot(const MPS<T> &bra) const;

      T normalize();

   private:

      //!dimension of the bonds
      int D;

      //!physical dimension
      int d_phys;
};

/**
 * output stream operator overloaded for MPS<T> 
 */
template<typename T>
ostream &operator<<(ostream &output,const MPS<T> &mps_p){

   for(int s = 0;s < mps_p.size();++s){

         output << std::endl;
         output << "Tensor on site (" << s << ")\t" << std::endl;
         output << std::endl;
         output << mps_p[s] << std::endl;

      }

   return output;

}

#endif

/* vim: set ts=3 sw=3 expandtab :*/
