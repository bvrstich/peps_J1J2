#ifndef PERM_H
#define PERM_H

#include <iostream>
#include <fstream>
#include <vector>

#include <btas/common/blas_cxx_interface.h>
#include <btas/common/TVector.h>
#include <btas/DENSE/TArray.h>

using std::ostream;
using std::vector;

/**
 * @author Brecht Verstichel
 * @data 24-12-2014\n\n
 * class used for permutations of tensors
 */
template<size_t N>
class Perm {

   public:

      Perm();

      Perm(const IVector<N> &,const IVector<N> &);

      Perm(const Perm &);

      virtual ~Perm();

      const IVector<N> &gorig() const;

      const IVector<N> &greorder() const;

      void permute(const DArray<N> &);

      const std::vector<int> &glist() const;

      const DArray<N> &gperm_tensor() const;

   private:

      //!shape of the original tensor
      IVector<N> orig;

      //!reordering
      IVector<N> reorder;

      //!list mapping old tensor to new
      std::vector<int> list;

      //!new permuted tensor, allocated at constructor level
      DArray<N> perm_tensor;

};

#endif
