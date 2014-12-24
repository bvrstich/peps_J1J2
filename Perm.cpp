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
template<size_t N>
Perm<N>::Perm(){ }

/** 
 * construct on given shape
 * @param orig_in input shape of N-leg tensor to be permuted
 * @param reorder_in permutation reordering, all indices have to be different
 */
template<size_t N>
Perm<N>::Perm(const IVector<N> &orig_in,const IVector<N> &reorder_in){

   orig = orig_in;
   reorder = reorder_in;

   int dim = 1;

   for(int i = 0;i < N;++i){

      dim *= orig_in[i];

   }

   list.resize(dim);

}

/**
 * copy constructor
 * @param perm_copy object to be copied
 */
template<size_t N>
Perm<N>::Perm(const Perm<N> &perm_copy){ 

   orig = perm_copy.gorig();
   reorder = perm_copy.greorder();

}

/**
 * pointless destructor
 */
template<size_t N>
Perm<N>::~Perm(){ }

/**
 * @return the shape vector for the permutation
 */
template<size_t N>
const IVector<N> &Perm<N>::gorig() const {

   return orig;

}

/**
 * @return the shape vector for the permutation
 */
template<size_t N>
const IVector<N> &Perm<N>::greorder() const {

   return reorder;

}

/**
 * finally: permute!
 * @param orig_tensor original tensor, const
 * @param perm_tensor on exit, contains permuted tensor
 */
template<>
void Perm<2>::permute(const DArray<2> &orig_tensor,DArray<2> &perm_tensor){

   IVector<2> index;

   //loop over the new array
   for(index[0] = 0;index[0] < orig[reorder[0]];++index[0])
      for(index[1] = 0;index[1] < orig[reorder[1]];++index[1])
         list[ index[0]*orig[reorder[1]] + index[1] ] = index[reorder[0]] * orig[1] + index[reorder[1]];

   perm_tensor.resize( shape(orig[reorder[0]],orig[reorder[1]]) );

   for(int i = 0;i < list.size();++i)
      perm_tensor.data()[i] = orig_tensor.data()[list[i]];

}

//specification for the dimensions to be used (up to 10?)... god why

//empty constructor
template Perm<2>::Perm();
template Perm<3>::Perm();
template Perm<4>::Perm();
template Perm<5>::Perm();
template Perm<6>::Perm();
template Perm<7>::Perm();
template Perm<8>::Perm();
template Perm<9>::Perm();
template Perm<10>::Perm();

//constructor taking in shape of tensor
template Perm<2>::Perm(const IVector<2> &orig_in,const IVector<2> &reorder_in);
template Perm<3>::Perm(const IVector<3> &orig_in,const IVector<3> &reorder_in);
template Perm<4>::Perm(const IVector<4> &orig_in,const IVector<4> &reorder_in);
template Perm<5>::Perm(const IVector<5> &orig_in,const IVector<5> &reorder_in);
template Perm<6>::Perm(const IVector<6> &orig_in,const IVector<6> &reorder_in);
template Perm<7>::Perm(const IVector<7> &orig_in,const IVector<7> &reorder_in);
template Perm<8>::Perm(const IVector<8> &orig_in,const IVector<8> &reorder_in);
template Perm<9>::Perm(const IVector<9> &orig_in,const IVector<9> &reorder_in);
template Perm<10>::Perm(const IVector<10> &orig_in,const IVector<10> &reorder_in);

//constructor taking in shape of tensor
template Perm<2>::Perm(const Perm<2> &perm_copy);
template Perm<3>::Perm(const Perm<3> &perm_copy);
template Perm<4>::Perm(const Perm<4> &perm_copy);
template Perm<5>::Perm(const Perm<5> &perm_copy);
template Perm<6>::Perm(const Perm<6> &perm_copy);
template Perm<7>::Perm(const Perm<7> &perm_copy);
template Perm<8>::Perm(const Perm<8> &perm_copy);
template Perm<9>::Perm(const Perm<9> &perm_copy);
template Perm<10>::Perm(const Perm<10> &perm_copy);

//destructor
template Perm<2>::~Perm();
template Perm<3>::~Perm();
template Perm<4>::~Perm();
template Perm<5>::~Perm();
template Perm<6>::~Perm();
template Perm<7>::~Perm();
template Perm<8>::~Perm();
template Perm<9>::~Perm();
template Perm<10>::~Perm();

//getter for original shape
template const IVector<2> &Perm<2>::gorig() const;
template const IVector<3> &Perm<3>::gorig() const;
template const IVector<4> &Perm<4>::gorig() const;
template const IVector<5> &Perm<5>::gorig() const;
template const IVector<6> &Perm<6>::gorig() const;
template const IVector<7> &Perm<7>::gorig() const;
template const IVector<8> &Perm<8>::gorig() const;
template const IVector<9> &Perm<9>::gorig() const;
template const IVector<10> &Perm<10>::gorig() const;

template const IVector<2> &Perm<2>::greorder() const;
template const IVector<3> &Perm<3>::greorder() const;
template const IVector<4> &Perm<4>::greorder() const;
template const IVector<5> &Perm<5>::greorder() const;
template const IVector<6> &Perm<6>::greorder() const;
template const IVector<7> &Perm<7>::greorder() const;
template const IVector<8> &Perm<8>::greorder() const;
template const IVector<9> &Perm<9>::greorder() const;
template const IVector<10> &Perm<10>::greorder() const;
