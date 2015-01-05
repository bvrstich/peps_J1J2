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
 * construct on given shape: specified to index N=2
 * @param orig_in input shape of N-leg tensor to be permuted
 * @param reorder_in permutation reordering, all indices have to be different
 */
template<>
Perm<2>::Perm(const IVector<2> &orig_in,const IVector<2> &reorder_in){

   orig = orig_in;
   reorder = reorder_in;

   int dim = 1;

   for(int i = 0;i < 2;++i)
      dim *= orig_in[i];

   //now make the inverse list, tells where the old index is with respect to the new/permuted one
   IVector<2> inverse;
   
   for(int i = 0;i < 2;++i)
      inverse[reorder[i]] = i;

   list.resize(dim);

   IVector<2> index;

   //loop over the new array
   for(index[0] = 0;index[0] < orig[reorder[0]];++index[0])
      for(index[1] = 0;index[1] < orig[reorder[1]];++index[1])
         list[ index[0]*orig[reorder[1]] + index[1] ] = index[inverse[0]] * orig[1] + index[inverse[1]];

   perm_tensor.resize( shape(orig[reorder[0]],orig[reorder[1]]) );


}

/** 
 * construct on given shape: specified to index N=3
 * @param orig_in input shape of N-leg tensor to be permuted
 * @param reorder_in permutation reordering, all indices have to be different
 */
template<>
Perm<3>::Perm(const IVector<3> &orig_in,const IVector<3> &reorder_in){

   orig = orig_in;
   reorder = reorder_in;

   int dim = 1;

   for(int i = 0;i < 3;++i)
      dim *= orig_in[i];
   
   //now make the inverse list, tells where the old index is with respect to the new/permuted one
   IVector<3> inverse;
   
   for(int i = 0;i < 3;++i)
      inverse[reorder[i]] = i;

   list.resize(dim);

   IVector<3> index;

   //loop over the new array
   for(index[0] = 0;index[0] < orig[reorder[0]];++index[0])
      for(index[1] = 0;index[1] < orig[reorder[1]];++index[1])
         for(index[2] = 0;index[2] < orig[reorder[2]];++index[2])
            list[ index[0]* orig[reorder[1]] * orig[reorder[2]] + index[1] * orig[reorder[2]] + index[2] ] = index[inverse[0]] * orig[1] * orig[2] + index[inverse[1]] * orig[2] + index[inverse[2]];

   perm_tensor.resize( shape(orig[reorder[0]],orig[reorder[1]],orig[reorder[2]]) );

}

/** 
 * construct on given shape: specified to index N=4
 * @param orig_in input shape of N-leg tensor to be permuted
 * @param reorder_in permutation reordering, all indices have to be different
 */
template<>
Perm<4>::Perm(const IVector<4> &orig_in,const IVector<4> &reorder_in){

   orig = orig_in;
   reorder = reorder_in;

   int dim = 1;

   for(int i = 0;i < 4;++i)
      dim *= orig_in[i];
   
   //now make the inverse list, tells where the old index is with respect to the new/permuted one
   IVector<4> inverse;
   
   for(int i = 0;i < 4;++i)
      inverse[reorder[i]] = i;

   list.resize(dim);

   IVector<4> index;

   //loop over the new array
   for(index[0] = 0;index[0] < orig[reorder[0]];++index[0])
      for(index[1] = 0;index[1] < orig[reorder[1]];++index[1])
         for(index[2] = 0;index[2] < orig[reorder[2]];++index[2])
            for(index[3] = 0;index[3] < orig[reorder[3]];++index[3]){

               list[ index[0] * orig[reorder[1]] * orig[reorder[2]] * orig[reorder[3]] + index[1] * orig[reorder[2]] * orig[reorder[3]] + index[2] * orig[reorder[3]] + index[3] ] 
                  
                  = index[inverse[0]] * orig[1] * orig[2] * orig[3] + index[inverse[1]] * orig[2] * orig[3] + index[inverse[2]] * orig[3] + index[inverse[3]];

            }

   perm_tensor.resize( shape(orig[reorder[0]],orig[reorder[1]],orig[reorder[2]],orig[reorder[3]]) );

}

/** 
 * construct on given shape: specified to index N=5
 * @param orig_in input shape of N-leg tensor to be permuted
 * @param reorder_in permutation reordering, all indices have to be different
 */
template<>
Perm<5>::Perm(const IVector<5> &orig_in,const IVector<5> &reorder_in){

   orig = orig_in;
   reorder = reorder_in;

   int dim = 1;

   for(int i = 0;i < 5;++i)
      dim *= orig_in[i];
   
   //now make the inverse list, tells where the old index is with respect to the new/permuted one
   IVector<5> inverse;
   
   for(int i = 0;i < 5;++i)
      inverse[reorder[i]] = i;

   list.resize(dim);

   IVector<5> index;

   //loop over the new array
   for(index[0] = 0;index[0] < orig[reorder[0]];++index[0])
      for(index[1] = 0;index[1] < orig[reorder[1]];++index[1])
         for(index[2] = 0;index[2] < orig[reorder[2]];++index[2])
            for(index[3] = 0;index[3] < orig[reorder[3]];++index[3])
               for(index[4] = 0;index[4] < orig[reorder[4]];++index[4]){

                  list[ index[0] * orig[reorder[1]] * orig[reorder[2]] * orig[reorder[3]] * orig[reorder[4]] + index[1] * orig[reorder[2]] * orig[reorder[3]] * orig[reorder[4]] 
                     
                     + index[2] * orig[reorder[3]] * orig[reorder[4]] + index[3] * orig[reorder[4]] + index[4] ] 

                     = index[inverse[0]] * orig[1] * orig[2] * orig[3] * orig[4] + index[inverse[1]] * orig[2] * orig[3] * orig[4] + index[inverse[2]] * orig[3] * orig[4] + index[inverse[3]] * orig[4] + index[inverse[4]];

               }

   perm_tensor.resize( shape(orig[reorder[0]],orig[reorder[1]],orig[reorder[2]],orig[reorder[3]],orig[reorder[4]]) );

}

/** 
 * construct on given shape: specified to index N=6
 * @param orig_in input shape of N-leg tensor to be permuted
 * @param reorder_in permutation reordering, all indices have to be different
 */
template<>
Perm<6>::Perm(const IVector<6> &orig_in,const IVector<6> &reorder_in){

   orig = orig_in;
   reorder = reorder_in;

   int dim = 1;

   for(int i = 0;i < 6;++i)
      dim *= orig_in[i];

   //now make the inverse list, tells where the old index is with respect to the new/permuted one
   IVector<6> inverse;

   for(int i = 0;i < 6;++i)
      inverse[reorder[i]] = i;

   list.resize(dim);

   IVector<6> index;

   //loop over the new array
   for(index[0] = 0;index[0] < orig[reorder[0]];++index[0])
      for(index[1] = 0;index[1] < orig[reorder[1]];++index[1])
         for(index[2] = 0;index[2] < orig[reorder[2]];++index[2])
            for(index[3] = 0;index[3] < orig[reorder[3]];++index[3])
               for(index[4] = 0;index[4] < orig[reorder[4]];++index[4])
                  for(index[5] = 0;index[5] < orig[reorder[5]];++index[5]){

                     list[ index[0] * orig[reorder[1]] * orig[reorder[2]] * orig[reorder[3]] * orig[reorder[4]] * orig[reorder[5]] + index[1] * orig[reorder[2]] * orig[reorder[3]] * orig[reorder[4]] * orig[reorder[5]] 

                        + index[2] * orig[reorder[3]] * orig[reorder[4]] * orig[reorder[5]] + index[3] * orig[reorder[4]] * orig[reorder[5]] + index[4] * orig[reorder[5]] + index[5] ] 

                        = index[inverse[0]] * orig[1] * orig[2] * orig[3] * orig[4] * orig[5] + index[inverse[1]] * orig[2] * orig[3] * orig[4] * orig[5]
                        
                        + index[inverse[2]] * orig[3] * orig[4] * orig[5] + index[inverse[3]] * orig[4] * orig[5] + index[inverse[4]] * orig[5] + index[inverse[5]];

                  }

   perm_tensor.resize( shape(orig[reorder[0]],orig[reorder[1]],orig[reorder[2]],orig[reorder[3]],orig[reorder[4]],orig[reorder[5]]) );

}

/** 
 * construct on given shape: specified to index N=7
 * @param orig_in input shape of N-leg tensor to be permuted
 * @param reorder_in permutation reordering, all indices have to be different
 */
template<>
Perm<7>::Perm(const IVector<7> &orig_in,const IVector<7> &reorder_in){

   orig = orig_in;
   reorder = reorder_in;

   int dim = 1;

   for(int i = 0;i < 7;++i)
      dim *= orig_in[i];

   //now make the inverse list, tells where the old index is with respect to the new/permuted one
   IVector<7> inverse;

   for(int i = 0;i < 7;++i)
      inverse[reorder[i]] = i;

   list.resize(dim);

   IVector<7> index;

   //loop over the new array
   for(index[0] = 0;index[0] < orig[reorder[0]];++index[0])
      for(index[1] = 0;index[1] < orig[reorder[1]];++index[1])
         for(index[2] = 0;index[2] < orig[reorder[2]];++index[2])
            for(index[3] = 0;index[3] < orig[reorder[3]];++index[3])
               for(index[4] = 0;index[4] < orig[reorder[4]];++index[4])
                  for(index[5] = 0;index[5] < orig[reorder[5]];++index[5])
                     for(index[6] = 0;index[6] < orig[reorder[6]];++index[6]){

                        list[ index[0] * orig[reorder[1]] * orig[reorder[2]] * orig[reorder[3]] * orig[reorder[4]] * orig[reorder[5]] * orig[reorder[6]] 
                           
                           + index[1] * orig[reorder[2]] * orig[reorder[3]] * orig[reorder[4]] * orig[reorder[5]] * orig[reorder[6]]

                           + index[2] * orig[reorder[3]] * orig[reorder[4]] * orig[reorder[5]] * orig[reorder[6]] + index[3] * orig[reorder[4]] * orig[reorder[5]] * orig[reorder[6]] + index[4] * orig[reorder[5]] * orig[reorder[6]]
                           
                           + index[5] * orig[reorder[6]] + index[6] ] 

                           = index[inverse[0]] * orig[1] * orig[2] * orig[3] * orig[4] * orig[5] * orig[6] + index[inverse[1]] * orig[2] * orig[3] * orig[4] * orig[5] * orig[6]

                           + index[inverse[2]] * orig[3] * orig[4] * orig[5] * orig[6] + index[inverse[3]] * orig[4] * orig[5] * orig[6] + index[inverse[4]] * orig[5] * orig[6] + index[inverse[5]] * orig[6] + index[inverse[6]];

                     }

   perm_tensor.resize( shape(orig[reorder[0]],orig[reorder[1]],orig[reorder[2]],orig[reorder[3]],orig[reorder[4]],orig[reorder[5]],orig[reorder[6]]) );

}

/** 
 * construct on given shape: specified to index N=8
 * @param orig_in input shape of N-leg tensor to be permuted
 * @param reorder_in permutation reordering, all indices have to be different
 */
template<>
Perm<8>::Perm(const IVector<8> &orig_in,const IVector<8> &reorder_in){

   orig = orig_in;
   reorder = reorder_in;

   int dim = 1;

   for(int i = 0;i < 8;++i)
      dim *= orig_in[i];

   //now make the inverse list, tells where the old index is with respect to the new/permuted one
   IVector<8> inverse;

   for(int i = 0;i < 8;++i)
      inverse[reorder[i]] = i;

   list.resize(dim);

   IVector<8> index;

   //loop over the new array
   for(index[0] = 0;index[0] < orig[reorder[0]];++index[0])
      for(index[1] = 0;index[1] < orig[reorder[1]];++index[1])
         for(index[2] = 0;index[2] < orig[reorder[2]];++index[2])
            for(index[3] = 0;index[3] < orig[reorder[3]];++index[3])
               for(index[4] = 0;index[4] < orig[reorder[4]];++index[4])
                  for(index[5] = 0;index[5] < orig[reorder[5]];++index[5])
                     for(index[6] = 0;index[6] < orig[reorder[6]];++index[6])
                        for(index[7] = 0;index[7] < orig[reorder[7]];++index[7]){

                           list[ index[0] * orig[reorder[1]] * orig[reorder[2]] * orig[reorder[3]] * orig[reorder[4]] * orig[reorder[5]] * orig[reorder[6]] * orig[reorder[7]] 

                              + index[1] * orig[reorder[2]] * orig[reorder[3]] * orig[reorder[4]] * orig[reorder[5]] * orig[reorder[6]] * orig[reorder[7]]

                              + index[2] * orig[reorder[3]] * orig[reorder[4]] * orig[reorder[5]] * orig[reorder[6]] * orig[reorder[7]] + index[3] * orig[reorder[4]] * orig[reorder[5]] * orig[reorder[6]] * orig[reorder[7]]
                              
                              + index[4] * orig[reorder[5]] * orig[reorder[6]]  * orig[reorder[7]] + index[5] * orig[reorder[6]] * orig[reorder[7]] + index[6] * orig[reorder[7]] + index[7] ] 

                              = index[inverse[0]] * orig[1] * orig[2] * orig[3] * orig[4] * orig[5] * orig[6] * orig[7] + index[inverse[1]] * orig[2] * orig[3] * orig[4] * orig[5] * orig[6] * orig[7]

                              + index[inverse[2]] * orig[3] * orig[4] * orig[5] * orig[6] * orig[7] + index[inverse[3]] * orig[4] * orig[5] * orig[6] * orig[7] + index[inverse[4]] * orig[5] * orig[6] * orig[7]
                              
                              + index[inverse[5]] * orig[6] * orig[7] + index[inverse[6]] * orig[7] + index[inverse[7]];

                        }

   perm_tensor.resize( shape(orig[reorder[0]],orig[reorder[1]],orig[reorder[2]],orig[reorder[3]],orig[reorder[4]],orig[reorder[5]],orig[reorder[6]],orig[reorder[7]]) );

}

/** 
 * construct on given shape: specified to index N=9
 * @param orig_in input shape of N-leg tensor to be permuted
 * @param reorder_in permutation reordering, all indices have to be different
 */
template<>
Perm<9>::Perm(const IVector<9> &orig_in,const IVector<9> &reorder_in){

   orig = orig_in;
   reorder = reorder_in;

   int dim = 1;

   for(int i = 0;i < 9;++i)
      dim *= orig_in[i];

   //now make the inverse list, tells where the old index is with respect to the new/permuted one
   IVector<9> inverse;

   for(int i = 0;i < 9;++i)
      inverse[reorder[i]] = i;

   list.resize(dim);

   IVector<9> index;

   //loop over the new array
   for(index[0] = 0;index[0] < orig[reorder[0]];++index[0])
      for(index[1] = 0;index[1] < orig[reorder[1]];++index[1])
         for(index[2] = 0;index[2] < orig[reorder[2]];++index[2])
            for(index[3] = 0;index[3] < orig[reorder[3]];++index[3])
               for(index[4] = 0;index[4] < orig[reorder[4]];++index[4])
                  for(index[5] = 0;index[5] < orig[reorder[5]];++index[5])
                     for(index[6] = 0;index[6] < orig[reorder[6]];++index[6])
                        for(index[7] = 0;index[7] < orig[reorder[7]];++index[7])
                           for(index[8] = 0;index[8] < orig[reorder[8]];++index[8]){

                              list[ index[0] * orig[reorder[1]] * orig[reorder[2]] * orig[reorder[3]] * orig[reorder[4]] * orig[reorder[5]] * orig[reorder[6]] * orig[reorder[7]] * orig[reorder[8]] 

                                 + index[1] * orig[reorder[2]] * orig[reorder[3]] * orig[reorder[4]] * orig[reorder[5]] * orig[reorder[6]] * orig[reorder[7]] * orig[reorder[8]]

                                 + index[2] * orig[reorder[3]] * orig[reorder[4]] * orig[reorder[5]] * orig[reorder[6]] * orig[reorder[7]] * orig[reorder[8]] 
                                 
                                 + index[3] * orig[reorder[4]] * orig[reorder[5]] * orig[reorder[6]] * orig[reorder[7]] * orig[reorder[8]] + index[4] * orig[reorder[5]] * orig[reorder[6]]  * orig[reorder[7]] * orig[reorder[8]]
                                 
                                 + index[5] * orig[reorder[6]] * orig[reorder[7]] * orig[reorder[8]] + index[6] * orig[reorder[7]] * orig[reorder[8]] + index[7] * orig[reorder[8]] + index[8] ] 

                                 = index[inverse[0]] * orig[1] * orig[2] * orig[3] * orig[4] * orig[5] * orig[6] * orig[7] * orig[8] + index[inverse[1]] * orig[2] * orig[3] * orig[4] * orig[5] * orig[6] * orig[7] * orig[8]

                                 + index[inverse[2]] * orig[3] * orig[4] * orig[5] * orig[6] * orig[7] * orig[8] + index[inverse[3]] * orig[4] * orig[5] * orig[6] * orig[7] * orig[8] 
                                 
                                 + index[inverse[4]] * orig[5] * orig[6] * orig[7] * orig[8] + index[inverse[5]] * orig[6] * orig[7] * orig[8] + index[inverse[6]] * orig[7] * orig[8] + index[inverse[7]] * orig[8] + index[inverse[8]];

                           }

   perm_tensor.resize( shape(orig[reorder[0]],orig[reorder[1]],orig[reorder[2]],orig[reorder[3]],orig[reorder[4]],orig[reorder[5]],orig[reorder[6]],orig[reorder[7]],orig[reorder[8]]) );

}

/**
 * copy constructor
 * @param perm_copy object to be copied
 */
template<size_t N>
Perm<N>::Perm(const Perm<N> &perm_copy){ 

   orig = perm_copy.gorig();
   reorder = perm_copy.greorder();
   perm_tensor = perm_copy.gperm_tensor();

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
 * @return the actual permuted tensor
 */
template<size_t N>
const DArray<N> &Perm<N>::gperm_tensor() const {

   return perm_tensor;

}

/**
 * finally: permute! just a single loop!
 * @param orig_tensor original tensor, const
 */
template<size_t N>
void Perm<N>::permute(const DArray<N> &orig_tensor){

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

//copy constructor
template Perm<2>::Perm(const Perm<2> &perm_copy);
template Perm<3>::Perm(const Perm<3> &perm_copy);
template Perm<4>::Perm(const Perm<4> &perm_copy);
template Perm<5>::Perm(const Perm<5> &perm_copy);
template Perm<6>::Perm(const Perm<6> &perm_copy);
template Perm<7>::Perm(const Perm<7> &perm_copy);
template Perm<8>::Perm(const Perm<8> &perm_copy);
template Perm<9>::Perm(const Perm<9> &perm_copy);

//destructor
template Perm<2>::~Perm();
template Perm<3>::~Perm();
template Perm<4>::~Perm();
template Perm<5>::~Perm();
template Perm<6>::~Perm();
template Perm<7>::~Perm();
template Perm<8>::~Perm();
template Perm<9>::~Perm();

//getter for original shape
template const IVector<2> &Perm<2>::gorig() const;
template const IVector<3> &Perm<3>::gorig() const;
template const IVector<4> &Perm<4>::gorig() const;
template const IVector<5> &Perm<5>::gorig() const;
template const IVector<6> &Perm<6>::gorig() const;
template const IVector<7> &Perm<7>::gorig() const;
template const IVector<8> &Perm<8>::gorig() const;
template const IVector<9> &Perm<9>::gorig() const;

//getter for the reordering vector
template const IVector<2> &Perm<2>::greorder() const;
template const IVector<3> &Perm<3>::greorder() const;
template const IVector<4> &Perm<4>::greorder() const;
template const IVector<5> &Perm<5>::greorder() const;
template const IVector<6> &Perm<6>::greorder() const;
template const IVector<7> &Perm<7>::greorder() const;
template const IVector<8> &Perm<8>::greorder() const;
template const IVector<9> &Perm<9>::greorder() const;

//getter for the permuted tensor
template const DArray<2> &Perm<2>::gperm_tensor() const;
template const DArray<3> &Perm<3>::gperm_tensor() const;
template const DArray<4> &Perm<4>::gperm_tensor() const;
template const DArray<5> &Perm<5>::gperm_tensor() const;
template const DArray<6> &Perm<6>::gperm_tensor() const;
template const DArray<7> &Perm<7>::gperm_tensor() const;
template const DArray<8> &Perm<8>::gperm_tensor() const;
template const DArray<9> &Perm<9>::gperm_tensor() const;

//actual permute function
template void Perm<2>::permute(const DArray<2> &orig_tensor);
template void Perm<3>::permute(const DArray<3> &orig_tensor);
template void Perm<4>::permute(const DArray<4> &orig_tensor);
template void Perm<5>::permute(const DArray<5> &orig_tensor);
template void Perm<6>::permute(const DArray<6> &orig_tensor);
template void Perm<7>::permute(const DArray<7> &orig_tensor);
template void Perm<8>::permute(const DArray<8> &orig_tensor);
template void Perm<9>::permute(const DArray<9> &orig_tensor);
