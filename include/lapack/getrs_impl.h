#ifndef __BTAS_LAPACK_GETRS_IMPL_H
#define __BTAS_LAPACK_GETRS_IMPL_H 1

#include <lapack/types.h>

namespace btas {
   namespace lapack {

      template<typename T>
         void getrs (
               const int& order,
               const char &trans,
               const size_t& N,
               const size_t &nrhs,
               T* A,
               const size_t& ldA,
               int* ipiv,
               T* b,
               const size_t &ldb)
         {
            BTAS_LAPACK_ASSERT(false, "getrs must be specialized.");
         }

      inline void getrs (
            const int& order,
            const char &trans,
            const size_t& N,
            const size_t &nrhs,
            float* A,
            const size_t& ldA,
            int* ipiv,
            float* b,
            const size_t &ldb)
      {
         LAPACKE_sgetrs(order, trans, N, nrhs, A, ldA, ipiv,b,ldb);
      }

      inline void getrs (
            const int& order,
            const char &trans,
            const size_t& N,
            const size_t &nrhs,
            double* A,
            const size_t& ldA,
            int* ipiv,
            double* b,
            const size_t &ldb)
      {
         LAPACKE_dgetrs(order, trans, N, nrhs, A, ldA, ipiv,b,ldb);

      }

      inline void getrs (
            const int& order,
            const char &trans,
            const size_t& N,
            const size_t &nrhs,
            std::complex<float>* A,
            const size_t& ldA,
            int* ipiv,
            std::complex<float>* b,
            const size_t &ldb)
      {
         LAPACKE_cgetrs(order, trans, N, nrhs, A, ldA, ipiv,b,ldb);
      }

      inline void getrs (
            const int& order,
            const char &trans,
            const size_t& N,
            const size_t &nrhs,
            std::complex<double>* A,
            const size_t& ldA,
            int* ipiv,
            std::complex<double>* b,
            const size_t &ldb)
      {
         LAPACKE_zgetrs(order, trans, N, nrhs, A, ldA, ipiv,b,ldb);
      }

   } // namespace lapack
} // namespace btas

#endif // __BTAS_LAPACK_GETRS_IMPL_H
