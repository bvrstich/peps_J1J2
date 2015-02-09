#ifndef __BTAS_LAPACK_SYTRS_IMPL_H
#define __BTAS_LAPACK_SYTRS_IMPL_H 1

#include <lapack/types.h>

namespace btas {
   namespace lapack {

      template<typename T>
         void sytrs (
               const int& order,
               const char &uplo,
               const size_t& N,
               const size_t& nrs,
               const T* A,
               const size_t& ldA,
               const int *ipiv,
               T* b,
               const size_t& ldb
               )
         {
            BTAS_LAPACK_ASSERT(false, "sytrs must be specialized.");
         }

      inline void sytrs (
            const int& order,
            const char &uplo,
            const size_t& N,
            const size_t& nrs,
            const float* A,
            const size_t& ldA,
            const int *ipiv,
            float* b,
            const size_t& ldb)
      {
         LAPACKE_ssytrs(order,uplo, N,nrs, A, ldA,ipiv,b,ldb);
      }

      inline void sytrs (
            const int& order,
            const char &uplo,
            const size_t& N,
            const size_t& nrs,
            const double* A,
            const size_t& ldA,
            const int *ipiv,
            double* b,
            const size_t& ldb)
      {
         LAPACKE_dsytrs(order,uplo, N,nrs, A, ldA,ipiv,b,ldb);

      }

   } // namespace lapack
} // namespace btas

#endif // __BTAS_LAPACK_SYTRS_IMPL_H
