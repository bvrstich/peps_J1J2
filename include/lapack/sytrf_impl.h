#ifndef __BTAS_LAPACK_SYTRF_IMPL_H
#define __BTAS_LAPACK_SYTRF_IMPL_H 1

#include <lapack/types.h>

namespace btas {
   namespace lapack {

      template<typename T>
         void sytrf (
               const int& order,
               const char &uplo,
               const size_t& N,
               T* A,
               const size_t& ldA,
               int *ipiv
               )
         {
            BTAS_LAPACK_ASSERT(false, "sytrf must be specialized.");
         }

      inline void sytrf (
            const int& order,
            const char &uplo,
            const size_t& N,
            float* A,
            const size_t& ldA,
            int *ipiv
            )
      {
         LAPACKE_ssytrf(order,uplo, N, A, ldA,ipiv);
      }

      inline void sytrf (
            const int& order,
            const char &uplo,
            const size_t& N,
            double* A,
            const size_t& ldA,
            int *ipiv
            )
      {
         LAPACKE_dsytrf(order,uplo, N, A, ldA,ipiv);

      }

   } // namespace lapack
} // namespace btas

#endif // __BTAS_LAPACK_SYTRF_IMPL_H
