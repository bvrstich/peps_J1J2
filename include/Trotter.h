#ifndef TROTTER_H
#define TROTTER_H

#include <iostream>
#include <fstream>
#include <vector>

#include <btas/common/blas_cxx_interface.h>
#include <btas/common/TVector.h>
#include <btas/DENSE/TArray.h>

using std::ostream;
using std::vector;

using namespace btas;

/**
 * @author Brecht Verstichel
 * @date 13-05-2014\n\n
 * This class Trotter contains the two-site gates for the imaginary time evolution in the trotter decomposition
 */
class Trotter {

   public:

      Trotter();

      Trotter(double tau);

      Trotter(const Trotter &);

      virtual ~Trotter();

      double gtau() const;

      const DArray<3> &gLO() const;

      const DArray<3> &gRO() const;

      const DArray<2> &geB() const;

   private:
      
      //!actual operators: left
      DArray<3> LO;

      //!actual operators: right
      DArray<3> RO;

      //!local propagator
      DArray<2> eB;

      //!timestep
      double tau;


};

#endif

/* vim: set ts=3 sw=3 expandtab :*/
