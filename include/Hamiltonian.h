#ifndef HAMILTONIAN_H
#define HAMILTONIAN_H

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
 * @date 24-10-2014\n\n
 * This class contains the nearest neighbour interaction operators which the Hamiltonian consists of
 */
class Hamiltonian {

   public:
      
      //empty constructor
      Hamiltonian();

      //copy constructur
      Hamiltonian(const Hamiltonian &);

      virtual ~Hamiltonian();

      void set_J1J2(bool);

      int gdelta() const;

      const DArray<2> &gL(int) const;

      const DArray<2> &gR(int) const;

      const std::vector< DArray<2> > &gL() const;

      const std::vector< DArray<2> > &gR() const;

      const double &gcoef_n(int) const;

      const std::vector<double> &gcoef_n() const;

      const double &gcoef_nn(int) const;

      const std::vector<double> &gcoef_nn() const;

   private:

      //!number of terms in the hamiltonian (i.e. 3 in Heisenberg, 2 in XY model, 1 for ising,...)
      int delta;

      //!coefficients to the nearest neighbour interaction
      std::vector<double> coef_n;

      //!coefficients to the next-nearest neighbour interaction
      std::vector<double> coef_nn;

      //!Left and right operators (array of 'delta' operators on the left side of a spin-pair)
      std::vector< DArray<2> > L;
      std::vector< DArray<2> > R;

};

#endif

/* vim: set ts=3 sw=3 expandtab :*/
