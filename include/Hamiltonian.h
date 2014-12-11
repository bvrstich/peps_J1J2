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

      void set_heisenberg(bool);

      void set_transverse_field_ising(double);

      int gdelta() const;

      bool gis_local() const;

      const DArray<2> &gL(int) const;

      const DArray<2> &gR(int) const;

      const std::vector< DArray<2> > &gL() const;

      const std::vector< DArray<2> > &gR() const;

      const DArray<2> &gB() const;

      const double &gcoef(int) const;

      const std::vector<double> &gcoef() const;

   private:

      //!number of terms in the hamiltonian (i.e. 3 in Heisenberg, 2 in XY model, 1 for ising,...)
      int delta;

      //!coefficients to the nearest neighbour interaction
      std::vector<double> coef;

      //!Left operators (array of 'delta' operators on the left side of the nn-pair)
      std::vector< DArray<2> > L;

      //!Right operators (array of 'delta' operators on the right side of the nn-pair)
      std::vector< DArray<2> > R;

      //!local (on-site) operator
      DArray<2> B;

      //!flag if there is local term or not
      bool is_local;

};

#endif

/* vim: set ts=3 sw=3 expandtab :*/
