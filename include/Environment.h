#ifndef ENVIRONMENT_H
#define ENVIRONMENT_H

#include <iostream>
#include <fstream>
#include <vector>

#include <btas/common/blas_cxx_interface.h>
#include <btas/common/TVector.h>
#include <btas/DENSE/TArray.h>

using std::ostream;
using std::vector;

template<typename T>
class PEPS;

template<typename T>
class MPO;

/**
 * @author Brecht Verstichel
 * @data 02-05-2014\n\n
 * Class used to calculate the enviroment of a peps. Needed for the calculation of expectation values and the update of tensors.
 */
class Environment {

   public:

      Environment();

      Environment(int,int,int);

      //copy constructor
      Environment(const Environment &);

      //destructor
      virtual ~Environment();

      void calc(const char,const PEPS<double> &);

      void add_layer(const char,int,const PEPS<double> &);

      void test();

      const MPO<double> &gl(int) const;
      MPO<double> &gl(int);

      const MPO<double> &gr(int) const;
      MPO<double> &gr(int);

      const MPO<double> &gt(int) const;
      MPO<double> &gt(int);

      const MPO<double> &gb(int) const;
      MPO<double> &gb(int);

      const vector< MPO<double> > &gl() const;
      const vector< MPO<double> > &gr() const;
      const vector< MPO<double> > &gt() const;
      const vector< MPO<double> > &gb() const;

      int gD() const;
      int gD_aux() const;

      int gcomp_sweeps() const;

      void sD(int);
      void sD_aux(int);

   private:

      //!stores an array environment MPO's for l(eft) , r(ight), t(op) and b(ottom)
      vector< MPO<double> > l;
      vector< MPO<double> > r;
      vector< MPO<double> > t;
      vector< MPO<double> > b;

      //!regular bond dimension of peps
      int D;

      //!Auxiliary dimension, for the contraction
      int D_aux;

      //!flags that tell if previous guess for environment is present
      bool flag_l;
      bool flag_r;
      bool flag_t;
      bool flag_b;

      //!nr of sweeps in compression
      int comp_sweeps;

};

#endif
