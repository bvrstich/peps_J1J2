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

      void calc(const char,PEPS<double> &);

      void add_layer(const char,int,PEPS<double> &);

      double cost_function(const char,int,int,const PEPS<double> &,const std::vector< DArray<4> > &);

      void test();

      const MPO<double> &gt(int) const;
      MPO<double> &gt(int);

      const MPO<double> &gb(int) const;
      MPO<double> &gb(int);

      const vector< MPO<double> > &gt() const;
      const vector< MPO<double> > &gb() const;

      int gD() const;
      int gD_aux() const;

      int gcomp_sweeps() const;

      void sD(int);
      void sD_aux(int);

      void init_svd(char,int,const PEPS<double> &);

   private:

      //!stores an array environment MPO's for t(op) and b(ottom)
      vector< MPO<double> > t;
      vector< MPO<double> > b;

      //!regular bond dimension of peps
      int D;

      //!Auxiliary dimension, for the contraction
      int D_aux;

      //!nr of sweeps in compression
      int comp_sweeps;

};

#endif
