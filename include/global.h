#ifndef GLOBAL_H
#define GLOBAL_H

#include <iostream>
#include <fstream>
#include <vector>
#include <complex>

#include <btas/common/blas_cxx_interface.h>
#include <btas/common/TVector.h>
#include <btas/DENSE/TArray.h>

using std::ostream;
using std::vector;
using std::complex;

using namespace btas;

class Environment;

namespace global {

   //!random number generator
   extern Random RN;

   //!x dimension of the lattice, nr of cols
   extern int Lx;

   //!y dimension of the lattice, nr of rows
   extern int Ly;

   //!physical dimension of sites
   extern int d;

   //!bond dimension of the tensors
   extern int D;

   //!auxiliary bond dimension of the tensors
   extern int D_aux;

   //!next-nearest neigbour coupling
   extern double J2;

   //!hamiltonian object, containing nn-interaction operators and coefficients
   extern Hamiltonian ham;

   //!nr of sweeps for MPO compression
   extern int comp_sweeps;

   //!physical unit matrix
   extern DArray<2> I;

   //!Environment object, needed for basically everything when calculating
   extern Environment env;

   //!initializer
   void init(int,int,int,int,int,int,double);

   //!set the timestep
   void stau(double);
   
   void sD(int);

   template<typename T>
      T rgen();

};

#endif

/* vim: set ts=3 sw=3 expandtab :*/
