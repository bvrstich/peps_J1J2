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

/**
 * empty constructor
 */
Hamiltonian::Hamiltonian(){ }

/**
 * copy constructor
 * @param ham_c object to copy
 */
Hamiltonian::Hamiltonian(const Hamiltonian &ham_c){

   this->delta = ham_c.gdelta();

   this->L = ham_c.gL();
   this->R = ham_c.gR();

   this->coef_n = ham_c.gcoef_n();
   this->coef_nn = ham_c.gcoef_nn();

}

/**
 * empty destructor
 */
Hamiltonian::~Hamiltonian(){ }

/**
 * @return the number of terms in the nn-interaction
 */
int Hamiltonian::gdelta() const {

   return delta;

}

/**
 * @param i index of the nn term
 * @return coefficient 'i' to the next-nearest-neighbour interaction
 */
const double &Hamiltonian::gcoef_n(int i) const {

   return coef_n[i];

}

/**
 * @param i index of the nn term
 * @return coefficient 'i' to the next-nearest-neighbour interaction
 */
const double &Hamiltonian::gcoef_nn(int i) const {

   return coef_nn[i];

}

/**
 * @return array of 'delta' coefficients to the nearest-neighbour interaction
 */
const std::vector<double> &Hamiltonian::gcoef_n() const {

   return coef_n;

}

/**
 * @return array of 'delta' coefficients to the next-nearest-neighbour interaction
 */
const std::vector<double> &Hamiltonian::gcoef_nn() const {

   return coef_nn;

}

/**
 * @param i index of operator to return
 * @return the left operator with index 'i'
 */
const DArray<2> &Hamiltonian::gL(int i) const {

   return L[i];

}

/**
 * @param i index of operator to return
 * @return the right operator with index 'i'
 */
const DArray<2> &Hamiltonian::gR(int i) const {

   return R[i];

}

/**
 * @return the array of 'delta' left operators
 */
const std::vector< DArray<2> > &Hamiltonian::gL() const {

   return L;

}

/**
 * @return the array of 'delta' right operators
 */
const std::vector< DArray<2> > &Hamiltonian::gR() const {

   return R;

}

/**
 * initialize the operators on the nn-Heisenberg model
 * @param ladder if true, use ladder operators (+,-,z), if false, use x,y,z
 */
void Hamiltonian::set_J1J2(bool ladder) {

   delta = 3;

   L.resize(delta);
   R.resize(delta);

   coef_n.resize(delta);
   coef_nn.resize(delta);

   if(ladder){

      L[0].resize(global::d,global::d); //S+
      L[1].resize(global::d,global::d); //S-
      L[2].resize(global::d,global::d); //Sz

      R[0].resize(global::d,global::d); //S-
      R[1].resize(global::d,global::d); //S+
      R[2].resize(global::d,global::d); //Sz

      L[0] = 0.0; L[1] = 0.0; L[2] = 0.0;
      R[0] = 0.0; R[1] = 0.0; R[2] = 0.0;

      L[0](0,0) = 1.0;
      L[0](1,1) = 1.0;

      L[1](0,0) = 1.0;
      L[1](1,1) = 1.0;

      L[2](0,0) = 1.0;
      L[2](1,1) = 1.0;

      R[0](0,0) = 1.0;
      R[0](1,1) = 1.0;

      R[1](0,0) = 1.0;
      R[1](1,1) = 1.0;

      R[2](0,0) = 1.0;
      R[2](1,1) = 1.0;

      //nearest-neigbour coefficients:
      coef_n[0] = 1.0;//minus sign because of the Marshall sign rule
      coef_n[1] = 1.0;//minus sign because of the Marshall sign rule
      coef_n[2] = 1.0;
      
      //next-nearest neigbour coefficients:
      coef_nn[0] = 0.5 * global::J2;
      coef_nn[1] = 0.5 * global::J2;
      coef_nn[2] = 1.0 * global::J2;

/*
      //S+
      L[0](1,0) = 1.0;
      R[1](1,0) = 1.0;

      //S-
      L[1](0,1) = 1.0;
      R[0](0,1) = 1.0;

      //Sz
      L[2](0,0) = -0.5;
      L[2](1,1) = 0.5;

      R[2](0,0) = -0.5;
      R[2](1,1) = 0.5;

      //nearest-neigbour coefficients:
      coef_n[0] = -0.5;//add minus sign to use the Marshall sign rule
      coef_n[1] = -0.5;//add minus sign to use the Marshall sign rule
      coef_n[2] = 1.0;
      
      //next-nearest neigbour coefficients:
      coef_nn[0] = 0.5 * global::J2;
      coef_nn[1] = 0.5 * global::J2;
      coef_nn[2] = 1.0 * global::J2;
  */
   }
   else{

      L[0].resize(global::d,global::d); //Sx
      L[1].resize(global::d,global::d); //iSy
      L[2].resize(global::d,global::d); //Sz

      R[0].resize(global::d,global::d); //Sx
      R[1].resize(global::d,global::d); //iSy
      R[2].resize(global::d,global::d); //Sz

      L[0] = 0.0; L[1] = 0.0; L[2] = 0.0;
      R[0] = 0.0; R[1] = 0.0; R[2] = 0.0;

      //Sx
      L[0](0,1) = 0.5;
      L[0](1,0) = 0.5;

      R[0](0,1) = 0.5;
      R[0](1,0) = 0.5;

      //iSy
      L[1](0,1) = -0.5;
      L[1](1,0) = 0.5;

      R[1](0,1) = -0.5;
      R[1](1,0) = 0.5;

      //Sz
      L[2](0,0) = -0.5;
      L[2](1,1) = 0.5;

      R[2](0,0) = -0.5;
      R[2](1,1) = 0.5;

      //nearest neigbour coefficients:
      coef_n[0] = -1.0;//minus sign because of the Marshall sign rule
      coef_n[1] = 1.0;//minus sign of the Marshall sign rule cancelled because of squared 'i'
      coef_n[2] = 1.0;

      //next-nearest-neigbour coefficients:
      coef_nn[0] = global::J2;
      coef_nn[1] = -global::J2;//minus sign ofbecause of squared 'i'
      coef_nn[2] = global::J2;

   }

}
