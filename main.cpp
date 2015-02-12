/**
 * @mainpage 
 * This is an implementation of an imaginary time evolution algorithm for the optimization of PEPS.
 * @author Brecht Verstichel
 * @date 25-03-2014
 */
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

using namespace btas;

int main(int argc,char *argv[]){

   cout.precision(15);

   int L = atoi(argv[1]);//dimension of the lattice: LxL
   int d = atoi(argv[2]);//physical dimension
   int D = atoi(argv[3]);//virtual dimension
   int D_aux = atoi(argv[4]);//virtual dimension
   int J2 = atoi(argv[5]);

   double tau = 0.01;

   //initialize some statics dimensions
   global::init(D,D_aux,d,L,L,J2,tau);

   PEPS<double> peps(D);
   peps.normalize();

   for(int i = 0;i < 1000;++i){

      propagate::step(peps,1);

      peps.rescale_tensors(1.0);
      peps.normalize();

      global::env.calc('A',peps); 
      cout << i << "\t" << peps.energy() << endl;

   }

   return 0;

}
