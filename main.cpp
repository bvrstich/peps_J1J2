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

   bool update = true;
   double tau = 0.01;
   int n_steps = 10;

   int mult = D_aux / (D*D);

   //initialize some statics dimensions
   global::init(D,D_aux,d,L,L,tau);

   PEPS<double> peps(D);
   peps.normalize();
   //peps.initialize_jastrow(0.74);

   for(int i= 0;i < 5000;++i){

      propagate::step(update,peps,10);

      if(i % 10 == 0){

         peps.normalize();

         global::env.calc('A',peps);
         cout << i << "\t" << peps.energy()/(L*L) << endl;

      }

   }

   tau /= 10.0;
   global::stau(tau);

   for(int i= 5000;i < 15000;++i){

      propagate::step(update,peps,10);

      if(i % 10 == 0){

         peps.normalize();

         global::env.calc('A',peps);
         cout << i << "\t" << peps.energy()/(L*L) << endl;

      }

   }

   tau /= 10.0;
   global::stau(tau);

   for(int i= 15000;i < 30000;++i){

      propagate::step(update,peps,10);

      if(i % 10 == 0){

         peps.normalize();

         global::env.calc('A',peps);
         cout << i << "\t" << peps.energy()/(L*L) << endl;

      }

   }

   cout << peps << endl;

   return 0;

}
