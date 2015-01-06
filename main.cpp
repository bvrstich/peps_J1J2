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
#include <chrono>
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

   bool update = true;
   double tau = 0.01;
   int n_steps = 10;

   //initialize some statics dimensions
   global::init(D,D_aux,d,L,L,J2,tau);

   DArray<9> tmp(4,16,4,4,4,4,4,4,4);
   Perm<9> perm(tmp.shape(),shape(0,1,2,3,5,4,6,7,8));

   DArray<9> tmpbis;

   double tijd = 0.0;

   for(int i = 0;i < 10000;++i){

      tmp.generate(global::rgen<double>);

      auto start = std::chrono::high_resolution_clock::now();
      //Permute(tmp,shape(0,1,2,3,4,5,6,8,7),tmpbis);
      perm.permute(tmp);
      auto end = std::chrono::high_resolution_clock::now();

      tijd += std::chrono::duration_cast<std::chrono::duration<double,std::ratio<1>>>(end-start).count();

   }

   cout << tijd << endl;

   return 0;

}
