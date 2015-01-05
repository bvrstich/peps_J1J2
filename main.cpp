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

   bool update = true;
   double tau = 0.01;
   int n_steps = 10;

   //initialize some statics dimensions
   global::init(D,D_aux,d,L,L,J2,tau);

   DArray<2> tmp2(5,5);
   tmp2.generate(global::rgen<double>);

   Perm<2> perm(tmp2.shape(),shape(1,0));

   perm.permute(tmp2);

   cout << perm.gperm_tensor() << endl;

   DArray<2> tmp2bis;

   Permute(tmp2,shape(1,0),tmp2bis);
   cout << tmp2bis << endl;
  
   return 0;

}
