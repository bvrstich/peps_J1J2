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

   DArray<6> tmp6(2,3,4,5,6,7);
   tmp6.generate(global::rgen<double>);

   Perm<6> perm(tmp6.shape(),shape(1,4,3,2,0,5));

   perm.permute(tmp6);

   cout << perm.gperm_tensor() << endl;

   DArray<6> tmp6bis;

   Permute(tmp6,shape(1,4,3,2,0,5),tmp6bis);
   ofstream out("old.out");
   out.precision(15);
   out << tmp6bis << endl;
  
   return 0;

}
