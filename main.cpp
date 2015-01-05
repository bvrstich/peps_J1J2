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

   DArray<9> tmp9(2,3,4,5,6,7,8,9,10);
   tmp9.generate(global::rgen<double>);

   Perm<9> perm(tmp9.shape(),shape(1,2,7,4,5,6,3,8,0));

   perm.permute(tmp9);

   cout << perm.gperm_tensor() << endl;

   DArray<9> tmp9bis;

   Permute(tmp9,shape(1,2,7,4,5,6,3,8,0),tmp9bis);
   ofstream out("old.out");
   out.precision(15);
   out << tmp9bis << endl;
    
   return 0;

}
