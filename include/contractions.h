#ifndef CONTRACTIONS_H
#define CONTRACTIONS_H

#include <iostream>
#include <iomanip>

#include <btas/common/blas_cxx_interface.h>
#include <btas/common/TVector.h>
#include <btas/DENSE/TArray.h>

using namespace btas;

//some much repeated contractions are put in separate functions here
namespace contractions {

   void init_ro(char option,const PEPS<double> &,vector< DArray<5> > &R);

   void init_ro(int row,const PEPS<double> &,vector< DArray<6> > &RO);

   double rescale_norm(int row,PEPS<double> &,vector< DArray<6> > &RO);

   void update_L(char option,int col,const PEPS<double> &,DArray<5> &L);

   void update_L(int row,int col,const PEPS<double> &,DArray<6> &LO);

}

#endif
