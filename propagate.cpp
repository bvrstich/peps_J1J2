#include <iostream>
#include <iomanip>
#include <fstream>
#include <complex>

using std::cout;
using std::endl;
using std::ostream;
using std::complex;
using std::ofstream;

#include "include.h"

using namespace btas;
using namespace global;

namespace propagate {

   /**
    * propagate the peps one imaginary time step
    * @param full if true do full update, if false do simple update...
    * @param peps the PEPS to be propagated
    * @param n_sweeps the number of sweeps performed for the solution of the linear problem
    */
   void step(bool full,PEPS<double> &peps,int n_sweeps){

      if(ham.gis_local())
         prop_local(peps);

      enum {i,j,k,l,m,n,o};

      // ########################################################### //
      // ########################################################### //
      // ##                                                       ## //
      // ## First propagate applying the gates from bottom to top ## //
      // ##                                                       ## //
      // ########################################################### //
      // ########################################################### //

      // --------------------------------------//
      // --- !!! (1) the bottom row (1) !!! ---// 
      // --------------------------------------//

      //containers for the renormalized operators
      vector< DArray<3> > R(Lx - 1);

      DArray<3> L;

      DArray<4> QL;
      DArray<3> a_L;

      DArray<4> QR;
      DArray<3> a_R;

      DArray<4> N_eff;

      //construct the full top environment:
      env.calc('T',peps);

      //and the bottom row environment
      env.gb(0).fill('b',peps);

      //initialize the right operators for the bottom row
      contractions::init_ro(false,'b',peps,R);

      //middle sites of the bottom row:
      for(int col = 0;col < Lx - 1;++col){

         //first construct the reduced tensors of the first pair to propagate
         construct_reduced_tensor('H','L',peps(0,col),QL,a_L);
         construct_reduced_tensor('H','R',peps(0,col+1),QR,a_R);

         //calculate the effective environment N_eff
         calc_N_eff('b',col,L,QL,R[col + 1],QR,N_eff);

         //make environment close to unitary before the update
         canonicalize(full,N_eff,a_L,QL,a_R,QR);

         //now do the update! Apply the gates!
         update(full,N_eff,a_L,a_R,n_sweeps);

         //and expand back to the full tensors
         Contract(1.0,QL,shape(i,j,k,o),a_L,shape(o,m,n),0.0,peps(0,col),shape(i,j,m,k,n));
         Contract(1.0,a_R,shape(i,j,k),QR,shape(k,o,m,n),0.0,peps(0,col+1),shape(i,o,j,m,n));

         contractions::update_L('b',col,peps,L);

      }

      //update the bottom row for the new peps
      env.gb(0).fill('b',peps);

      // ---------------------------------------------------//
      // --- !!! (2) the middle rows (1 -> Ly-2) (2) !!! ---// 
      // ---------------------------------------------------//

      //renormalized operators for the middle sites
      vector< DArray<4> > RO(Lx - 1);
      DArray<4> LO;

      for(int row = 1;row < Ly-1;++row){

         //first create right renormalized operator
         contractions::init_ro(false,'H',row,peps,RO);

         for(int col = 0;col < Lx - 1;++col){

            //first construct the reduced tensors of the first pair to propagate
            construct_reduced_tensor('H','L',peps(row,col),QL,a_L);
            construct_reduced_tensor('H','R',peps(row,col+1),QR,a_R);

            //calculate the effective environment N_eff
            calc_N_eff('H',row,col,LO,QL,RO[col + 1],QR,N_eff);

            //make environment close to unitary before the update
            canonicalize(full,N_eff,a_L,QL,a_R,QR);

            //now do the update! Apply the gates!
            update(full,N_eff,a_L,a_R,n_sweeps);

            //and expand back to the full tensors
            Contract(1.0,QL,shape(i,j,k,o),a_L,shape(o,m,n),0.0,peps(row,col),shape(i,j,m,k,n));
            Contract(1.0,a_R,shape(i,j,k),QR,shape(k,o,m,n),0.0,peps(row,col+1),shape(i,o,j,m,n));

            //first construct a double layer object for the newly updated bottom 
            contractions::update_L('H',row,col,peps,LO);

         }

         //finally update the 'bottom' environment for the row
         env.add_layer('b',row,peps);

      }

      // ------------------------------------------//
      // --- !!! (3) the top row (Ly-1) (3) !!! ---// 
      // ------------------------------------------//

      //make the right operators
      contractions::init_ro(false,'t',peps,R);

      for(int col = 0;col < Lx - 1;++col){

         //first construct the reduced tensors of the first pair to propagate
         construct_reduced_tensor('H','L',peps(Ly-1,col),QL,a_L);
         construct_reduced_tensor('H','R',peps(Ly-1,col + 1),QR,a_R);

         //calculate the effective environment N_eff
         calc_N_eff('t',col,L,QL,R[col + 1],QR,N_eff);

         //make environment close to unitary before the update
         canonicalize(full,N_eff,a_L,QL,a_R,QR);

         //now do the update! Apply the gates!
         update(full,N_eff,a_L,a_R,n_sweeps);

         //and expand back to the full tensors
         Contract(1.0,QL,shape(i,j,k,o),a_L,shape(o,m,n),0.0,peps(Ly-1,col),shape(i,j,m,k,n));
         Contract(1.0,a_R,shape(i,j,k),QR,shape(k,o,m,n),0.0,peps(Ly-1,col+1),shape(i,o,j,m,n));

         contractions::update_L('t',col,peps,L);

      }

      //get the norm matrix
      contractions::update_L('t',Lx-1,peps,L);

      //scale the peps
      peps.scal(1.0/sqrt(L(0,0,0)));

      // ########################################################## //
      // ########################################################## //
      // ##                                                      ## //
      // ## Then propagate applying the gates from right to left ## //
      // ##                                                      ## //
      // ########################################################## //
      // ########################################################## //

      // ----------------------------------------//
      // --- !!! (1) the right column (1) !!! ---// 
      // ----------------------------------------//

      //construct the full left environment:
      env.calc('L',peps);

      //and the rightmost column environment
      env.gr(Ly-2).fill('r',peps);

      //initialize the right operators for the right column
      contractions::init_ro(false,'r',peps,R);

      //middle sites of the bottom row:
      for(int row = 0;row < Ly - 1;++row){

         //first construct the reduced tensors of the first pair to propagate
         construct_reduced_tensor('V','L',peps(row,Lx-1),QL,a_L);
         construct_reduced_tensor('V','R',peps(row+1,Lx-1),QR,a_R);

         //calculate the effective environment N_eff
         calc_N_eff('r',row,L,QL,R[row + 1],QR,N_eff);

         //make environment close to unitary before the update
         canonicalize(full,N_eff,a_L,QL,a_R,QR);

         //now do the update! Apply the gates!
         update(full,N_eff,a_L,a_R,n_sweeps);

         //and expand back to the full tensors
         Contract(1.0,QL,shape(i,j,k,o),a_L,shape(o,m,n),0.0,peps(row,Lx-1),shape(j,n,m,i,k));
         Contract(1.0,a_R,shape(i,j,k),QR,shape(k,o,m,n),0.0,peps(row+1,Lx-1),shape(o,n,j,i,m));

         contractions::update_L('r',row,peps,L);

      }

      //update the bottom row for the new peps
      env.gr(Lx-2).fill('r',peps);

      // ---------------------------------------------------//
      // --- !!! (2) the middle cols (Lx-2 -> 1) (2) !!! ---// 
      // ---------------------------------------------------//

      for(int col = Lx - 2;col > 0;--col){

         //first create right renormalized operator
         contractions::init_ro(false,'V',col,peps,RO);

         for(int row = 0;row < Ly - 1;++row){

            //first construct the reduced tensors of the first pair to propagate
            construct_reduced_tensor('V','L',peps(row,col),QL,a_L);
            construct_reduced_tensor('V','R',peps(row+1,col),QR,a_R);

            //calculate the effective environment N_eff
            calc_N_eff('V',row,col,LO,QL,RO[row + 1],QR,N_eff);

            //make environment close to unitary before the update
            canonicalize(full,N_eff,a_L,QL,a_R,QR);

            //now do the update! Apply the gates!
            update(full,N_eff,a_L,a_R,n_sweeps);

            //and expand back to the full tensors
            Contract(1.0,QL,shape(i,j,k,o),a_L,shape(o,m,n),0.0,peps(row,col),shape(j,n,m,i,k));
            Contract(1.0,a_R,shape(i,j,k),QR,shape(k,o,m,n),0.0,peps(row+1,col),shape(o,n,j,i,m));

            //first construct a double layer object for the newly updated bottom 
            contractions::update_L('V',row,col,peps,LO);

         }

         //finally update the 'bottom' environment for the row
         env.add_layer('r',col-1,peps);

      }

      // -----------------------------------------------//
      // --- !!! (3) the leftmost column (0) (3) !!! ---// 
      // -----------------------------------------------//

      //make the right operators
      contractions::init_ro(false,'l',peps,R);

      for(int row = 0;row < Ly - 1;++row){

         //first construct the reduced tensors of the first pair to propagate
         construct_reduced_tensor('V','L',peps(row,0),QL,a_L);
         construct_reduced_tensor('V','R',peps(row + 1,0),QR,a_R);

         //calculate the effective environment N_eff
         calc_N_eff('l',row,L,QL,R[row + 1],QR,N_eff);

         //make environment close to unitary before the update
         canonicalize(full,N_eff,a_L,QL,a_R,QR);

         //now do the update! Apply the gates!
         update(full,N_eff,a_L,a_R,n_sweeps);

         //and expand back to the full tensors
         Contract(1.0,QL,shape(i,j,k,o),a_L,shape(o,m,n),0.0,peps(row,0),shape(j,n,m,i,k));
         Contract(1.0,a_R,shape(i,j,k),QR,shape(k,o,m,n),0.0,peps(row+1,0),shape(o,n,j,i,m));

         contractions::update_L('l',row,peps,L);

      }

      //get the norm matrix
      contractions::update_L('l',Ly-1,peps,L);

      //scale the peps
      peps.scal(1.0/sqrt(L(0,0,0)));

      if(ham.gis_local())
         prop_local(peps);

   }

   /**
    * construct the left or right reduced tensor form of a peps element by performing QR or LQ decomposition
    * @param 'H'orizontal or 'V'ertical
    * @param L == left, R == right
    */
   void construct_reduced_tensor(char hv,char option,const DArray<5> &peps,DArray<4> &Q,DArray<3> &red){

      if(hv == 'H'){

         if(option == 'L'){

            DArray<5> tmp;
            Permute(peps,shape(0,1,3,2,4),tmp);

            int nrows = tmp.shape(0) * tmp.shape(1) * tmp.shape(2);
            int ncols = tmp.shape(3) * tmp.shape(4);

            int min = std::min(nrows,ncols);

            double* tau = new double [min];

            lapack::geqrf(CblasRowMajor,nrows,ncols, tmp.data(), ncols, tau);

            red.resize(shape(min,tmp.shape(3),tmp.shape(4)));
            Q.resize(tmp.shape(0),tmp.shape(1),tmp.shape(2),min);

            red = 0.0;

            //r is in the upper diagonal part of tmp on exit of geqrf:
            for(int i = 0;i < min;++i)
               for(int j = i;j < ncols;++j)
                  red.data()[i*ncols + j] = tmp.data()[i*ncols + j];

            //get the input for the Q construction function: lower diagonal part of the matrix
            for(int i = 0;i < nrows;++i)
               for(int j = 0;j < min;++j)
                  Q.data()[i*min + j] = tmp.data()[i*ncols + j];

            //now get the Q matrix out
            if(nrows < ncols)
               lapack::orgqr(CblasRowMajor, nrows, nrows, min,Q.data(), nrows, tau);
            else
               lapack::orgqr(CblasRowMajor, nrows, ncols, min,Q.data(), ncols, tau);

            delete [] tau;

         }
         else{//LQ

            DArray<5> tmp;
            Permute(peps,shape(0,2,1,3,4),tmp);

            int nrows = tmp.shape(0) * tmp.shape(1);
            int ncols = tmp.shape(2) * tmp.shape(3) * tmp.shape(4);

            int min = std::min(nrows,ncols);

            double* tau = new double [min];

            lapack::gelqf(CblasRowMajor,nrows,ncols,tmp.data(),ncols,tau);

            red.resize(shape(tmp.shape(0),tmp.shape(1),min));
            Q.resize(min,tmp.shape(2),tmp.shape(3),tmp.shape(4));

            red = 0.0;

            //l is in the lower diagonal part of tmp on exit of gelqf:
            for(int j = 0;j < ncols;++j)
               for(int i = j;i < nrows;++i)
                  red.data()[i*min + j] = tmp.data()[i*ncols + j];

            //get the input for the Q construction function: upper diagonal part of the matrix
            for(int i = 0;i < min;++i)
               for(int j = 0;j < ncols;++j)
                  Q.data()[i*ncols + j] = tmp.data()[i*ncols + j];

            //now get the Q matrix out
            if(nrows > ncols)
               lapack::orglq(CblasRowMajor,ncols,ncols,min,Q.data(),ncols,tau);
            else
               lapack::orglq(CblasRowMajor,nrows,ncols,min,Q.data(),ncols,tau);

            delete [] tau;

         }

      }
      else{//Vertical gates

         if(option == 'L'){

            DArray<5> tmp;
            Permute(peps,shape(3,0,4,2,1),tmp);

            int nrows = tmp.shape(0) * tmp.shape(1) * tmp.shape(2);
            int ncols = tmp.shape(3) * tmp.shape(4);

            int min = std::min(nrows,ncols);

            double* tau = new double [min];

            lapack::geqrf(CblasRowMajor,nrows,ncols, tmp.data(), ncols, tau);

            red.resize(shape(min,tmp.shape(3),tmp.shape(4)));
            Q.resize(tmp.shape(0),tmp.shape(1),tmp.shape(2),min);

            red = 0.0;

            //r is in the upper diagonal part of tmp on exit of geqrf:
            for(int i = 0;i < min;++i)
               for(int j = i;j < ncols;++j)
                  red.data()[i*ncols + j] = tmp.data()[i*ncols + j];

            //get the input for the Q construction function: lower diagonal part of the matrix
            for(int i = 0;i < nrows;++i)
               for(int j = 0;j < min;++j)
                  Q.data()[i*min + j] = tmp.data()[i*ncols + j];

            //now get the Q matrix out
            if(nrows < ncols)
               lapack::orgqr(CblasRowMajor, nrows, nrows, min,Q.data(), nrows, tau);
            else
               lapack::orgqr(CblasRowMajor, nrows, ncols, min,Q.data(), ncols, tau);

            delete [] tau;

         }
         else{//LQ

            DArray<5> tmp;
            Permute(peps,shape(3,2,0,4,1),tmp);

            int nrows = tmp.shape(0) * tmp.shape(1);
            int ncols = tmp.shape(2) * tmp.shape(3) * tmp.shape(4);

            int min = std::min(nrows,ncols);

            double* tau = new double [min];

            lapack::gelqf(CblasRowMajor,nrows,ncols,tmp.data(),ncols,tau);

            red.resize(shape(tmp.shape(0),tmp.shape(1),min));
            Q.resize(min,tmp.shape(2),tmp.shape(3),tmp.shape(4));

            red = 0.0;

            //l is in the lower diagonal part of tmp on exit of gelqf:
            for(int j = 0;j < ncols;++j)
               for(int i = j;i < nrows;++i)
                  red.data()[i*min + j] = tmp.data()[i*ncols + j];

            //get the input for the Q construction function: upper diagonal part of the matrix
            for(int i = 0;i < min;++i)
               for(int j = 0;j < ncols;++j)
                  Q.data()[i*ncols + j] = tmp.data()[i*ncols + j];

            //now get the Q matrix out
            if(nrows > ncols)
               lapack::orglq(CblasRowMajor,ncols,ncols,min,Q.data(),ncols,tau);
            else
               lapack::orglq(CblasRowMajor,nrows,ncols,min,Q.data(),ncols,tau);

            delete [] tau;

         }

      }

   }

   /** 
    * wrapper function invert square general matrix DArray<2>.
    * @param A both input as output matrix: on input A, on output A^{-1}
    */
   void invert(DArray<2> &A){

      int *ipiv = new int [A.shape(0)];

      lapack::getrf(CblasRowMajor,A.shape(0),A.shape(1), A.data(), A.shape(1), ipiv);

      lapack::getri(CblasRowMajor,A.shape(0), A.data(), A.shape(1), ipiv);

      delete [] ipiv;

   }

   /** 
    * wrapper function solve positive symmetric linear system: N_eff * x = b
    * @param N_eff input matrix
    * @param b right hand side input
    * @param x output
    */
   void solve(DArray<4> &N_eff,DArray<3> &b){

      int n = N_eff.shape(0) * N_eff.shape(1);

      lapack::potrf(CblasRowMajor,'U',n, N_eff.data(), n);

      lapack::potrs(CblasRowMajor,'U',n,d, N_eff.data(), n,b.data(),d);

   }

   /**
    * calculate the effective environment of a pair with the left tensor on site (row,col)
    * @param option 't'op ,'b'ottom row or 'l'eft, 'r'ight column
    * @param rc column or row index of the left tensor
    * @param L left environment matrix
    * @param QL left unitary matrix coming out of the reduced tensor construction
    * @param R left environment matrix
    * @param QR right unitary matrix coming out of the reduced tensor construction
    * @param N_eff output DArray<4> object containing the effective norm environment on exit
    */
   void calc_N_eff(char option,int rc,const DArray<3> &L,const DArray<4> &QL,const DArray<3> &R, const DArray<4> &QR,DArray<4> &N_eff){

      if(option == 'b'){

         if(rc == 0){//left edge

            //for this one only top contraction is needed:
            DArray<6> tmp6;
            Contract(1.0,env.gt(0)[0],shape(1),QL,shape(1),0.0,tmp6);

            DArray<8> tmp8;
            Contract(1.0,tmp6,shape(1),QL,shape(1),0.0,tmp8);

            //construct the 'Left' eff environment
            DArray<3> L_env = tmp8.reshape_clear(shape(tmp8.shape(1),tmp8.shape(4),tmp8.shape(7)));

            //contract with right renormalized operator:
            DArray<5> tmp5;
            Contract(1.0,env.gt(0)[1],shape(3),R,shape(0),0.0,tmp5);

            //to construct the R_environment
            DArray<5> tmp5bis;
            Contract(1.0,tmp5,shape(1,3),QR,shape(1,3),0.0,tmp5bis);

            tmp5.clear();
            Contract(1.0,tmp5bis,shape(1,2),QR,shape(1,3),0.0,tmp5);

            //construct the 'Right' eff environment
            DArray<3> R_env = tmp5.reshape_clear(shape(tmp5.shape(0),tmp5.shape(1),tmp5.shape(3)));

            //now contract left and right environment to form N_eff
            enum {i,j,k,m,n};

            N_eff.clear();
            Contract(1.0,L_env,shape(i,j,k),R_env,shape(i,n,m),0.0,N_eff,shape(j,n,k,m));

         }
         else if(rc == Lx - 2){//right edge

            //contraction with left renormalized operator
            DArray<5> tmp5;
            Contract(1.0,L,shape(0),env.gt(0)[rc],shape(0),0.0,tmp5);

            DArray<5> tmp5bis;
            Contract(1.0,tmp5,shape(0,2),QL,shape(0,1),0.0,tmp5bis);

            tmp5.clear();
            Contract(1.0,tmp5bis,shape(0,1),QL,shape(0,1),0.0,tmp5);

            //construct the 'Left' eff environment
            DArray<3> L_env = tmp5.reshape_clear( shape(tmp5.shape(0),tmp5.shape(2),tmp5.shape(4)) );

            //Right

            //only attach to top to construct R_env
            DArray<6> tmp6;
            Contract(1.0,env.gt(0)[Lx-1],shape(1),QR,shape(1),0.0,tmp6);

            DArray<8> tmp8;
            Contract(1.0,tmp6,shape(1),QR,shape(1),0.0,tmp8);

            DArray<3> R_env = tmp8.reshape_clear( shape( tmp8.shape(0),tmp8.shape(2),tmp8.shape(5) ));

            //construct effective environment
            enum {i,j,k,m,n};

            N_eff.clear();
            Contract(1.0,L_env,shape(i,j,k),R_env,shape(i,m,n),0.0,N_eff,shape(j,m,k,n));

         }
         else{//middle

            //contraction with left renormalized operator
            DArray<5> tmp5;
            Contract(1.0,L,shape(0),env.gt(0)[rc],shape(0),0.0,tmp5);

            DArray<5> tmp5bis;
            Contract(1.0,tmp5,shape(0,2),QL,shape(0,1),0.0,tmp5bis);

            tmp5.clear();
            Contract(1.0,tmp5bis,shape(0,1),QL,shape(0,1),0.0,tmp5);

            //construct the 'Left' eff environment
            DArray<3> L_env = tmp5.reshape_clear( shape(tmp5.shape(0),tmp5.shape(2),tmp5.shape(4)) );

            //contract with right renormalized operator:
            tmp5.clear();
            Contract(1.0,env.gt(0)[rc+1],shape(3),R,shape(0),0.0,tmp5);

            //to construct the R_environment
            tmp5bis.clear();
            Contract(1.0,tmp5,shape(1,3),QR,shape(1,3),0.0,tmp5bis);

            tmp5.clear();
            Contract(1.0,tmp5bis,shape(1,2),QR,shape(1,3),0.0,tmp5);

            //construct the 'Right' eff environment
            DArray<3> R_env = tmp5.reshape_clear(shape(tmp5.shape(0),tmp5.shape(1),tmp5.shape(3)));

            //now contract left and right environment to form N_eff
            enum {i,j,k,l,m};

            N_eff.clear();
            Contract(1.0,L_env,shape(i,j,k),R_env,shape(i,l,m),0.0,N_eff,shape(j,l,k,m));

         }

      }
      else if(option == 't'){//top row!

         if(rc == 0){

            //Left

            //for this one only bottom contraction is needed:
            DArray<6> tmp6;
            Contract(1.0,QL,shape(2),env.gb(Ly-2)[0],shape(2),0.0,tmp6);

            //construct the 'Left' eff environment
            DArray<8> tmp8;
            Contract(1.0,QL,shape(2),tmp6,shape(4),0.0,tmp8);

            DArray<3> L_env = tmp8.reshape_clear(shape( tmp8.shape(2),tmp8.shape(5),tmp8.shape(7) ) );

            //Right

            //contract with right renormalized operator:
            DArray<5> tmp5;
            Contract(1.0,env.gb(Ly-2)[1],shape(3),R,shape(2),0.0,tmp5);

            DArray<5> tmp5bis;
            Contract(1.0,QR,shape(2,3),tmp5,shape(2,4),0.0,tmp5bis);

            tmp5.clear();
            Contract(1.0,QR,shape(2,3),tmp5bis,shape(3,4),0.0,tmp5);

            DArray<3> R_env = tmp5.reshape_clear( shape(tmp5.shape(0),tmp5.shape(2),tmp5.shape(4)) );

            //now contract left and right environment to form N_eff
            enum {i,j,k,l,m};

            N_eff.clear();
            Contract(1.0,L_env,shape(i,j,k),R_env,shape(l,m,k),0.0,N_eff,shape(i,l,j,m));

         }
         else if(rc == Lx - 2){//right edge of top row

            //Left
            DArray<5> tmp5;
            Contract(1.0,L,shape(2),env.gb(Ly-2)[Lx - 2],shape(0),0.0,tmp5);

            DArray<5> tmp5bis;
            Contract(1.0,QL,shape(0,2),tmp5,shape(1,3),0.0,tmp5bis);

            tmp5.clear();
            Contract(1.0,QL,shape(0,2),tmp5bis,shape(2,3),0.0,tmp5);

            DArray<3> L_env = tmp5.reshape_clear( shape(tmp5.shape(1),tmp5.shape(3),tmp5.shape(4)) );

            //Right

            //for this one only bottom contraction is needed:
            DArray<6> tmp6;
            Contract(1.0,QR,shape(2),env.gb(Ly-2)[Lx - 1],shape(2),0.0,tmp6);

            //construct the 'Left' eff environment
            DArray<8> tmp8;
            Contract(1.0,QR,shape(2),tmp6,shape(4),0.0,tmp8);

            DArray<3> R_env = tmp8.reshape_clear(shape( tmp8.shape(0),tmp8.shape(3),tmp8.shape(6) ) );

            //construct effective environment
            enum {i,j,k,l,m};

            N_eff.clear();
            Contract(1.0,L_env,shape(i,j,k),R_env,shape(l,m,k),0.0,N_eff,shape(i,l,j,m));

         }
         else{//middle columns

            //Left
            DArray<5> tmp5;
            Contract(1.0,L,shape(2),env.gb(Ly-2)[rc],shape(0),0.0,tmp5);

            DArray<5> tmp5bis;
            Contract(1.0,QL,shape(0,2),tmp5,shape(1,3),0.0,tmp5bis);

            tmp5.clear();
            Contract(1.0,QL,shape(0,2),tmp5bis,shape(2,3),0.0,tmp5);

            DArray<3> L_env = tmp5.reshape_clear( shape(tmp5.shape(1),tmp5.shape(3),tmp5.shape(4)) );

            //Right

            //contract with right renormalized operator:
            tmp5.clear();
            Contract(1.0,env.gb(Ly-2)[rc + 1],shape(3),R,shape(2),0.0,tmp5);

            tmp5bis.clear();
            Contract(1.0,QR,shape(2,3),tmp5,shape(2,4),0.0,tmp5bis);

            tmp5.clear();
            Contract(1.0,QR,shape(2,3),tmp5bis,shape(3,4),0.0,tmp5);

            DArray<3> R_env = tmp5.reshape_clear( shape(tmp5.shape(0),tmp5.shape(2),tmp5.shape(4)) );

            //now contract left and right environment to form N_eff
            enum {i,j,k,l,m};

            N_eff.clear();
            Contract(1.0,L_env,shape(i,j,k),R_env,shape(l,m,k),0.0,N_eff,shape(i,l,j,m));

         }

      }
      else if(option == 'l'){//left

         if(rc == 0){//left edge

            //Left

            //for this one only right contraction is needed:
            DArray<6> tmp6;
            Contract(1.0,QL,shape(2),env.gr(0)[0],shape(2),0.0,tmp6);

            //construct the 'Left' eff environment
            DArray<8> tmp8;
            Contract(1.0,QL,shape(2),tmp6,shape(4),0.0,tmp8);

            DArray<3> L_env = tmp8.reshape_clear(shape( tmp8.shape(2),tmp8.shape(5),tmp8.shape(7) ) );

            //Right

            //contract with right renormalized operator:
            DArray<5> tmp5;
            Contract(1.0,env.gr(0)[1],shape(3),R,shape(2),0.0,tmp5);

            DArray<5> tmp5bis;
            Contract(1.0,QR,shape(2,3),tmp5,shape(2,4),0.0,tmp5bis);

            tmp5.clear();
            Contract(1.0,QR,shape(2,3),tmp5bis,shape(3,4),0.0,tmp5);

            DArray<3> R_env = tmp5.reshape_clear( shape(tmp5.shape(0),tmp5.shape(2),tmp5.shape(4)) );

            //now contract left and right environment to form N_eff
            enum {i,j,k,l,m};

            N_eff.clear();
            Contract(1.0,L_env,shape(i,j,k),R_env,shape(l,m,k),0.0,N_eff,shape(i,l,j,m));

         }
         else if(rc == Ly - 2){//'right' edge of left column

            //Left
            DArray<5> tmp5;
            Contract(1.0,L,shape(2),env.gr(0)[Ly - 2],shape(0),0.0,tmp5);

            DArray<5> tmp5bis;
            Contract(1.0,QL,shape(0,2),tmp5,shape(1,3),0.0,tmp5bis);

            tmp5.clear();
            Contract(1.0,QL,shape(0,2),tmp5bis,shape(2,3),0.0,tmp5);

            DArray<3> L_env = tmp5.reshape_clear( shape(tmp5.shape(1),tmp5.shape(3),tmp5.shape(4)) );

            //Right

            //for this one only right contraction is needed:
            DArray<6> tmp6;
            Contract(1.0,QR,shape(2),env.gr(0)[Ly - 1],shape(2),0.0,tmp6);

            //construct the 'Left' eff environment
            DArray<8> tmp8;
            Contract(1.0,QR,shape(2),tmp6,shape(4),0.0,tmp8);

            DArray<3> R_env = tmp8.reshape_clear(shape( tmp8.shape(0),tmp8.shape(3),tmp8.shape(6) ) );

            //construct effective environment
            enum {i,j,k,l,m};

            N_eff.clear();
            Contract(1.0,L_env,shape(i,j,k),R_env,shape(l,m,k),0.0,N_eff,shape(i,l,j,m));

         }
         else{//middle

            //Left
            DArray<5> tmp5;
            Contract(1.0,L,shape(2),env.gr(0)[rc],shape(0),0.0,tmp5);

            DArray<5> tmp5bis;
            Contract(1.0,QL,shape(0,2),tmp5,shape(1,3),0.0,tmp5bis);

            tmp5.clear();
            Contract(1.0,QL,shape(0,2),tmp5bis,shape(2,3),0.0,tmp5);

            DArray<3> L_env = tmp5.reshape_clear( shape(tmp5.shape(1),tmp5.shape(3),tmp5.shape(4)) );

            //Right

            //contract with right renormalized operator:
            tmp5.clear();
            Contract(1.0,env.gr(0)[rc + 1],shape(3),R,shape(2),0.0,tmp5);

            tmp5bis.clear();
            Contract(1.0,QR,shape(2,3),tmp5,shape(2,4),0.0,tmp5bis);

            tmp5.clear();
            Contract(1.0,QR,shape(2,3),tmp5bis,shape(3,4),0.0,tmp5);

            DArray<3> R_env = tmp5.reshape_clear( shape(tmp5.shape(0),tmp5.shape(2),tmp5.shape(4)) );

            //now contract left and right environment to form N_eff
            enum {i,j,k,l,m};

            N_eff.clear();
            Contract(1.0,L_env,shape(i,j,k),R_env,shape(l,m,k),0.0,N_eff,shape(i,l,j,m));

         }

      }
      else{//rightmost column

         if(rc == 0){

            //for this one only left contraction is needed:
            DArray<6> tmp6;
            Contract(1.0,env.gl(Lx-2)[0],shape(1),QL,shape(1),0.0,tmp6);

            DArray<8> tmp8;
            Contract(1.0,tmp6,shape(1),QL,shape(1),0.0,tmp8);

            //construct the 'Left' eff environment
            DArray<3> L_env = tmp8.reshape_clear(shape(tmp8.shape(1),tmp8.shape(4),tmp8.shape(7)));

            //contract with right renormalized operator:
            DArray<5> tmp5;
            Contract(1.0,env.gl(Lx-2)[1],shape(3),R,shape(0),0.0,tmp5);

            //to construct the R_environment
            DArray<5> tmp5bis;
            Contract(1.0,tmp5,shape(1,3),QR,shape(1,3),0.0,tmp5bis);

            tmp5.clear();
            Contract(1.0,tmp5bis,shape(1,2),QR,shape(1,3),0.0,tmp5);

            //construct the 'Right' eff environment
            DArray<3> R_env = tmp5.reshape_clear(shape(tmp5.shape(0),tmp5.shape(1),tmp5.shape(3)));

            //now contract left and right environment to form N_eff
            enum {i,j,k,m,n};

            N_eff.clear();
            Contract(1.0,L_env,shape(i,j,k),R_env,shape(i,n,m),0.0,N_eff,shape(j,n,k,m));

         }
         else if(rc == Lx-2){//right edge of top row

            //contraction with left renormalized operator
            DArray<5> tmp5;
            Contract(1.0,L,shape(0),env.gl(Lx-2)[rc],shape(0),0.0,tmp5);

            DArray<5> tmp5bis;
            Contract(1.0,tmp5,shape(0,2),QL,shape(0,1),0.0,tmp5bis);

            tmp5.clear();
            Contract(1.0,tmp5bis,shape(0,1),QL,shape(0,1),0.0,tmp5);

            //construct the 'Left' eff environment
            DArray<3> L_env = tmp5.reshape_clear( shape(tmp5.shape(0),tmp5.shape(2),tmp5.shape(4)) );

            //Right

            //only attach to top to construct R_env
            DArray<6> tmp6;
            Contract(1.0,env.gl(Lx-2)[Ly-1],shape(1),QR,shape(1),0.0,tmp6);

            DArray<8> tmp8;
            Contract(1.0,tmp6,shape(1),QR,shape(1),0.0,tmp8);

            DArray<3> R_env = tmp8.reshape_clear( shape( tmp8.shape(0),tmp8.shape(2),tmp8.shape(5) ));

            //construct effective environment
            enum {i,j,k,m,n};

            N_eff.clear();
            Contract(1.0,L_env,shape(i,j,k),R_env,shape(i,m,n),0.0,N_eff,shape(j,m,k,n));

         }
         else{//middle columns

            //contraction with left renormalized operator
            DArray<5> tmp5;
            Contract(1.0,L,shape(0),env.gl(Lx-2)[rc],shape(0),0.0,tmp5);

            DArray<5> tmp5bis;
            Contract(1.0,tmp5,shape(0,2),QL,shape(0,1),0.0,tmp5bis);

            tmp5.clear();
            Contract(1.0,tmp5bis,shape(0,1),QL,shape(0,1),0.0,tmp5);

            //construct the 'Left' eff environment
            DArray<3> L_env = tmp5.reshape_clear( shape(tmp5.shape(0),tmp5.shape(2),tmp5.shape(4)) );

            //contract with right renormalized operator:
            tmp5.clear();
            Contract(1.0,env.gl(Lx-2)[rc + 1],shape(3),R,shape(0),0.0,tmp5);

            //to construct the R_environment
            tmp5bis.clear();
            Contract(1.0,tmp5,shape(1,3),QR,shape(1,3),0.0,tmp5bis);

            tmp5.clear();
            Contract(1.0,tmp5bis,shape(1,2),QR,shape(1,3),0.0,tmp5);

            //construct the 'Right' eff environment
            DArray<3> R_env = tmp5.reshape_clear(shape(tmp5.shape(0),tmp5.shape(1),tmp5.shape(3)));

            //now contract left and right environment to form N_eff
            enum {i,j,k,m,n};

            N_eff.clear();
            Contract(1.0,L_env,shape(i,j,k),R_env,shape(i,n,m),0.0,N_eff,shape(j,n,k,m));

         }

      }

   }

   /**
    * make the environment as 'canonical' as possible so that the svd for the pair update is as optimal as possible.
    * Do this by contructing the best positive symmetric approximation to the environment and taking a QR and RQ of the resultatant X^+ X = N_eff
    * @param full boolean if true, full update, if false simple
    * @param X (XX^T) is the best positive approximation to the environment of the pair
    * @param a_L left reduced tensor, will be multiplied with the L of the environment
    * @param QL unitary part of the left tensor reduction, will be multiplied with the inverse of the L of the environment
    * @param a_R right reduced tensor, will be multiplied with the R of the environment
    * @param QR unitary part of the right tensor reduction, will be multiplied with the inverse of the R of the environment
    */
   void canonicalize(bool full,DArray<4> &N_eff,DArray<3> &a_L,DArray<4> &QL,DArray<3> &a_R,DArray<4> &QR){

      int matdim = N_eff.shape(0)*N_eff.shape(1);

      //symmetrize
      for(int i = 0;i < matdim;++i)
         for(int j = i + 1;j < matdim;++j){

            N_eff.data()[i*matdim + j] = 0.5 * (N_eff.data()[i*matdim + j]  + N_eff.data()[j*matdim + i]);
            N_eff.data()[j*matdim + i] = N_eff.data()[i*matdim + j];

         }

      DArray<1> eig(matdim);

      //diagonalize
      lapack::syev(CblasRowMajor, 'V','U', matdim, N_eff.data(), matdim, eig.data());

      DArray<3> X(N_eff.shape(0),matdim,N_eff.shape(1));

      //get the square root of the positive approximant:
      for(int iL = 0;iL < N_eff.shape(0);++iL)
         for(int iR = 0;iR < N_eff.shape(1);++iR)
            for(int kL = 0;kL < N_eff.shape(0);++kL)
               for(int kR = 0;kR < N_eff.shape(1);++kR){

                  if(eig(kL*N_eff.shape(1) + kR) > 0.0)
                     X(iL,kL*N_eff.shape(1) + kR,iR) = sqrt( eig(kL*N_eff.shape(1) + kR) ) * N_eff(iL,iR,kL,kR);

               }

      //now QR and LQ the X matrix and paste it on the aR and aL
      DArray<3> X_copy(X);

      //QR
      DArray<2> tmp2;
      Geqrf(X_copy,tmp2);

      //first paste it on the reduced tensor: a_R * R
      DArray<3> tmp3;
      Contract(1.0,a_R,shape(2),tmp2,shape(1),0.0,tmp3);

      a_R = std::move(tmp3);

      //paste the inverse to the environment tensor: R^{-1} * QR
      invert(tmp2);

      DArray<4> tmp4;
      Contract(1.0,tmp2,shape(0),QR,shape(0),0.0,tmp4);

      QR = std::move(tmp4);

      //LQ
      Gelqf(tmp2,X);

      //first paste it on the reduced tensor: L * a_L
      tmp3.clear();
      Contract(1.0,tmp2,shape(0),a_L,shape(0),0.0,tmp3);

      a_L = std::move(tmp3);

      //paste the inverse to the environment tensor: QL * L^{-1}
      invert(tmp2);

      tmp4.clear();
      Contract(1.0,QL,shape(3),tmp2,shape(1),0.0,tmp4);

      QL = std::move(tmp4);

      if(full){//update N_eff

         X.clear();
         Contract(1.0,tmp2,shape(1),X_copy,shape(0),0.0,X);

         //finally get the best positive approximation to N_eff
         Contract(1.0,X,shape(1),X,shape(1),0.0,N_eff);

      }

   }

   /**
    * update a tensor pair by applying a trotter gate over their mutual bond
    * after which a svd is performed over the bond and the dimensions are set back to D
    * @param full if true, full update, if false, simple update
    */
   void update(bool full,const DArray<4> &N_eff,DArray<3> &a_L,DArray<3> &a_R,int n_iter){

      enum {i,j,k,l,m,n};

      //left
      DArray<4> tmp4;
      Contract(1.0,a_L,shape(i,j,k),trot.gLO(),shape(n,m,j),0.0,tmp4,shape(i,n,m,k));

      a_L = tmp4.reshape_clear(shape(a_L.shape(0),a_L.shape(1),a_L.shape(2)*trot.gLO().shape(1)));

      //right
      tmp4.clear();
      Contract(1.0,trot.gRO(),shape(i,j,k),a_R,shape(n,k,m),0.0,tmp4,shape(j,n,i,m));

      a_R = tmp4.reshape_clear(shape(a_R.shape(0)*trot.gRO().shape(1),a_R.shape(1),a_R.shape(2)));

      //now create 'two-site' object
      tmp4.clear();
      Contract(1.0,a_L,shape(2),a_R,shape(0),0.0,tmp4);

      //create 'b' object of linear system: N_eff A = b
      DArray<4> b;
      Contract(1.0,N_eff,shape(i,j,k,l),tmp4,shape(k,m,n,l),0.0,b,shape(i,m,n,j));

      //svd the fucker
      DArray<1> S;
      Gesvd ('S','S', tmp4, S,a_L,a_R,D);

      //take the square root of the sv's
      for(int i = 0;i < S.size();++i)
         S(i) = sqrt(S(i));

      //and multiply it left and right to the tensors
      Dimm(S,a_R);
      Dimm(a_L,S);

      if(full){

         //start sweeping
         int iter = 0;

         while(iter < n_iter){

            //construct right hand side 
            DArray<3> tmp3;
            Contract(1.0,b,shape(i,j,k,l),a_R,shape(m,k,l),0.0,tmp3,shape(i,m,j));

            //construct effective environment for left site
            DArray<5> tmp5;
            Contract(1.0,N_eff,shape(3),a_R,shape(2),0.0,tmp5);

            tmp4.clear();
            Contract(1.0,tmp5,shape(i,j,k,l,m),a_R,shape(n,m,j),0.0,tmp4,shape(i,n,k,l));

            solve(tmp4,tmp3);

            Permute(tmp3,shape(0,2,1),a_L);

            //Right block

            //construct right hand side
            tmp3.clear();
            Contract(1.0,b,shape(i,j,k,l),a_L,shape(i,j,m),0.0,tmp3,shape(m,l,k));

            //construct effective environment for left site
            tmp5.clear();
            Contract(1.0,a_L,shape(0),N_eff,shape(0),0.0,tmp5);

            tmp4.clear();
            Contract(1.0,tmp5,shape(i,j,k,l,m),a_L,shape(l,i,n),0.0,tmp4,shape(j,k,n,m));

            solve(tmp4,tmp3);

            Permute(tmp3,shape(0,2,1),a_R);

            ++iter;

         }

         //When converged, put both objects on equal footing
         tmp4.clear();
         Contract(1.0,a_L,shape(2),a_R,shape(0),0.0,tmp4);

         DArray<1> S;
         Gesvd ('S','S', tmp4, S,a_L,a_R,D);

         //take the square root of the sv's
         for(int i = 0;i < S.size();++i)
            S(i) = sqrt(S(i));

         //and multiply it left and right to the tensors
         Dimm(S,a_R);
         Dimm(a_L,S);

      }

   }

   /**
    * calculate the effective environment of a pair with the left tensor on site (row,col)
    * @param option 'H'orizontal or 'V'ertical
    * @param row index 
    * @param col index
    * @param LO left environment matrix
    * @param QL left unitary matrix coming out of the reduced tensor construction
    * @param RO left environment matrix
    * @param QR right unitary matrix coming out of the reduced tensor construction
    * @param N_eff output DArray<4> object containing the effective norm environment on exit
    */
   void calc_N_eff(char option,int row,int col,const DArray<4> &LO,const DArray<4> &QL,const DArray<4> &RO, const DArray<4> &QR,DArray<4> &N_eff){

      if(option == 'H'){

         if(col == 0){//left edge

            //Left

            //construct the 'Left' eff environment
            DArray<6> tmp6;
            Contract(1.0,env.gt(row)[0],shape(1),QL,shape(1),0.0,tmp6);

            DArray<8> tmp8;
            Contract(1.0,tmp6,shape(1),QL,shape(1),0.0,tmp8);

            DArray<8> tmp8bis;
            Contract(1.0,tmp8,shape(3,6),env.gb(row - 1)[0],shape(1,2),0.0,tmp8bis);

            DArray<4> LO_env = tmp8bis.reshape_clear( shape(tmp8bis.shape(1),tmp8bis.shape(3),tmp8bis.shape(5),tmp8bis.shape(7)) );
            //Right

            //contract with right renormalized operator:
            tmp6.clear();
            Contract(1.0,env.gt(row)[1],shape(3),RO,shape(0),0.0,tmp6);

            //to construct the R_environment
            DArray<6> tmp6bis;
            Contract(1.0,tmp6,shape(1,3),QR,shape(1,3),0.0,tmp6bis);

            tmp6.clear();
            Contract(1.0,tmp6bis,shape(1,2),QR,shape(1,3),0.0,tmp6);

            tmp6bis.clear();
            Permute(tmp6,shape(0,2,4,3,5,1),tmp6bis);

            //construct the 'Right' eff environment
            DArray<4> RO_env;
            Gemm(CblasNoTrans,CblasTrans,1.0,tmp6bis,env.gb(row - 1)[1],0.0,RO_env);

            //construct effective environment
            enum {i,j,k,l,m,n,o};

            N_eff.clear();
            Contract(1.0,LO_env,shape(i,j,k,l),RO_env,shape(i,m,n,l),0.0,N_eff,shape(j,m,k,n));

         }
         else if(col == Lx - 2){//right edge

            //Left
            //first attach top to left unity
            DArray<6> tmp6;
            Contract(1.0,env.gt(row)[col],shape(0),LO,shape(0),0.0,tmp6);

            DArray<6> tmp6bis;
            Contract(1.0,tmp6,shape(3,0),QL,shape(0,1),0.0,tmp6bis);

            tmp6.clear();
            Contract(1.0,tmp6bis,shape(2,0),QL,shape(0,1),0.0,tmp6);

            tmp6bis.clear();
            Permute(tmp6,shape(0,3,5,1,2,4),tmp6bis);

            DArray<4> LO_env;
            Gemm(CblasNoTrans,CblasNoTrans,1.0,tmp6bis,env.gb(row - 1)[col],0.0,LO_env);

            //Right

            //construct the 'Right' eff environment
            tmp6.clear();
            Contract(1.0,env.gt(row)[Lx-1],shape(1),QR,shape(1),0.0,tmp6);

            DArray<8> tmp8;
            Contract(1.0,tmp6,shape(1),QR,shape(1),0.0,tmp8);

            DArray<8> tmp8bis;
            Contract(1.0,tmp8,shape(3,6),env.gb(row - 1)[Lx - 1],shape(1,2),0.0,tmp8bis);

            DArray<4> RO_env = tmp8bis.reshape_clear( shape( tmp8bis.shape(0),tmp8bis.shape(2),tmp8bis.shape(4),tmp8bis.shape(6) ) );

            //construct effective environment
            enum {i,j,k,l,m,n,o};

            N_eff.clear();
            Contract(1.0,LO_env,shape(i,j,k,l),RO_env,shape(i,m,n,l),0.0,N_eff,shape(j,m,k,n));

         }
         else{//middle

            //Left
            //first attach top to left unity
            DArray<6> tmp6;
            Contract(1.0,env.gt(row)[col],shape(0),LO,shape(0),0.0,tmp6);

            DArray<6> tmp6bis;
            Contract(1.0,tmp6,shape(3,0),QL,shape(0,1),0.0,tmp6bis);

            tmp6.clear();
            Contract(1.0,tmp6bis,shape(2,0),QL,shape(0,1),0.0,tmp6);

            tmp6bis.clear();
            Permute(tmp6,shape(0,3,5,1,2,4),tmp6bis);

            DArray<4> LO_env;
            Gemm(CblasNoTrans,CblasNoTrans,1.0,tmp6bis,env.gb(row - 1)[col],0.0,LO_env);

            //contract with right renormalized operator:
            tmp6.clear();
            Contract(1.0,env.gt(row)[col + 1],shape(3),RO,shape(0),0.0,tmp6);

            //to construct the R_environment
            tmp6bis.clear();
            Contract(1.0,tmp6,shape(1,3),QR,shape(1,3),0.0,tmp6bis);

            tmp6.clear();
            Contract(1.0,tmp6bis,shape(1,2),QR,shape(1,3),0.0,tmp6);

            tmp6bis.clear();
            Permute(tmp6,shape(0,2,4,3,5,1),tmp6bis);

            //construct the 'Right' eff environment
            DArray<4> RO_env;
            Gemm(CblasNoTrans,CblasTrans,1.0,tmp6bis,env.gb(row - 1)[col + 1],0.0,RO_env);

            //construct effective environment
            enum {i,j,k,l,m,n,o};

            N_eff.clear();
            Contract(1.0,LO_env,shape(i,j,k,l),RO_env,shape(i,m,n,l),0.0,N_eff,shape(j,m,k,n));

         }

      }
      else{//calc environment for vertical pairs

         if(row == 0){//left edge

            //Left
            //construct the 'Left' eff environment
            DArray<6> tmp6;
            Contract(1.0,env.gl(col - 1)[0],shape(1),QL,shape(1),0.0,tmp6);

            DArray<8> tmp8;
            Contract(1.0,tmp6,shape(1),QL,shape(1),0.0,tmp8);

            DArray<8> tmp8bis;
            Contract(1.0,tmp8,shape(3,6),env.gr(col)[0],shape(1,2),0.0,tmp8bis);

            DArray<4> LO_env = tmp8bis.reshape_clear( shape(tmp8bis.shape(1),tmp8bis.shape(3),tmp8bis.shape(5),tmp8bis.shape(7)) );

            //Right

            //contract with right renormalized operator:
            tmp6.clear();
            Contract(1.0,env.gl(col - 1)[1],shape(3),RO,shape(0),0.0,tmp6);

            //to construct the R_environment
            DArray<6> tmp6bis;
            Contract(1.0,tmp6,shape(1,3),QR,shape(1,3),0.0,tmp6bis);

            tmp6.clear();
            Contract(1.0,tmp6bis,shape(1,2),QR,shape(1,3),0.0,tmp6);

            tmp6bis.clear();
            Permute(tmp6,shape(0,2,4,3,5,1),tmp6bis);

            //construct the 'Right' eff environment
            DArray<4> RO_env;
            Gemm(CblasNoTrans,CblasTrans,1.0,tmp6bis,env.gr(col)[1],0.0,RO_env);

            //construct effective environment
            enum {i,j,k,l,m,n,o};

            N_eff.clear();
            Contract(1.0,LO_env,shape(i,j,k,l),RO_env,shape(i,m,n,l),0.0,N_eff,shape(j,m,k,n));

         }
         else if(row == Ly - 2){//right edge

            //Left

            //first attach top to left unity
            DArray<6> tmp6;
            Contract(1.0,env.gl(col - 1)[row],shape(0),LO,shape(0),0.0,tmp6);

            DArray<6> tmp6bis;
            Contract(1.0,tmp6,shape(3,0),QL,shape(0,1),0.0,tmp6bis);

            tmp6.clear();
            Contract(1.0,tmp6bis,shape(2,0),QL,shape(0,1),0.0,tmp6);

            tmp6bis.clear();
            Permute(tmp6,shape(0,3,5,1,2,4),tmp6bis);

            DArray<4> LO_env;
            Gemm(CblasNoTrans,CblasNoTrans,1.0,tmp6bis,env.gr(col)[row],0.0,LO_env);

            //Right

            //construct the 'Right' eff environment
            tmp6.clear();
            Contract(1.0,env.gl(col - 1)[Ly-1],shape(1),QR,shape(1),0.0,tmp6);

            DArray<8> tmp8;
            Contract(1.0,tmp6,shape(1),QR,shape(1),0.0,tmp8);

            DArray<8> tmp8bis;
            Contract(1.0,tmp8,shape(3,6),env.gr(col)[Ly - 1],shape(1,2),0.0,tmp8bis);

            DArray<4> RO_env = tmp8bis.reshape_clear( shape( tmp8bis.shape(0),tmp8bis.shape(2),tmp8bis.shape(4),tmp8bis.shape(6) ) );

            //construct effective environment
            enum {i,j,k,l,m,n,o};

            N_eff.clear();
            Contract(1.0,LO_env,shape(i,j,k,l),RO_env,shape(i,m,n,l),0.0,N_eff,shape(j,m,k,n));

         }
         else{//middle

            //Left
            //first attach top to left unity
            DArray<6> tmp6;
            Contract(1.0,env.gl(col - 1)[row],shape(0),LO,shape(0),0.0,tmp6);

            DArray<6> tmp6bis;
            Contract(1.0,tmp6,shape(3,0),QL,shape(0,1),0.0,tmp6bis);

            tmp6.clear();
            Contract(1.0,tmp6bis,shape(2,0),QL,shape(0,1),0.0,tmp6);

            tmp6bis.clear();
            Permute(tmp6,shape(0,3,5,1,2,4),tmp6bis);

            DArray<4> LO_env;
            Gemm(CblasNoTrans,CblasNoTrans,1.0,tmp6bis,env.gr(col)[row],0.0,LO_env);

            //contract with right renormalized operator:
            tmp6.clear();
            Contract(1.0,env.gl(col - 1)[row + 1],shape(3),RO,shape(0),0.0,tmp6);

            //to construct the R_environment
            tmp6bis.clear();
            Contract(1.0,tmp6,shape(1,3),QR,shape(1,3),0.0,tmp6bis);

            tmp6.clear();
            Contract(1.0,tmp6bis,shape(1,2),QR,shape(1,3),0.0,tmp6);

            tmp6bis.clear();
            Permute(tmp6,shape(0,2,4,3,5,1),tmp6bis);

            //construct the 'Right' eff environment
            DArray<4> RO_env;
            Gemm(CblasNoTrans,CblasTrans,1.0,tmp6bis,env.gr(col)[row + 1],0.0,RO_env);

            //construct effective environment
            enum {i,j,k,l,m,n,o};

            N_eff.clear();
            Contract(1.0,LO_env,shape(i,j,k,l),RO_env,shape(i,m,n,l),0.0,N_eff,shape(j,m,k,n));

         }

      }

   }

   /**
    * check if the effective norm, when contracted with left and right reduced tensors, effectively returns the norm
    */
   void check_N_eff(const DArray<4> &N_eff,const DArray<3> &a_L,const DArray<3> &a_R){

      enum {i,j,k,l,m,n,o};

      DArray<5> tmp5;
      Contract(1.0,N_eff,shape(i,j,k,l),a_R,shape(m,n,j),0.0,tmp5,shape(i,m,k,n,l));

      DArray<4> tmp4;
      Contract(1.0,tmp5,shape(i,j,k,n,l),a_R,shape(m,n,l),0.0,tmp4,shape(i,j,k,m));

      DArray<3> tmp3;
      Contract(1.0,tmp4,shape(i,j,k,l),a_L,shape(i,n,j),0.0,tmp3,shape(k,n,l));

   }

   /**
    * check the value of the cost function for the full update
    */
   double cost_function(const DArray<4> &N_eff,const DArray<4> &b,const DArray<3> &a_L,const DArray<3> &a_R){

      DArray<4> tmp4;
      Contract(1.0,a_L,shape(2),a_R,shape(0),0.0,tmp4);

      double val = -2.0 * Dot(tmp4,b);

      DArray<4> tmp4bis;
      Contract(1.0,tmp4,shape(1,2),tmp4,shape(1,2),0.0,tmp4bis);

      val += Dot(tmp4bis,N_eff);

      return val;

   }

   /**
    * propagate a step with a local action
    */
    void prop_local(PEPS<double> &peps){

       DArray<5> peps_eB;

       for(int row = 0;row < Ly;++row)
          for(int col = 0;col < Lx;++col){

             peps_eB.clear();
             Contract(1.0,trot.geB(),shape(1),peps(row,col),shape(2),0.0,peps_eB);

             Permute(peps_eB,shape(1,2,0,3,4),peps(row,col));

          }

    }

}
