#include <iostream>
#include <iomanip>
#include <fstream>
#include <complex>
#include <chrono>

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
    * @param peps the PEPS to be propagated
    * @param n_sweeps the number of sweeps performed for the solution of the linear problem
    */
   void step(PEPS<double> &peps,int n_sweeps){

      enum {i,j,k,l,m,n,o};

      // -------------------------------------------//
      // --- !!! (1) the bottom two rows (1) !!! ---// 
      // -------------------------------------------//

      //containers for the renormalized operators
      vector< DArray<5> > R(Lx - 1);

      DArray<5> L;

      DArray<4> QL;
      DArray<3> a_L;

      DArray<4> QR;
      DArray<3> a_R;

      DArray<4> N_eff_n;

      //construct the full top environment:
      env.calc('T',peps);

      //and the bottom row environment
      env.gb(0).fill('b',peps);

      //initialize the right operators for the bottom row
      contractions::init_ro('b',peps,R);

      //middle sites of the bottom row:
     // for(int col = 0;col < Lx - 1;++col){
      int col = 0;

      //--- (1) update the vertical pair on column 'col' ---
      construct_reduced_tensor('V','L',peps(0,col),QL,a_L);
      construct_reduced_tensor('V','R',peps(1,col),QR,a_R);

      calc_vertical_N_eff('b',col,L,QL,R[col],QR,N_eff_n);

      canonicalize_n(N_eff_n,a_L,QL,a_R,QR);

      //now do the update! Apply the gates!
      update_n(N_eff_n,a_L,a_R,n_sweeps);
      
      //and expand back to the full tensors:
      Contract(1.0,QL,shape(i,j,k,o),a_L,shape(o,m,n),0.0,peps(0,col),shape(j,n,m,i,k));
      Contract(1.0,a_R,shape(i,j,k),QR,shape(k,o,m,n),0.0,peps(1,col),shape(o,n,j,i,m));

      //--- (2) update the horizontal pair on column 'col'-'col+1' ---
      
      //construct the reduced tensors of the pair to propagate
      construct_reduced_tensor('H','L',peps(0,col),QL,a_L);
      construct_reduced_tensor('H','R',peps(0,col+1),QR,a_R);

      //calculate the effective environment N_eff
      calc_horizontal_N_eff('b',col,peps,L,QL,R[col + 1],QR,N_eff_n);

      //make environment close to unitary before the update
      canonicalize_n(N_eff_n,a_L,QL,a_R,QR);

      //now do the update! Apply the gates!
      update_n(N_eff_n,a_L,a_R,n_sweeps);
     /*
      //and expand back to the full tensors
      Contract(1.0,QL,shape(i,j,k,o),a_L,shape(o,m,n),0.0,peps(0,col),shape(i,j,m,k,n));
      Contract(1.0,a_R,shape(i,j,k),QR,shape(k,o,m,n),0.0,peps(0,col+1),shape(i,o,j,m,n));

      contractions::update_L('b',col,peps,L);

      //}

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

      */
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

   void calc_vertical_N_eff(char option,int col,const DArray<5> &L,const DArray<4> &QL,const DArray<5> &R,const DArray<4> &QR,DArray<4> &N_eff_n){

      enum {i,j,k,l,m,n,o,p};

      if(option == 'b'){//bottom two rows

         if(col == 0){

            DArray<7> tmp7;
            Gemm(CblasNoTrans,CblasNoTrans,1.0,env.gt(0)[0],R,0.0,tmp7);

            DArray<7> tmp7bis;
            Permute(tmp7,shape(0,3,1,2,4,5,6),tmp7bis);

            //add QR to tmp7bis
            DArray<5> tmp5;
            Gemm(CblasNoTrans,CblasNoTrans,1.0,QR,tmp7bis,0.0,tmp5);

            DArray<5> tmp5bis;
            Permute(tmp5,shape(2,1,0,3,4),tmp5bis);

            //and another
            int M = QR.shape(0);
            int N = tmp5bis.shape(2) * tmp5bis.shape(3) * tmp5bis.shape(4);
            int K = tmp5bis.shape(0) * tmp5bis.shape(1);

            DArray<4> tmp4(QR.shape(0),QR.shape(0),D,D);
            blas::gemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, M, N, K, 1.0,QR.data(),K,tmp5bis.data(),N,0.0,tmp4.data(),N);

            DArray<6> tmp6;
            Contract(1.0,tmp4,shape(3),QL,shape(2),0.0,tmp6);

            N_eff_n.clear();
            Contract(1.0,QL,shape(i,j,k,l),tmp6,shape(m,n,k,i,j,o),0.0,N_eff_n,shape(l,n,o,m));

         }
         else{

         }

      }
      else{//top two rows

      }

   }

   void calc_horizontal_N_eff(char option,int col,const PEPS<double> &peps,const DArray<5> &L,const DArray<4> &QL,const DArray<5> &R,const DArray<4> &QR,DArray<4> &N_eff_n){

      enum {i,j,k,l,m,n,o,p};

      if(option == 'b'){//bottom two rows

         if(col == 0){

            //paste bottom QR to R
            DArray<7> tmp7;
            Gemm(CblasNoTrans,CblasTrans,1.0,QR,R,0.0,tmp7);

            //paste other (top) QR on
            int M = QR.shape(0) * QR.shape(1);
            int N = tmp7.shape(0) * tmp7.shape(1) * tmp7.shape(3) * tmp7.shape(4) * tmp7.shape(5);
            int K = QR.shape(3);

            DArray<7> tmp7bis(d*D,D,d*D,D,tmp7.shape(3),D,D);
            blas::gemm(CblasRowMajor, CblasNoTrans, CblasTrans, M, N, K, 1.0,QR.data(),K,tmp7.data(),K,0.0,tmp7bis.data(),N);

            tmp7.clear();
            Permute(tmp7bis,shape(3,6,0,1,2,4,5),tmp7);

            //add the tensors on row = 1
            DArray<8> tmp8;
            Gemm(CblasNoTrans,CblasNoTrans,1.0,peps(1,1),tmp7,0.0,tmp8);

            DArray<8> tmp8bis;
            Permute(tmp8,shape(2,4,7,0,1,3,5,6),tmp8bis);

            tmp7.clear();
            Gemm(CblasNoTrans,CblasNoTrans,1.0,peps(1,1),tmp8bis,0.0,tmp7);

            Permute(tmp7,shape(1,3,6,0,2,4,5),tmp7bis);

            //finally get 'right hand side': R_eff
            DArray<5> R_eff;
            Gemm(CblasNoTrans,CblasNoTrans,1.0,env.gt(0)[1],tmp7bis,0.0,R_eff);

            //left side, connect top to row 1 peps
            DArray<5> tmp5;
            Gemm(CblasTrans,CblasNoTrans,1.0,env.gt(0)[0],peps(1,0),0.0,tmp5);

            DArray<5> tmp5bis;
            Permute(tmp5,shape(1,3,4,0,2),tmp5bis);
            
            //add another PEPS
            M = tmp5bis.shape(0) * tmp5bis.shape(1) * tmp5bis.shape(2);
            N = D*D;
            K = tmp5bis.shape(3) * tmp5bis.shape(4);

            tmp5.resize( shape(env.gt(0)[0].shape(3),D,D,D,D) );
            blas::gemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, M, N, K, 1.0,tmp5bis.data(),K,peps(1,0).data(),N,0.0,tmp5.data(),N);

            tmp5bis.clear();
            Permute(tmp5,shape(0,2,4,3,1),tmp5bis);

            //now contract with QL
            M = tmp5bis.shape(0) * tmp5bis.shape(1) * tmp5bis.shape(2) * tmp5bis.shape(3);
            N = D;
            K = D;

            tmp5.resize(env.gt(0)[0].shape(3),D,D,D,D);
            blas::gemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, M, N, K, 1.0,tmp5bis.data(),K,QL.data(),N,0.0,tmp5.data(),N);
 
            //one final permuation
            tmp5bis.clear();
            Permute(tmp5,shape(0,1,2,4,3),tmp5bis);
            
            //now contract with (bottom) QL
            M = tmp5bis.shape(0) * tmp5bis.shape(1) * tmp5bis.shape(2) * tmp5bis.shape(3);
            N = D;
            K = D;

            DArray<5> L_eff(env.gt(0)[0].shape(3),D,D,D,D);
            blas::gemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, M, N, K, 1.0,tmp5bis.data(),K,QL.data(),N,0.0,L_eff.data(),N);

            //now construct N_eff:
            DArray<4> tmp4(D,D,d*D,d*D);
            Gemm(CblasTrans,CblasNoTrans,1.0,L_eff,R_eff,0.0,tmp4);

            N_eff_n.clear();
            Permute(tmp4,shape(0,2,1,3),N_eff_n);
 
         }
         else{

         }

      }
      else{//top two rows

      }

   }


   /**
    * make the environment as 'canonical' as possible so that the svd for the pair update is as optimal as possible. For a nearest neighbour update!
    * Do this by contructing the best positive symmetric approximation to the environment and taking a QR and RQ of the resultant X^+ X = N_eff
    * @param X (XX^T) is the best positive approximation to the environment of the pair
    * @param a_L left reduced tensor, will be multiplied with the L of the environment
    * @param QL unitary part of the left tensor reduction, will be multiplied with the inverse of the L of the environment
    * @param a_R right reduced tensor, will be multiplied with the R of the environment
    * @param QR unitary part of the right tensor reduction, will be multiplied with the inverse of the R of the environment
    */
   void canonicalize_n(DArray<4> &N_eff_n,DArray<3> &a_L,DArray<4> &QL,DArray<3> &a_R,DArray<4> &QR){

      int matdim = N_eff_n.shape(0)*N_eff_n.shape(1);

      //symmetrize
      for(int i = 0;i < matdim;++i)
         for(int j = i + 1;j < matdim;++j){

            N_eff_n.data()[i*matdim + j] = 0.5 * (N_eff_n.data()[i*matdim + j]  + N_eff_n.data()[j*matdim + i]);
            N_eff_n.data()[j*matdim + i] = N_eff_n.data()[i*matdim + j];

         }

      DArray<1> eig(matdim);

      //diagonalize
      lapack::syev(CblasRowMajor, 'V','U', matdim, N_eff_n.data(), matdim, eig.data());

      DArray<3> X(N_eff_n.shape(0),matdim,N_eff_n.shape(1));

      //get the square root of the positive approximant:
      for(int iL = 0;iL < N_eff_n.shape(0);++iL)
         for(int iR = 0;iR < N_eff_n.shape(1);++iR)
            for(int kL = 0;kL < N_eff_n.shape(0);++kL)
               for(int kR = 0;kR < N_eff_n.shape(1);++kR){

                  if(eig(kL*N_eff_n.shape(1) + kR) > 0.0)
                     X(iL,kL*N_eff_n.shape(1) + kR,iR) = sqrt( eig(kL*N_eff_n.shape(1) + kR) ) * N_eff_n(iL,iR,kL,kR);

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

      X.clear();
      Contract(1.0,tmp2,shape(1),X_copy,shape(0),0.0,X);

      //finally get the best positive approximation to N_eff
      Contract(1.0,X,shape(1),X,shape(1),0.0,N_eff_n);

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
    * update a nearest-neigbour tensor pair by applying a trotter gate over their mutual bond
    * after which a svd is performed over the bond and the dimensions are set back to D
    */
   void update_n(const DArray<4> &N_eff_n,DArray<3> &a_L,DArray<3> &a_R,int n_iter){

      enum {i,j,k,l,m,n};

      //left
      DArray<4> tmp4;
      Contract(1.0,a_L,shape(i,j,k),trot.gLO_n(),shape(n,m,j),0.0,tmp4,shape(i,n,m,k));

      a_L = tmp4.reshape_clear(shape(a_L.shape(0),a_L.shape(1),a_L.shape(2)*trot.gLO_n().shape(1)));

      //right
      tmp4.clear();
      Contract(1.0,trot.gRO_n(),shape(i,j,k),a_R,shape(n,k,m),0.0,tmp4,shape(j,n,i,m));

      a_R = tmp4.reshape_clear(shape(a_R.shape(0)*trot.gRO_n().shape(1),a_R.shape(1),a_R.shape(2)));

      //now create 'two-site' object
      tmp4.clear();
      Contract(1.0,a_L,shape(2),a_R,shape(0),0.0,tmp4);

      //create 'b' object of linear system: N_eff A = b
      DArray<4> b;
      Contract(1.0,N_eff_n,shape(i,j,k,l),tmp4,shape(k,m,n,l),0.0,b,shape(i,m,n,j));

      //svd the fucker
      DArray<1> S;
      Gesvd ('S','S', tmp4, S,a_L,a_R,D);

      //take the square root of the sv's
      for(int i = 0;i < S.size();++i)
         S(i) = sqrt(S(i));

      //and multiply it left and right to the tensors
      Dimm(S,a_R);
      Dimm(a_L,S);

      //start sweeping
      int iter = 0;

      while(iter < n_iter){

         //construct right hand side 
         DArray<3> tmp3;
         Contract(1.0,b,shape(i,j,k,l),a_R,shape(m,k,l),0.0,tmp3,shape(i,m,j));

         //construct effective environment for left site
         DArray<5> tmp5;
         Contract(1.0,N_eff_n,shape(3),a_R,shape(2),0.0,tmp5);

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
         Contract(1.0,a_L,shape(0),N_eff_n,shape(0),0.0,tmp5);

         tmp4.clear();
         Contract(1.0,tmp5,shape(i,j,k,l,m),a_L,shape(l,i,n),0.0,tmp4,shape(j,k,n,m));

         solve(tmp4,tmp3);

         Permute(tmp3,shape(0,2,1),a_R);

         ++iter;

      }

      //When converged, put both objects on equal footing
      tmp4.clear();
      Contract(1.0,a_L,shape(2),a_R,shape(0),0.0,tmp4);

      Gesvd ('S','S', tmp4, S,a_L,a_R,D);

      //take the square root of the sv's
      for(int i = 0;i < S.size();++i)
         S(i) = sqrt(S(i));

      //and multiply it left and right to the tensors
      Dimm(S,a_R);
      Dimm(a_L,S);

   }

   /**
    * check the value of the cost function for the nearest neigbour update!
    */
   double cost_function_n(const DArray<4> &N_eff_n,const DArray<4> &b,const DArray<3> &a_L,const DArray<3> &a_R){

      DArray<4> tmp4;
      Contract(1.0,a_L,shape(2),a_R,shape(0),0.0,tmp4);

      double val = -2.0 * Dot(tmp4,b);

      DArray<4> tmp4bis;
      Contract(1.0,tmp4,shape(1,2),tmp4,shape(1,2),0.0,tmp4bis);

      val += Dot(tmp4bis,N_eff_n);

      return val;

   }

   /** 
    * wrapper function solve positive symmetric linear system: N_eff * x = b
    * @param N_eff input matrix
    * @param b right hand side input and x output
    */
   void solve(DArray<4> &N_eff,DArray<3> &b){

      int n = N_eff.shape(0) * N_eff.shape(1);

      lapack::potrf(CblasRowMajor,'U',n, N_eff.data(), n);

      lapack::potrs(CblasRowMajor,'U',n,d, N_eff.data(), n,b.data(),d);

   }

}

