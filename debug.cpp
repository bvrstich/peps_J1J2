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
using namespace propagate;
using namespace global;

namespace debug {

   /**
    * evaluate the cost function of the linear system for two vertically connected PEPS: -2 <\Psi|\Psi'> + <\Psi'|\Psi'> where \Psi full PEPS 
    * with operator and \Psi is old PEPS.
    * for top or bottom rows, i.e. with L and R of order 5 and intermediates of order 7
    * @param dir vertical, horizontal,diagonal lurd or diagonal ldru update
    * @param row , the row index of the bottom site
    * @param col column index of the vertical column
    * @param peps, full PEPS object before update
    * @param L Left environment contraction
    * @param R Right environment contraction
    * @param LI7 left intermediate object
    * @param RI7 left intermediate object
    * @param b_L left intermediate object
    * @param b_R right intermediate object
    */
   template<>
      double cost_function(const PROP_DIR &dir,int row,int col,PEPS<double> &peps,const DArray<6> &lop,const DArray<6> &rop,

            const DArray<5> &L,const DArray<5> &R,const DArray<7> &LI7,const DArray<7> &RI7,const DArray<7> &b_L,const DArray<7> &b_R){ 

         if(dir == VERTICAL){

            if(row == 0){

               if(col == 0){

                  // --- (1) calculate overlap of approximated state:

                  //paste bottom peps to right intermediate
                  DArray<10> tmp10;
                  Gemm(CblasNoTrans,CblasTrans,1.0,RI7,peps(row,col),0.0,tmp10);

                  DArray<7> tmp7 = tmp10.reshape_clear( shape(D,D,D,D,D,D,d) );

                  //another bottom peps to this one
                  DArray<8> tmp8;
                  Contract(1.0,tmp7,shape(6,4),peps(row,col),shape(2,4),0.0,tmp8);

                  DArray<6> tmp6 = tmp8.reshape_clear( shape(D,D,D,D,D,D) );

                  //add upper peps
                  DArray<5> tmp5;
                  Contract(1.0,tmp6,shape(0,5,2),peps(row+1,col),shape(1,3,4),0.0,tmp5);

                  DArray<5> tmp5bis;
                  Permute(tmp5,shape(3,0,4,2,1),tmp5bis);

                  double val = Dot(tmp5bis,peps(row+1,col));

                  // --- (2) calculate 'b' part of overlap

                  //right hand side: add left operator to tmp7
                  DArray<9> tmp9;
                  Contract(1.0,tmp7,shape(6,4),lop,shape(2,5),0.0,tmp9);

                  //remove the dimension-one legs
                  tmp7 = tmp9.reshape_clear( shape(D,D,D,D,D,D,global::trot.gLO_n().shape(1)) );

                  tmp5.clear();
                  Contract(1.0,tmp7,shape(0,6,5,2),rop,shape(1,3,4,5),0.0,tmp5);

                  tmp5bis.clear();
                  Permute(tmp5,shape(3,0,4,2,1),tmp5bis);

                  val -= 2.0 * Dot(tmp5bis,peps(row+1,col));

                  return val;

               }
               else if(col < Lx - 1){//col != 0

                  // --- (1) calculate overlap of approximated state:

                  //add upper peps to LI7
                  DArray<8> tmp8;
                  Contract(1.0,LI7,shape(1,5),peps(row+1,col),shape(0,1),0.0,tmp8);

                  //and another
                  DArray<7> tmp7;
                  Contract(1.0,tmp8,shape(0,3,5),peps(row+1,col),shape(0,1,2),0.0,tmp7);

                  DArray<8> tmp8bis;
                  Contract(1.0,tmp7,shape(0,5),peps(row,col),shape(0,1),0.0,tmp8bis);

                  DArray<5> tmp5;
                  Contract(1.0,tmp8bis,shape(0,2,5,6),peps(row,col),shape(0,1,2,3),0.0,tmp5);

                  DArray<5> tmp5bis;
                  Permute(tmp5,shape(0,2,1,3,4),tmp5bis);

                  double val = Dot(tmp5bis,R);

                  // --- (2) calculate 'b' part of overlap

                  //add right operator to tmp8
                  tmp8bis.clear();
                  Contract(1.0,tmp8,shape(0,3,5),rop,shape(0,1,2),0.0,tmp8bis);

                  //then add left operator
                  tmp8.clear();
                  Contract(1.0,tmp8bis,shape(0,6,5),lop,shape(0,1,3),0.0,tmp8);

                  //finally add lop
                  tmp5.clear();
                  Contract(1.0,tmp8,shape(0,2,5,6),peps(row,col),shape(0,1,2,3),0.0,tmp5);

                  tmp5bis.clear();
                  Permute(tmp5,shape(0,2,1,3,4),tmp5bis);

                  val -= 2.0 * Dot(tmp5bis,R);

                  return val;

               }
               else{//col == Lx - 1

                  DArray<8> tmp8;
                  Contract(1.0,LI7,shape(1,5),peps(row+1,col),shape(0,1),0.0,tmp8);

                  //and another
                  DArray<7> tmp7;
                  Contract(1.0,tmp8,shape(0,3,5),peps(row+1,col),shape(0,1,2),0.0,tmp7);

                  DArray<8> tmp8bis;
                  Contract(1.0,tmp7,shape(0,5),peps(row,col),shape(0,1),0.0,tmp8bis);

                  DArray<5> tmp5 = tmp8bis.reshape_clear( shape(D,D,d,1,1) );

                  double val =  Dot(tmp5,peps(row,col));

                  // --- (2) calculate 'b' part of overlap

                  //add right operator to tmp8
                  tmp8bis.clear();
                  Contract(1.0,tmp8,shape(0,3,5),rop,shape(0,1,2),0.0,tmp8bis);

                  //then add left operator
                  tmp8.clear();
                  Contract(1.0,tmp8bis,shape(0,6,5),lop,shape(0,1,3),0.0,tmp8);

                  //finally add lop
                  tmp5 = tmp8.reshape_clear( shape(D,D,d,1,1) );

                  val -= 2.0 * Dot(tmp5,peps(row,col));

                  return val;

               }

            }
            else{//row == Lx - 2 

               if(col == 0){

                  // (1) construct N_eff

                  //paste bottom peps to right intermediate
                  DArray<8> tmp8;
                  Contract(1.0,peps(row,col),shape(3,4),RI7,shape(2,6),0.0,tmp8);

                  //and another bottom peps to tmp8
                  DArray<7> tmp7;
                  Contract(1.0,peps(row,col),shape(2,3,4),tmp8,shape(2,4,7),0.0,tmp7);

                  //upper peps
                  DArray<8> tmp8bis;
                  Contract(1.0,peps(row+1,col),shape(3,4),tmp7,shape(3,6),0.0,tmp8bis);

                  DArray<5> tmp5 = tmp8bis.reshape_clear( shape(1,1,d,D,D) );

                  double val = Dot(tmp5,peps(row+1,col));

                  //(2) right hand side:

                  //add left operator to tmp8
                  DArray<6> tmp6;
                  Contract(1.0,lop,shape(0,2,4,5),tmp8,shape(0,2,4,7),0.0,tmp6);

                  //add right operator
                  DArray<6> tmp6bis;
                  Contract(1.0,rop,shape(3,4,5),tmp6,shape(1,0,4),0.0,tmp6bis);

                  val -= 2.0 * blas::dot(tmp6bis.size(), tmp6bis.data(), 1, peps(row+1,col).data(), 1);

                  return val;

               }
               else if(col < Lx - 1){

                  // (1) construct N_eff

                  //paste bottom peps to right intermediate
                  DArray<8> tmp8;
                  Contract(1.0,peps(row,col),shape(3,4),RI7,shape(2,6),0.0,tmp8);

                  //and another bottom peps to tmp8
                  DArray<7> tmp7;
                  Contract(1.0,peps(row,col),shape(2,3,4),tmp8,shape(2,4,7),0.0,tmp7);

                  //add top peps
                  DArray<8> tmp8bis;
                  Contract(1.0,peps(row+1,col),shape(3,4),tmp7,shape(3,6),0.0,tmp8bis);

                  //final top peps
                  DArray<5> tmp5;
                  Contract(1.0,peps(row+1,col),shape(1,2,3,4),tmp8bis,shape(1,2,4,7),0.0,tmp5);

                  //contact with left hand side
                  double val = Dot(tmp5,L);

                  //(2) right hand side:

                  //add left operator to tmp8
                  tmp8bis.clear();
                  Contract(1.0,lop,shape(2,4,5),tmp8,shape(2,4,7),0.0,tmp8bis);

                  //add right operator
                  tmp8.clear();
                  Contract(1.0,rop,shape(3,4,5),tmp8bis,shape(2,1,6),0.0,tmp8);

                  //add top peps
                  tmp5.clear();
                  Contract(1.0,peps(row+1,col),shape(1,2,3,4),tmp8,shape(1,2,5,7),0.0,tmp5);

                  DArray<5> tmp5bis;
                  Permute(tmp5,shape(1,0,2,3,4),tmp5bis);

                  //and contract with left hand side
                  val -= 2.0 * Dot(tmp5bis,L);

                  return val;

               }
               else{//col == Lx - 1

                  // (1) construct N_eff

                  //paste bottom peps to left intermediate
                  DArray<6> tmp6;
                  Contract(1.0,peps(row,col),shape(0,3,4),LI7,shape(3,5,6),0.0,tmp6);

                  //and another
                  DArray<5> tmp5;
                  Contract(1.0,peps(row,col),shape(0,2,3),tmp6,shape(4,1,5),0.0,tmp5);

                  //add top peps to it
                  DArray<4> tmp4;
                  Contract(1.0,peps(row+1,col),shape(0,3,4),tmp5,shape(4,2,1),0.0,tmp4);

                  DArray<4> tmp4bis;
                  Permute(tmp4,shape(3,0,1,2),tmp4bis);

                  tmp5 = tmp4bis.reshape_clear( shape(D,1,d,D,1));

                  double val = Dot(tmp5,peps(row+1,col));

                  //(2) right hand side:

                  //add left operator to tmp6
                  DArray<6> tmp6bis;
                  Contract(1.0,lop,shape(0,2,4),tmp6,shape(4,1,5),0.0,tmp6bis);

                  //add right operator
                  tmp4.clear();
                  Contract(1.0,rop,shape(0,3,4,5),tmp6bis,shape(4,1,0,2),0.0,tmp4);

                  tmp4bis.clear();
                  Permute(tmp4,shape(3,0,1,2),tmp4bis);

                  tmp5 = tmp4bis.reshape_clear( shape(D,1,d,D,1));

                  val -=  2.0 * Dot(tmp5,peps(row+1,col));

                  return val;

               }

            }

         }//end vertical
         else if(dir == HORIZONTAL){

            if(row == 0){

               // --- (1) calculate overlap of approximated state:

               //add peps to right side
               DArray<8> tmp8;
               Contract(1.0,RI7,shape(4,6),peps(row,col+1),shape(1,4),0.0,tmp8);

               //add second peps
               DArray<5> tmp5;
               Contract(1.0,tmp8,shape(3,6,7,4),peps(row,col+1),shape(1,2,3,4),0.0,tmp5);

               //add peps to left side
               DArray<8> tmp8bis;
               Contract(1.0,LI7,shape(6,4),peps(row,col),shape(0,1),0.0,tmp8bis);

               DArray<5> tmp5bis;
               Contract(1.0,tmp8bis,shape(4,3,5,6),peps(row,col),shape(0,1,2,3),0.0,tmp5bis);

               double val = Dot(tmp5bis,tmp5);

               // --- (2) calculate 'b' part of overlap

               //add right operator to tmp8
               DArray<6> tmp6;
               Contract(1.0,tmp8,shape(3,6,7,4),rop,shape(1,2,4,5),0.0,tmp6);

               //attach LI7 to right side
               DArray<7> tmp7;
               Gemm(CblasTrans,CblasNoTrans,1.0,LI7,tmp6,0.0,tmp7);

               //now paste left operator in
               tmp5.clear();
               Contract(1.0,tmp7,shape(2,0,6,5),lop,shape(0,1,3,5),0.0,tmp5);

               //and contract with peps(row,col)
               tmp5bis.clear();
               Permute(tmp5,shape(1,0,3,4,2),tmp5bis);

               val -= 2.0 * Dot(tmp5bis,peps(row,col));

               return val;

            }
            else if(row == Ly - 2){

               //add left peps to LI7
               DArray<8> tmp8_left;
               Contract(1.0,LI7,shape(5,3),peps(row,col),shape(0,1),0.0,tmp8_left);

               //and another peps
               DArray<7> tmp7;
               Contract(1.0,tmp8_left,shape(3,2,5),peps(row,col),shape(0,1,2),0.0,tmp7);

               //add bottom environment
               DArray<5> tmp5_left;
               Contract(1.0,tmp7,shape(2,5,3),env.gb(Ly-3)[col],shape(0,1,2),0.0,tmp5_left);

               //add right peps to intermediate
               DArray<8> tmp8_right;
               Contract(1.0,RI7,shape(3,5),peps(row,col+1),shape(1,4),0.0,tmp8_right);

               //add second peps
               tmp7.clear();
               Contract(1.0,tmp8_right,shape(2,6,3),peps(row,col+1),shape(1,2,4),0.0,tmp7);

               //add bottom environment
               DArray<5> tmp5;
               Contract(1.0,tmp7,shape(6,4,2),env.gb(Ly-3)[col+1],shape(1,2,3),0.0,tmp5);

               double val = Dot(tmp5,tmp5_left);

               //right hand side

               //add left operator
               DArray<8> tmp8bis;
               Contract(1.0,tmp8_left,shape(3,2,5),lop,shape(0,1,2),0.0,tmp8bis);

               //add bottom environment
               DArray<6> tmp6_left;
               Contract(1.0,tmp8bis,shape(2,6,3),env.gb(Ly-3)[col],shape(0,1,2),0.0,tmp6_left);

               //add right operator to tmp8
               tmp8bis.clear();
               Contract(1.0,tmp8_right,shape(2,6,3),rop,shape(1,2,5),0.0,tmp8bis);

               //add bottom environment
               DArray<6> tmp6_right;
               Contract(1.0,env.gb(Ly-3)[col+1],shape(1,2,3),tmp8bis,shape(7,4,2),0.0,tmp6_right);

               DArray<6> tmp6bis;
               Permute(tmp6_right,shape(1,2,3,5,4,0),tmp6bis);

               val -= 2.0 * Dot(tmp6_left,tmp6bis);

               return val;

            }
            else{//row == Ly - 1

               //add right peps to intermediate
               DArray<8> tmp8_right;
               Contract(1.0,peps(row,col+1),shape(3,4),RI7,shape(3,1),0.0,tmp8_right);

               //add second peps
               DArray<5> tmp5_right;
               Contract(1.0,peps(row,col+1),shape(1,2,3,4),tmp8_right,shape(1,2,4,3),0.0,tmp5_right);

               //add left to intermediate
               DArray<8> tmp8_left;
               Contract(1.0,LI7,shape(1,3),peps(row,col),shape(0,3),0.0,tmp8_left);

               //and another
               DArray<5> tmp5_left;
               Contract(1.0,tmp8_left,shape(0,5,6,1),peps(row,col),shape(0,1,2,3),0.0,tmp5_left);

               DArray<5> tmp5bis;
               Permute(tmp5_left,shape(4,3,0,1,2),tmp5bis);

               //norm
               double val = Dot(tmp5bis,tmp5_right);

               //add right operator to tmp8_right
               DArray<6> tmp6_right;
               Contract(1.0,rop,shape(1,2,4,5),tmp8_right,shape(1,2,4,3),0.0,tmp6_right);

               //add left operator
               DArray<6> tmp6_left;
               Contract(1.0,tmp8_left,shape(0,5,6,1),lop,shape(0,1,2,4),0.0,tmp6_left);

               DArray<6> tmp6bis;
               Permute(tmp6_left,shape(5,4,3,0,1,2),tmp6bis);

               val -= 2.0 * Dot(tmp6_right,tmp6bis);

               return val;

            }

         }//end horizontal
         else if(dir == DIAGONAL_LURD){

            if(row == 0){

               // --- (1) calculate overlap of approximated state:

               //right side
               DArray<8> tmp8_right;
               Contract(1.0,RI7,shape(4,6),peps(row,col+1),shape(1,4),0.0,tmp8_right);

               //and again
               DArray<5> tmp5_right;
               Contract(1.0,tmp8_right,shape(3,6,7,4),peps(row,col+1),shape(1,2,3,4),0.0,tmp5_right);

               //left side
               DArray<8> tmp8_left;
               Contract(1.0,LI7,shape(2,4),peps(row+1,col),shape(0,3),0.0,tmp8_left);

               //add another
               DArray<7> tmp7;
               Contract(1.0,tmp8_left,shape(1,6,2),peps(row+1,col),shape(0,2,3),0.0,tmp7);

               //add top environment to intermediate
               DArray<5> tmp5_left;
               Contract(1.0,tmp7,shape(0,5,3),env.gt(row)[col],shape(0,1,2),0.0,tmp5_left);

               DArray<5> tmp5;
               Permute(tmp5_left,shape(4,3,2,1,0),tmp5);

               double val = Dot(tmp5,tmp5_right);

               // --- (2) calculate 'b' part of overlap

               //add right operator to tmp8_right
               DArray<6> tmp6_right;
               Contract(1.0,tmp8_right,shape(3,6,7,4),rop,shape(1,2,4,5),0.0,tmp6_right);

               //add top peps to b_L
               tmp8_left.clear();
               Contract(1.0,b_L,shape(2,4),peps(row+1,col),shape(0,3),0.0,tmp8_left);

               DArray<8> tmp8bis;
               Contract(1.0,tmp8_left,shape(1,6,2),lop,shape(0,2,4),0.0,tmp8bis);

               //add top environment to intermediate
               DArray<6> tmp6_left;
               Contract(1.0,tmp8bis,shape(0,5,3),env.gt(row)[col],shape(0,1,2),0.0,tmp6_left);

               DArray<6> tmp6bis;
               Permute(tmp6_left,shape(5,4,2,1,0,3),tmp6bis);

               val -= 2.0 * Dot(tmp6bis,tmp6_right);

               return val;

            }
            else{//row == Ly - 2

               // --- (1) calculate overlap of approximated state:

               //add right peps to intermediate
               DArray<8> tmp8_right;
               Contract(1.0,RI7,shape(3,5),peps(row,col+1),shape(1,4),0.0,tmp8_right);

               //add second peps
               DArray<7> tmp7;
               Contract(1.0,tmp8_right,shape(2,6,3),peps(row,col+1),shape(1,2,4),0.0,tmp7);

               //add bottom environment
               DArray<5> tmp5_right;
               Contract(1.0,env.gb(Ly-3)[col+1],shape(1,2,3),tmp7,shape(6,4,2),0.0,tmp5_right);

               //add upper-left peps to intermediate left
               DArray<8> tmp8_left;
               Contract(1.0,LI7,shape(1,3),peps(row+1,col),shape(0,3),0.0,tmp8_left);

               //und again
               DArray<5> tmp5;
               Contract(1.0,tmp8_left,shape(0,5,6,1),peps(row+1,col),shape(0,1,2,3),0.0,tmp5);

               DArray<5> tmp5_left;
               Permute(tmp5,shape(2,4,3,1,0),tmp5_left);

               double val = Dot(tmp5_left,tmp5_right);

               // --- (2) operator part of cost function

               //add right operator to tmp8_right
               DArray<8> tmp8bis;
               Contract(1.0,tmp8_right,shape(2,6,3),rop,shape(1,2,5),0.0,tmp8bis);

               //add bottom environment
               DArray<6> tmp6_right;
               Contract(1.0,env.gb(row-1)[col+1],shape(1,2,3),tmp8bis,shape(7,4,2),0.0,tmp6_right);

               //add upper-left peps to b_L
               tmp8_left.clear();
               Contract(1.0,b_L,shape(1,3),peps(row+1,col),shape(0,3),0.0,tmp8_left);

               //next add left operator to tmp8
               DArray<6> tmp6;
               Contract(1.0,tmp8_left,shape(0,5,6,1),lop,shape(0,1,2,4),0.0,tmp6);

               DArray<6> tmp6_left;
               Permute(tmp6,shape(2,5,3,1,0,4),tmp6_left);

               val -= 2.0 * Dot(tmp6_left,tmp6_right);

               return val;

            }


         }
         else{//diagonal ldru

            if(row == 0){

               // --- (1) calculate overlap of approximated state:

               //RIGHT
               DArray<8> tmp8_right;
               Contract(1.0,peps(row+1,col+1),shape(3,4),RI7,shape(4,2),0.0,tmp8_right);

               //and again
               DArray<7> tmp7;
               Contract(1.0,peps(row+1,col+1),shape(2,3,4),tmp8_right,shape(2,5,4),0.0,tmp7);

               //and add top environment to intermediate
               DArray<5> tmp5_right;
               Contract(1.0,env.gt(row)[col+1],shape(1,2,3),tmp7,shape(1,3,4),0.0,tmp5_right);

               //LEFT
               DArray<8> tmp8_left;
               Contract(1.0,LI7,shape(6,4),peps(row,col),shape(0,1),0.0,tmp8_left);

               //and again
               DArray<5> tmp5;
               Contract(1.0,tmp8_left,shape(4,3,5,6),peps(row,col),shape(0,1,2,3),0.0,tmp5);

               DArray<5> tmp5_left;
               Permute(tmp5,shape(0,1,2,4,3),tmp5_left);

               double val = Dot(tmp5_left,tmp5_right);

               // --- (2) calculate 'b' part of overlap

               //RIGHT

               //add top right peps to b_R
               tmp8_right.clear();
               Contract(1.0,peps(row+1,col+1),shape(3,4),b_R,shape(4,2),0.0,tmp8_right);

               //add right operator to tmp8
               DArray<8> tmp8bis;
               Contract(1.0,rop,shape(2,4,5),tmp8_right,shape(2,5,4),0.0,tmp8bis);

               //add top environment to intermediate
               DArray<6> tmp6_right;
               Contract(1.0,env.gt(row)[col+1],shape(1,2,3),tmp8bis,shape(1,4,5),0.0,tmp6_right);

               //LEFT

               //add left operator to tmp8_left
               DArray<6> tmp6;
               Contract(1.0,tmp8_left,shape(4,3,5,6),lop,shape(0,1,2,4),0.0,tmp6);

               DArray<6> tmp6_left;
               Permute(tmp6,shape(0,1,4,2,5,3),tmp6_left);

               val -= 2.0 * Dot(tmp6_left,tmp6_right);

               return val;

            }
            else{//row == Ly - 2

               // --- (1) calculate overlap of approximated state:

               //right

               DArray<8> tmp8_right;
               Contract(1.0,peps(row+1,col+1),shape(3,4),RI7,shape(3,1),0.0,tmp8_right);

               //and again
               DArray<5> tmp5_right;
               Contract(1.0,peps(row+1,col+1),shape(1,2,3,4),tmp8_right,shape(1,2,4,3),0.0,tmp5_right);

               //left

               //add bottom peps to LI7
               DArray<8> tmp8_left;
               Contract(1.0,LI7,shape(5,3),peps(row,col),shape(0,1),0.0,tmp8_left);

               //and again
               DArray<7> tmp7;
               Contract(1.0,tmp8_left,shape(3,2,5),peps(row,col),shape(0,1,2),0.0,tmp7);

               //now add bottom environment
               DArray<5> tmp5;
               Contract(1.0,tmp7,shape(2,5,3),env.gb(row-1)[col],shape(0,1,2),0.0,tmp5);

               DArray<5> tmp5_left;
               Permute(tmp5,shape(0,1,3,2,4),tmp5_left);

               double val = Dot(tmp5_left,tmp5_right);

               // --- (2) calculate 'b' part of overlap

               //RIGHT
               tmp8_right.clear();
               Contract(1.0,peps(row+1,col+1),shape(3,4),b_R,shape(3,1),0.0,tmp8_right);

               //and add right operator
               DArray<6> tmp6_right;
               Contract(1.0,rop,shape(1,2,4,5),tmp8_right,shape(1,2,4,3),0.0,tmp6_right);

               //LEFT

               //add left operator to intermediate
               DArray<8> tmp8bis;
               Contract(1.0,tmp8_left,shape(3,2,5),lop,shape(0,1,2),0.0,tmp8bis);

               //now add bottom environment
               DArray<6> tmp6;
               Contract(1.0,tmp8bis,shape(2,6,3),env.gb(row-1)[col],shape(0,1,2),0.0,tmp6);

               DArray<6> tmp6_left;
               Permute(tmp6,shape(0,3,1,4,2,5),tmp6_left);

               val -= 2.0 * Dot(tmp6_left,tmp6_right);

               return val;

            }

         }

      }

   /**
    * evaluate the cost function of the linear system for two vertically connected PEPS: -2 <\Psi|\Psi'> + <\Psi'|\Psi'> where \Psi full
    * PEPS with operator and \Psi is old PEPS
    * for middle rows, i.e. with R and L order 6 and intermediates of order 8
    * @param dir vertical, horizontal,diagonal lurd or diagonal ldru update
    * @param row , the row index of the bottom site
    * @param col column index of the vertical column
    * @param peps, full PEPS object before update
    * @param LO Left environment contraction
    * @param RO Right environment contraction
    * @param LI8 left intermediate object
    * @param RI8 right intermediate object
    */
   template<>
      double cost_function(const PROP_DIR &dir,int row,int col,PEPS<double> &peps,const DArray<6> &lop,const DArray<6> &rop,

            const DArray<6> &LO,const DArray<6> &RO,const DArray<8> &LI8,const DArray<8> &RI8,const DArray<8> &b_L,const DArray<8> &b_R){

         if(dir == VERTICAL){

            if(col == 0){

               //add bottom peps  to intermediate
               DArray<9> tmp9;
               Contract(1.0,peps(row,col),shape(3,4),RI8,shape(7,5),0.0,tmp9);

               //and another
               DArray<8> tmp8;
               Contract(1.0,peps(row,col),shape(2,3,4),tmp9,shape(2,8,7),0.0,tmp8);

               //and top one
               DArray<5> tmp5;
               Contract(1.0,peps(row+1,col),shape(0,1,3,4),tmp8,shape(0,4,1,6),0.0,tmp5);

               DArray<5> tmp5bis;
               Permute(tmp5,shape(1,3,0,2,4),tmp5bis);

               double val = Dot(tmp5bis,peps(row+1,col));

               // (2) right hand side

               //add left operator to intermediate
               DArray<9> tmp9bis;
               Contract(1.0,lop,shape(2,4,5),tmp9,shape(2,8,7),0.0,tmp9bis);

               //and right operator
               tmp5.clear();
               Contract(1.0,rop,shape(0,1,3,4,5),tmp9bis,shape(0,5,2,1,7),0.0,tmp5);

               tmp5bis.clear(); 
               Permute(tmp5,shape(1,3,0,2,4),tmp5bis);

               val -= 2.0 * Dot(tmp5bis,peps(row+1,col));

               return val;

            }
            else if(col < Lx -  1){//col != 0

               //add bottom peps  to intermediate right
               DArray<9> tmp9;
               Contract(1.0,peps(row,col),shape(3,4),RI8,shape(2,7),0.0,tmp9);

               //and another
               DArray<8> tmp8;
               Contract(1.0,peps(row,col),shape(2,3,4),tmp9,shape(2,4,8),0.0,tmp8);

               DArray<8> rn;
               Permute(tmp8,shape(5,6,7,1,3,0,2,4),rn);

               //add upper peps to LI8
               DArray<9> tmp9tris;
               Contract(1.0,LI8,shape(1,6),peps(row+1,col),shape(0,1),0.0,tmp9tris);

               //and another
               tmp8.clear();
               Contract(1.0,tmp9tris,shape(0,4,6),peps(row+1,col),shape(0,1,2),0.0,tmp8);

               DArray<8> ln;
               Permute(tmp8,shape(3,7,5,6,4,0,1,2),ln);

               double val = Dot(ln,rn);

               // (2) right hand side

               //add left operator to intermediate right
               DArray<9> tmp9bis;
               Contract(1.0,lop,shape(2,4,5),tmp9,shape(2,4,8),0.0,tmp9bis);

               //and right operator
               tmp9.clear();
               Contract(1.0,rop,shape(3,4,5),tmp9bis,shape(2,1,7),0.0,tmp9);

               //contract with left hand side
               DArray<5> tmp5;
               Contract(1.0,LI8,shape(7,5,0,2,3,4),tmp9,shape(7,1,0,3,4,6),0.0,tmp5);

               val -= 2.0 * Dot(tmp5,peps(row+1,col));

               return val;

            }
            else{//col = Lx - 1

               //add bottom peps  to intermediate
               DArray<9> tmp9;
               Contract(1.0,LI8,shape(3,5),peps(row,col),shape(0,3),0.0,tmp9);

               //and another
               DArray<8> tmp8;
               Contract(1.0,tmp9,shape(2,7,3),peps(row,col),shape(0,2,3),0.0,tmp8);

               //upper peps
               DArray<5> tmp5;
               Contract(1.0,tmp8,shape(0,2,6,7),peps(row+1,col),shape(0,1,3,4),0.0,tmp5);

               DArray<5> tmp5bis;
               Permute(tmp5,shape(0,1,4,2,3),tmp5bis);

               double val = Dot(tmp5bis,peps(row+1,col));

               //add left operator to intermediate
               DArray<9> tmp9bis;
               Contract(1.0,tmp9,shape(2,7,3),lop,shape(0,2,4),0.0,tmp9bis);

               //and right operator
               tmp5.clear();
               Contract(1.0,tmp9bis,shape(0,2,7,6,8),rop,shape(0,1,3,4,5),0.0,tmp5);

               tmp5bis.clear();;
               Permute(tmp5,shape(0,1,4,2,3),tmp5bis);

               val -= 2.0 * Dot(tmp5bis,peps(row+1,col));

               return val;

            }

         }//end vertical
         else if(dir == HORIZONTAL){

            //(1) overlap term

            //add right peps to intermediate
            DArray<9> tmp9;
            Contract(1.0,RI8,shape(4,6),peps(row,col+1),shape(1,4),0.0,tmp9);

            //add second peps
            DArray<8> tmp8;
            Contract(1.0,tmp9,shape(3,7,4),peps(row,col+1),shape(1,2,4),0.0,tmp8);

            //add bottom environment
            DArray<6> tmp6;
            Contract(1.0,tmp8,shape(7,5,3),env.gb(row-1)[col+1],shape(1,2,3),0.0,tmp6);

            //add left peps to LI8
            DArray<9> tmp9bis;
            Contract(1.0,LI8,shape(6,4),peps(row,col),shape(0,1),0.0,tmp9bis);

            //and another peps
            tmp8.clear();
            Contract(1.0,tmp9bis,shape(4,3,6),peps(row,col),shape(0,1,2),0.0,tmp8);

            //now contract with bottom environment
            DArray<6> tmp6bis;
            Contract(1.0,tmp8,shape(3,6,4),env.gb(row-1)[col],shape(0,1,2),0.0,tmp6bis);

            double val = Dot(tmp6,tmp6bis);

            // (2) construct right hand side

            //add right operator
            DArray<9> tmp9op;
            Contract(1.0,tmp9,shape(3,7,4),rop,shape(1,2,5),0.0,tmp9op);

            //add bottom environment
            DArray<7> tmp7;
            Contract(1.0,tmp9op,shape(8,5,3),env.gb(row-1)[col+1],shape(1,2,3),0.0,tmp7);

            tmp9op.clear();
            Contract(1.0,tmp9bis,shape(4,3,6),lop,shape(0,1,2),0.0,tmp9op);

            //now contract with bottom environment
            DArray<7> tmp7bis;
            Contract(1.0,tmp9op,shape(3,7,4),env.gb(row-1)[col],shape(0,1,2),0.0,tmp7bis);

            DArray<7> perm7;
            Permute(tmp7bis,shape(0,1,2,3,5,4,6),perm7) ;

            val -= 2.0 * Dot(perm7,tmp7);

            return val;

         }//end horizontal
         else if(dir == DIAGONAL_LURD){

            //(1) overlap term

            //RIGHT

            //add right peps to intermediate
            DArray<9> tmp9_right;
            Contract(1.0,RI8,shape(4,6),peps(row,col+1),shape(1,4),0.0,tmp9_right);

            //add second peps
            DArray<8> tmp8;
            Contract(1.0,tmp9_right,shape(3,7,4),peps(row,col+1),shape(1,2,4),0.0,tmp8);

            //add bottom environment
            DArray<6> tmp6_right;
            Contract(1.0,tmp8,shape(7,5,3),env.gb(row-1)[col+1],shape(1,2,3),0.0,tmp6_right);

            //LEFT

            //add left-up peps to LI8
            DArray<9> tmp9_left;
            Contract(1.0,LI8,shape(2,4),peps(row+1,col),shape(0,3),0.0,tmp9_left);

            //and another peps
            tmp8.clear();
            Contract(1.0,tmp9_left,shape(1,7,2),peps(row+1,col),shape(0,2,3),0.0,tmp8);

            //now contract with top environment
            DArray<6> tmp6;
            Contract(1.0,tmp8,shape(0,6,4),env.gt(row)[col],shape(0,1,2),0.0,tmp6);

            DArray<6> tmp6_left;
            Permute(tmp6,shape(5,4,3,1,0,2),tmp6_left);

            double val = Dot(tmp6_left,tmp6_right);

            // (2) operator term

            //add right operator
            DArray<9> tmp9bis;
            Contract(1.0,tmp9_right,shape(3,7,4),rop,shape(1,2,5),0.0,tmp9bis);

            //add bottom environment
            DArray<7> tmp7_right;
            Contract(1.0,tmp9bis,shape(8,5,3),env.gb(row-1)[col+1],shape(1,2,3),0.0,tmp7_right);

            //add left-up peps to b_L
            tmp9_left.clear();
            Contract(1.0,b_L,shape(2,4),peps(row+1,col),shape(0,3),0.0,tmp9_left);

            //add left operator to intermediate
            tmp9bis.clear();
            Contract(1.0,tmp9_left,shape(1,7,2),lop,shape(0,2,4),0.0,tmp9bis);

            //now contract with top environment
            DArray<7> tmp7;
            Contract(1.0,tmp9bis,shape(0,6,4),env.gt(row)[col],shape(0,1,2),0.0,tmp7);

            DArray<7> tmp7_left;
            Permute(tmp7,shape(6,5,3,1,0,4,2),tmp7_left);

            val -= 2.0 * Dot(tmp7_left,tmp7_right);

            return val;

         }
         else{//diagonal ldru

            //(1) overlap term

            //RIGHT

            //add upper right peps to intermediate RI8
            DArray<9> tmp9_right;
            Contract(1.0,peps(row+1,col+1),shape(3,4),RI8,shape(4,2),0.0,tmp9_right);

            //and another
            DArray<8> tmp8;
            Contract(1.0,peps(row+1,col+1),shape(2,3,4),tmp9_right,shape(2,5,4),0.0,tmp8);

            //add top environment
            DArray<6> tmp6_right;
            Contract(1.0,env.gt(row)[col+1],shape(1,2,3),tmp8,shape(1,3,4),0.0,tmp6_right);

            //LEFT

            //add bottom left peps to intermediate LI8
            DArray<9> tmp9_left;
            Contract(1.0,LI8,shape(6,4),peps(row,col),shape(0,1),0.0,tmp9_left);

            //and another
            tmp8.clear();
            Contract(1.0,tmp9_left,shape(4,3,6),peps(row,col),shape(0,1,2),0.0,tmp8);

            //add bottom environment
            DArray<6> tmp6;
            Contract(1.0,tmp8,shape(3,6,4),env.gb(row-1)[col],shape(0,1,2),0.0,tmp6);

            DArray<6> tmp6_left;
            Permute(tmp6,shape(0,1,2,4,3,5),tmp6_left);

            double val = Dot(tmp6_left,tmp6_right);

            // (2) operator term

            //add upper right peps to intermediate b_R
            tmp9_right.clear();
            Contract(1.0,peps(row+1,col+1),shape(3,4),b_R,shape(4,2),0.0,tmp9_right);

            //and right operator to intermediate
            DArray<9> tmp9bis;
            Contract(1.0,rop,shape(2,4,5),tmp9_right,shape(2,5,4),0.0,tmp9bis);

            //add top environment
            DArray<7> tmp7_right;
            Contract(1.0,env.gt(row)[col+1],shape(1,2,3),tmp9bis,shape(1,4,5),0.0,tmp7_right);

            //add left operator to tmp9
            tmp9bis.clear();
            Contract(1.0,tmp9_left,shape(4,3,6),lop,shape(0,1,2),0.0,tmp9bis);

            //add bottom environment
            DArray<7> tmp7;
            Contract(1.0,tmp9bis,shape(3,7,4),env.gb(row-1)[col],shape(0,1,2),0.0,tmp7);

            DArray<7> tmp7_left;
            Permute(tmp7,shape(0,1,4,2,5,3,6),tmp7_left);

            val -= 2.0 * Dot(tmp7_left,tmp7_right);

            return val;

         }

      }

   /**
    * test for the performance of the SVD approximation 
    * for top or bottom rows, i.e. with L and R of order 5 and intermediates of order 7
    * @param dir vertical, horizontal,diagonal lurd or diagonal ldru update
    * @param row , the row index of the bottom site
    * @param col column index of the vertical column
    * @param peps, full PEPS object before update
    * @param L Left environment contraction
    * @param R Right environment contraction
    * @param LI7 left intermediate object
    * @param RI7 left intermediate object
    * @param b_L left intermediate object
    * @param b_R right intermediate object
    */
   template<>
      void svd_test(const PROP_DIR &dir,int row,int col,PEPS<double> &peps,const DArray<6> &lop,const DArray<6> &rop,

            const DArray<5> &L,const DArray<5> &R,const DArray<7> &LI7,const DArray<7> &RI7,int Dr){

         if(dir == VERTICAL){

            if(row == 0){

               if(col == 0){

                  DArray<8> tmp8;
                  Contract(1.0,lop,shape(1,3),rop,shape(4,3),0.0,tmp8);

                  //now get the "cost function limit"
                  DArray<6> tmp6;
                  Contract(1.0,tmp8,shape(0,1,2,4,6),tmp8,shape(0,1,2,4,6),0.0,tmp6);

                  DArray<6> tmp6bis;
                  Permute(tmp6,shape(1,4,2,5,0,3),tmp6bis);

                  cout << -blas::dot(tmp6bis.size(),tmp6bis.data(),1,RI7.data(),1) << endl;

                  //svd the fucker
                  DArray<5> UL;//left unitary
                  DArray<5> VR;//right unitary

                  DArray<1> S;
                  Gesvd ('S','S', tmp8, S,UL,VR,0);
                  cout <<  S << endl;

                  //take the square root of the sv's
                  for(int i = 0;i < S.size();++i)
                     S(i) = sqrt(S(i));

                  //and multiply it left and right to the tensors
                  Dimm(S,VR);
                  Dimm(UL,S);

                  DArray<8> tmp8_reduced;
                  Contract(1.0,UL,shape(4),VR,shape(0),0.0,tmp8_reduced);

                  //now get the "reduced value for the cost function"
                  Contract(1.0,tmp8_reduced,shape(0,1,2,4,6),tmp8_reduced,shape(0,1,2,4,6),0.0,tmp6);

                  Permute(tmp6,shape(1,4,2,5,0,3),tmp6bis);

                  double overlap = blas::dot(tmp6bis.size(),tmp6bis.data(),1,RI7.data(),1);

                  Contract(1.0,tmp8_reduced,shape(0,1,2,4,6),tmp8,shape(0,1,2,4,6),0.0,tmp6);

                  Permute(tmp6,shape(1,4,2,5,0,3),tmp6bis);

                  cout << overlap - 2.0 * blas::dot(tmp6bis.size(),tmp6bis.data(),1,RI7.data(),1) << endl;

               }
            }
         }
      }

   /**
    * test for the performance of the SVD approximation 
    * for top or bottom rows, i.e. with L and R of order 5 and intermediates of order 7
    * @param dir vertical, horizontal,diagonal lurd or diagonal ldru update
    * @param row , the row index of the bottom site
    * @param col column index of the vertical column
    * @param peps, full PEPS object before update
    * @param L Left environment contraction
    * @param R Right environment contraction
    * @param LI7 left intermediate object
    * @param RI7 left intermediate object
    * @param b_L left intermediate object
    * @param b_R right intermediate object
    */
   template<>
      void svd_test(const PROP_DIR &dir,int row,int col,PEPS<double> &peps,const DArray<6> &lop,const DArray<6> &rop,

            const DArray<6> &L,const DArray<6> &R,const DArray<8> &LI8,const DArray<8> &RI8,int Dr){

         if(dir == VERTICAL){

            if(row == 0){

               if(col == 0){

                  cout << "fuckityfuckityfuck" << endl;

               }
            }
         }
      }


} 
