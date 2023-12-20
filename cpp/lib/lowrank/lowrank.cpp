/* =====================================================================================

                      The UQ Toolkit (UQTk) version @UQTKVERSION@
                          Copyright (@UQTKYEAR@) NTESS
                        https://www.sandia.gov/UQToolkit/
                        https://github.com/sandialabs/UQTk

     Copyright @UQTKYEAR@ National Technology & Engineering Solutions of Sandia, LLC (NTESS).
     Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government
     retains certain rights in this software.

     This file is part of The UQ Toolkit (UQTk)

     UQTk is open source software: you can redistribute it and/or modify
     it under the terms of BSD 3-Clause License

     UQTk is distributed in the hope that it will be useful,
     but WITHOUT ANY WARRANTY; without even the implied warranty of
     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
     BSD 3 Clause License for more details.

     You should have received a copy of the BSD 3 Clause License
     along with UQTk. If not, see https://choosealicense.com/licenses/bsd-3-clause/.

     Questions? Contact the UQTk Developers at https://github.com/sandialabs/UQTk/discussions
     Sandia National Laboratories, Livermore, CA, USA
===================================================================================== */
/// \file lowrank.cpp
/// \author Prashan Rai  2016 -
/// \brief Low Rank representation class


#include <math.h>
#include <cfloat>
#include <iostream>


#include "lowrank.h"
#include "PCBasis.h"
#include "error_handlers.h"
#include "gen_defs.h"
//#include "alltools.h"
#include "arrayio.h"


#include "tools.h"
#include "lbfgs_routines.h"

#include <assert.h>



// Constructor
CanonicalTensor::CanonicalTensor(){
  core_ = Array1D<double>(1,0.0);
  size_ = Array1D<int>(1,0);
  order_ = 1.0;
  space_.PushBack(Array2D<double>(1,1,0));
}

// Constructor that initializes core and space
CanonicalTensor::CanonicalTensor(Array1D<double>& core, Array1D< Array2D<double> > space)
{
  core_ = core;
  space_ = space;
  order_ = space.XSize();
  Array1D<int> size_(order_,0);
  for(int i = 0; i<order_; i++){
    size_(i) = space_[i].XSize();
  }
}

// Constructor for customized initialization
void CanonicalTensor::create(string generator, int rank, Array1D<int>& size)
{
  size_ = size;
  order_ = size_.Length();
  Array1D<double> core(rank,1.0);
  space_.Clear();
  if(!generator.compare("ones")){
    for (int i=0;i<order_;i++)
      space_.PushBack(Array2D<double>(size(i),rank,1));
  }
  else if(!generator.compare("zeros")){
    for (int i=0;i<order_;i++)
      space_.PushBack(Array2D<double>(size(i),rank,0));
  }
  else{ //Need to implement random and greedy initilaization later
    cout<<"Not yet programmed"<<"\n";
  }
  core_ = core;
}

// Get the coefficient vector of univariate function in dimension dim for a given rank
Array1D<double> CanonicalTensor::getSpaceVec(int rank, int dim){
  Array1D<double> spaceVec;
  getCol(space_[dim],rank,spaceVec);
  return spaceVec;
}

// Set the coefficient vector of univariate function in dimension dim for a given rank
void CanonicalTensor::setSpaceVec(int rank, int dim, Array1D<double>& column){
  space_[dim].eraseCol(rank);
  space_[dim].insertCol(column, rank);
}

// Display the Tensor
void CanonicalTensor::display() {
  int order = core_.Length();
  cout << "====================================================" << endl;
  cout << "Canonical Tensor " << endl;
  cout << "----------------------------------------------------" << endl;
  cout<<"Order: "<<order_<< "\n";
  cout<<"Rank: "<<core_.Length()<<"\n";
  cout<<"Size:";
  cout << "[( ";
  for (int mu = 0; mu < size_.Length(); mu++){
      cout << size_(mu) << ", ";
  }
  cout << ")]" << endl;
  cout << "----------------------------------------------------" << endl;
}

// Scalar product of two canonical tensors (useful for finding norm of a canonical tensor)
double CanonicalTensor::prodscal(CanonicalTensor& u, CanonicalTensor& v){
  //todo check if the tensors are of the same size
  double dotprod = 0.0;
  for (int r1=0;r1<u.Rank();r1++){
    for (int r2=0;r2<v.Rank();r2++){
      double prod = 1.0;
      for (int k=0;k<u.order_;k++){
        Array1D<double> w1_k, w2_k;
        w1_k = u.getSpaceVec(r1,k);
        w2_k = v.getSpaceVec(r2,k);
        prod = prod * dot(w1_k,w2_k);
      }
      dotprod = dotprod + u.core_(r1)*v.core_(r2)*prod;
    }
  }
  return dotprod;
}

// Norm of a canonical tensor
double CanonicalTensor::norm()
{
   return sqrt(CanonicalTensor::prodscal(*this,*this));
}

// To be implemented

/*void CanonicalTensor::printTensor(){
  char* tfile;
  tfile = "space.dat";
  write_datafile(space_[1], tfile);
}

*/

// Polynomial Chaos basis parameters
PCBases::PCBases(Array1D<string> pctype, Array1D<int> order)
{
  pctype_ = pctype;
  order_ = order;
}

// Evaluate 1D PC basis at xval(:,dim) and return in 2D array 'bb'
void PCBases::functionEval(Array2D<double>& xval, int dim, Array2D<double>& bb){
  PCBasis currPCBasis;
  Array1D<double> xcol(xval.XSize(),getOrder(dim));
  getCol(xval,dim,xcol);
  currPCBasis.Eval1dBasisAtCustPoints(bb,order_(dim),xcol);
}

// Evaluate 1D Polynomial basis at xval(:,dim) and return in 2D array 'bb'
void PLBases::functionEval(Array2D<double>& xval, int dim, Array2D<double>& bb){
  int nx = xval.XSize();
  int nbas = getOrder(dim)+1;
  bb.Resize(nx,nbas,0.e0);
  Array1D<double> xcol(nx,getOrder(dim));
  getCol(xval,dim,xcol);
  for(int is=0;is<nx;is++)
    for(int ipc=0;ipc<=getOrder(dim);ipc++){
      bb(is,ipc)=pow(xcol(is),ipc);
    }
}

// Default functional tensor constructor (Tensor with basis information)
FunctionalTensor::FunctionalTensor(){
  tensor_ = CanonicalTensor();
  bases_ = Array1D<FunctionalBases* >(tensor_.getOrder(),NULL);
  fdims_ = 0;
}

// Functional tensor constructor (that initializes coefficient tensor only)
FunctionalTensor::FunctionalTensor(CanonicalTensor& f){
  tensor_ = f;
  bases_ = Array1D<FunctionalBases* >(tensor_.getOrder(),NULL);
  fdims_ = f.getOrder();
}

// Functional tensor constructor (initializes both coefficient tensor and basis)
FunctionalTensor::FunctionalTensor(CanonicalTensor& f, Array1D<FunctionalBases* >& bases){
  tensor_ = f;
  bases_ = bases;
  fdims_ = f.getOrder();
}

// Functional tensor constructor (initializes both coefficient tensor and basis)
Array1D<double> FunctionalTensor::tensorEval(Array2D<double>& xval){
  int N = xval.XSize();
  Array1D< Array2D<double> > H;
  if (bases_(0)!=NULL){
    Array2D<double> tempH;
    for(int i=0;i<bases_.Length();i++){
      for(int j=0;j<(bases_(i))->getNumOrder();j++){
        (bases_(i))->functionEval(xval,j,tempH);
        H.PushBack(tempH);
      }
     }
   }
   else{
    cout<<"Bases not set for the functional tensor\n";
    exit(1);
  }
  Array1D<double> yval;
  yval = tensorEval(H,N);
  return yval;
}

// Evaluated functional tensor (H is basis evaluated at input samples)
Array1D<double> FunctionalTensor::tensorEval(Array1D< Array2D<double> >& H, int N){
  Array1D<double> yval(N), alpha = tensor_.getCore(); // get singular values
  int order = H.XSize(); // isotropic order of basis
  //assert();
  Array1D< Array2D<double> > fH(order);
  Array2D<double> tempSpace, A;
  for (int mu=0; mu<fdims_; mu++){  // Dimension-wise multiplication of basis evals and function coefficients
    tempSpace = tensor_.getSpace(mu);
    prodAlphaMatMat(H[mu],tempSpace,1.0,fH[mu]);
  }
  int rank = alpha.Length();
  A.Resize(N,rank,1.0);
  for(int mu=0; mu<fdims_; mu++){
    A = dotmult(A,fH[mu]);
  }
  prodAlphaMatVec(A, alpha, 1.0, yval);
  return yval;
}
// ALS constructor
CanonicalTensorALSLeastSquares::CanonicalTensorALSLeastSquares(){
  Array1D<FunctionalBases* > bases_(1,NULL);
  alg_ = "direct"; // two types of algo: greedy and direct (only direct is implementes)
  tol_ = 1.0e-08;
  initializationType_ = "ones";
  rank_ = 1;
  maxIterations_ = 50;
  lambda_ = 0.0;
}

// Main function for direct ALS approximation
CanonicalTensor CanonicalTensorALSLeastSquares::solveDirect(Array2D<double>& xval,Array1D<double>& yval,bool vflag){
  Array1D< Array2D<double> > H; // H{0} ... H{dim} are matrices (Nsamples x basis order) of basis evaluated at sample points
  int N = yval.Length(); // Number of samples
  Array2D<double> tempH;
  if(bases_(0)!=NULL){  //i.e. if the user has provided the bases_ of type FunctionalBases
     for(int i=0;i<bases_.Length();i++){
      for(int j=0;j<(bases_(i))->getNumOrder();j++){
        (bases_(i))->functionEval(xval,j,tempH);
        H.PushBack(tempH);
      }
     }
     basesEval_ = H;
   }
   else{
    H = basesEval_;
  }
  int order = H.XSize(); // order of tensor i.e. dimension of surrogate

  Array1D<int> size(order,0.0);
  for (int i=0; i<order; i++)
    size(i) = H[i].YSize();
//  initializationType_ = "ones";
  //assert(!initializationType_.compare("none"));
  CanonicalTensor f, f0;
  if (!initializationType_.compare("ones"))
    f.create("ones", rank_, size);
  else if(!initializationType_.compare("zeros"))
    f.create("zeros", rank_, size);
  else{
    cout<<"Initialization not yet coded"<<"\n";
    exit(1);
  }

  //f.display();

  if (maxIterations_ == 0) // Number of ALS iterations
    return f;


  Array1D< Array2D<double> > fH(order); // fH is coefficient multiplied with basis evals
  Array2D<double> tempSpace;
  for (int mu=0; mu<order; mu++){
    tempSpace = f.getSpace(mu);
    prodAlphaMatMat(H[mu],tempSpace,1.0,fH[mu]);
  }
  Array2D<double> A(N,0);
  for(int k=0; k<maxIterations_; k++){ // ALS iteration begins
    Array1D<double> alpha, BTb;
    Array2D<double> BTB, invBTB;
    f0 = f;
    for (int mu=0; mu<order; mu++){ // Dimension-wise optimization loop
      Array2D<double> fHmu(N,rank_,1.0);
      for (int nomu=0; nomu<order; nomu++){ // loop for multiplying all dimension univariate function except mu
        if (mu!=nomu){
          fHmu = dotmult(fHmu,fH[nomu]);
        }
      }
      Array2D<double> Ar(N,size(mu)); //Ar is measurement matrix for 1<=r<=rank
      A.Resize(N,0);
      for(int r=0; r<rank_; r++){ //rank loop
        Array2D<double> Ap(N,0);
        Array1D<double> fixPhi;
        getCol(fHmu,r,fixPhi);
        for(int p=0; p<size(mu); p++){ //Sp loop
          Ap.insertCol(fixPhi,p); //
        }
        Ar = dotmult(Ap,H[mu]);
        A.insertCol(Ar,r*size(mu));
      }
        Array1D<double> a, ATb, yres = yval; //declarations to solve Least-squares problem
        Array2D<double> ATA, invATA, ATALambdaI, LambdaI(A.YSize(),A.YSize(),0.0);
        for(int l=0;l<A.YSize(); l++){
          LambdaI(l,l) = lambda_;
        }
        prodAlphaMatTMat(A,A,1.0,ATA); // sequence of code for solving A'A=A'b
        ATALambdaI = add(ATA,LambdaI);
        invATA = INV(ATALambdaI);
        prodAlphaMatTVec(A,yres,1.0,ATb);
        prodAlphaMatVec(invATA, ATb, 1.0, a); // solves the least squares problem
        //printarray(a);
        Array1D<double> spZeros(a.Length(),0.0), vecZeros(size(mu),0.0);
        if (is_equal(a, spZeros)){ // check if the coefficient vector is not zero
          cout<<"All zero coefficient returned, returning previous iterate"<<"\n";
          exit (1);
        }
        Array2D<double> updSpace(size(mu),rank_);
        fold_1dto2d_colfirst(a,updSpace);
        Array1D<double> newSpaceVec(size(mu)), norma(rank_);
        for (int r=0; r<rank_; r++){
          getCol(updSpace,r,newSpaceVec);
          if (is_equal(newSpaceVec,vecZeros)){ // check for all zero coefficients returned from least-squares solve
            cout<<"degenerate case: one factor is zero";
            exit(1);
          }
          else{
            norma(r) = norm(newSpaceVec); // calculating singular values (i.e. norm of all univariate function coefficients for a given rank)
            newSpaceVec = scale(newSpaceVec,(1.0/norma(r)));
            f.setSpaceVec(r, mu, newSpaceVec);
            updSpace.eraseCol(r);
            updSpace.insertCol(newSpaceVec,r); // update space
          }
        }
        f.setCore(norma); // set core
        prodAlphaMatMat(H[mu],updSpace,1.0,fH[mu]);
        // Array1D<double> test;
        // test = f.getCore();
        // printarray(test);
    }

    Array2D<double> B(N,rank_,1.0); // Applying least-squares approximation for updating singular values
                                    // considering each rank component as a basis function
    for(int nu=0; nu<order; nu++){
      B = dotmult(B,fH[nu]);
    }
    prodAlphaMatTMat(B, B, 1.0, BTB) ;
    invBTB = INV(BTB);
    prodAlphaMatTVec(B, yval, 1.0, BTb);
    prodAlphaMatVec(invBTB, BTb, 1.0, alpha);
    f.setCore(alpha);

    double nf = f.norm();
    double nf0 = f0.norm();
    double stagnation = 2*fabs(nf-nf0)/(nf+nf0);
    if(vflag){
      if((k+1)%10==0){
        printf("Iteration %3d - Stagnation = %e\n",k+1,stagnation); // find error between iterates
      }
    }
    f0 = f;
  }

  FunctionalTensor fun(f);

  return f;
}


