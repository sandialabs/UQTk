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
/// \file kldecompuni.cpp
/// \author B. Debusschere, K. Sargsyan, C. Safta, 2007 -
/// \brief Karhunen-Loeve decomposition class

#include <math.h>
#include "kle.h"
#include "deplapack.h"
#include "error_handlers.h"

KLDecompUni::KLDecompUni(const Array1D<double>& tSamples)
{
  // Initializations
  this->Init() ;

  const size_t n_t = tSamples.XSize();	// Number of time samples, which defines overall dimensions

  // Properly size arrays for matrix and weights, which are independent of the
  // number of KL modes requested.
  whcwh_.Resize(n_t,n_t,0.e0);
  w_.    Resize(n_t,    0.e0);
  wh_.   Resize(n_t,    0.e0);

  // Compute the weights for Nystrom's method for Fredholm integral equation solution
  // using the trapezoidal rule as the quadrature rule
  w_(0) = (tSamples(1)-tSamples(0))/2.e0;
  for(size_t i=1; i < n_t-1; i++){
    w_(i) = (tSamples(i+1)-tSamples(i-1))/2.e0;
  }
  w_(n_t-1) = (tSamples(n_t-1)-tSamples(n_t-2))/2.e0;

  // Get the square roots of these weights
  for( size_t i=0; i < n_t; i++ ){
    wh_(i) = sqrt(w_(i));
  }

}

KLDecompUni::KLDecompUni()
{

  // Initializations
  this -> Init() ;

}

void KLDecompUni::Init()
{

  // Initializations
  decomposed_ = false;

  // settings for the LAPACK eigensolver
  jobz_     = 'V' ;  // get eigenvalues and eigenvectors
  eigRange_ = 'I' ;  // get eigenvalues within a given range of indices (from smaller to larger)
  uplo_     = 'U' ;  // upper triangular storage
  vl_       = 0.e0;  // not used
  vu_       = 1.e0;  // not used
  absTol_   = 0.e0;  // let routine determine appropriate tolerance

}

void KLDecompUni::SetWeights(const Array1D<double> &weights)
{

  const int npts = weights.XSize();	// Number of time samples, which defines overall dimensions
  this->SetWeights(weights.GetConstArrayPointer(), npts);

  return ;

}

void KLDecompUni::SetWeights(const double *weights, const int npts)
{

  const size_t n_t = npts;	// Number of time samples, which defines overall dimensions

  w_.Resize(npts);
  for ( size_t i=0;i<npts;i++) w_[i] = weights[i];

  // Properly size arrays for matrix and weights, which are independent of the
  // number of KL modes requested.
  whcwh_.Resize(n_t,n_t,0.e0);
  wh_.   Resize(n_t,    0.e0);

  // Get the square roots of these weights
  for( size_t i=0; i < n_t; i++ ){
    wh_(i) = sqrt(w_(i));
  }

  return ;

}

int KLDecompUni::decompose(const Array2D<double>& corr, const int &nKL) {

  return (this->decompose(corr.GetConstArrayPointer(), nKL));

}

int KLDecompUni::decompose(const double *corr, const int &nKL) {

  const size_t n_t = whcwh_.XSize();  // dimension (rank) of the matrix, which determines
                                      // max # of eigenvalues

  // Set the range of indices of the eigenvalues we look for (in ascending order)
  // to get the largest nKL eigenvalues
  const size_t one = 1;
  il_ = max(one, n_t - nKL + one);
  iu_ = n_t;

  // Populate the upper triangular part of the matrix, with each entry being
  // \sqrt{w(i)} C(i,j) \sqrt{w(j)}
  // This needs to be done every time we call the eigenvalue solver as this solver
  // destroys the matrix entries.
  for(size_t i=0; i < n_t; i++){
    for(size_t j = i; j < n_t; j++){
      whcwh_(i,j) = wh_(i)*corr[j*n_t+i]*wh_(j);
    }
  }

  // Set the dimensions of the arrays that will hold the eigenvalues and
  // eigenvectors and reset all elements to 0.e0
  eig_values_.Resize(n_t,    0.e0);
  KL_modes_.  Resize(n_t,nKL,0.e0);

  // Work arrays
  Array1D<double> work_eigval(n_t,    0.e0);  // Temporary storage of eigenvalues
  Array2D<double> work_eigvec(n_t,nKL,0.e0);  // Temporary storage of eigenvectors

  int  lwork = 10 * n_t;             // size of work array (see LAPACK routine dsyevx.f for more info)
  Array1D<double> work(lwork,0.e0);  // double precision work array
  Array1D<int>    iwork(5*n_t,0);    // integer work array

  ifail_.Resize(n_t,0);  // on output: contains indices of eigenvectors that failed to converge
  eig_info_ = 0;         // info on success of the eigenvector solutions

  int n_eig = 0;         // on output: has number of eigenvalues that were obtained
  int n_t_int=n_t;

  // Call the eigenproblem solver
  FTN_NAME(dsyevx)(&jobz_, &eigRange_, &uplo_, &n_t_int, whcwh_.GetArrayPointer(),
                   &n_t_int, &vl_, &vu_, &il_, &iu_,
                   &absTol_, &n_eig, work_eigval.GetArrayPointer(),
                   work_eigvec.GetArrayPointer(), &n_t_int,
                   work.GetArrayPointer(), &lwork, iwork.GetArrayPointer(),
                   ifail_.GetArrayPointer(), &eig_info_);

  // Set the decomposed flag if successful eigenvalue solve
  if ( eig_info_ == 0 ){
    decomposed_ = true;
  } else {
    std::string err_message = "Something in the eigensolve went wrong. Check error code: ";
    err_message += eig_info_;
    err_message += " in the LAPACK routine dsyevx.f";
    throw Tantrum(err_message);
  }

  // Reverse the order in the array with eigenvalues so that they are ordered
  // in descending order (LAPACK dsyevx returns them in ascending order)
  for(size_t i=0; i <(size_t) n_eig; i++){
    eig_values_(i) = work_eigval(n_eig-1-i);
  }

  // Rearrange the eigenvectors accordingly and scale them with the square root of the
  // integration weights to get the eigenmodes of the autocorrelation matrix.
  for(size_t i=0; i < n_t; i++){
    for(size_t j=0; j < (size_t) n_eig; j++){
      KL_modes_(i,j) = work_eigvec(i,n_eig-1-j)/wh_(i);
    }
  }

  // Return the number of obtained eigenvalues
  return n_eig;

}

void KLDecompUni::KLproject(const Array2D<double> &realiz, Array2D<double> &xi)
{

  // Get dimensions
  int nKL = xi.XSize();
  const size_t n_t = realiz.XSize();
  const size_t n_r = realiz.YSize();

  // dimension check
  if ( n_r != xi.YSize() ) {
    printf("KLproject(): dimension error"); exit(1);
  }

  Array1D<double> mean_realiz(n_t,0.e0);
  this->meanRealiz(realiz,mean_realiz);

  for(int ikl=0; ikl < nKL; ikl++)
    for(int k=0; k < (int) n_r; k++) {
      xi(ikl,k) = 0.0 ;
      for(size_t i=0; i < n_t; i++)
        xi(ikl,k) += KL_modes_(i,ikl) * ( realiz(i,k) - mean_realiz(i) )
                    * w_(i) / sqrt(eig_values_(ikl)) ;
    }

  return;

}

const Array1D<double>& KLDecompUni::eigenvalues() const
{

  if(!decomposed_){
    std::string err_message = "Eigenvalues are not yet available";
    throw Tantrum(err_message);
  } else {
    return eig_values_;
  }

}

void KLDecompUni::eigenvalues(const int nEIG, double *eigs) const
{

  if(!decomposed_){
    std::string err_message = "Eigenvalues are not yet available";
    throw Tantrum(err_message);
  } else {
    if (nEIG > eig_values_.XSize() ) {
      std::cout<<nEIG<<" "<<eig_values_.XSize() <<std::endl;
      std::string err_message = "eigenvalues(): number of eigenvalues more than available";
      throw Tantrum(err_message);
    }
    for(size_t i=0; i < nEIG; i++)
      eigs[i] = eig_values_(i);
  }

  return ;

}

const Array2D<double>& KLDecompUni::KLmodes() const
{

  if(!decomposed_){
    std::string err_message = "KL modes are not yet available";
    throw Tantrum(err_message);
  } else {
    return KL_modes_;
  }

}

void KLDecompUni::KLmodes(const int npts, const int nKL, double *klModes ) const
{

  if(!decomposed_){
    std::string err_message = "KL modes are not yet available";
    throw Tantrum(err_message);
  } else {
    if (npts != KL_modes_.XSize() ) {
      std::string err_message = "KLmodes(): number of grid points does not match";
      throw Tantrum(err_message);
    }
    if (nKL != KL_modes_.YSize() ) {
      std::string err_message = "KLmodes(): number of modes does not match";
      throw Tantrum(err_message);
    }
    for(size_t j=0; j < nKL; j++)
      for(size_t i=0; i < npts; i++)
        klModes[j*npts+i] = KL_modes_(i,j);
  }

  return ;

}

void KLDecompUni::meanRealiz(const Array2D<double> &realiz, Array1D<double> &mean_realiz)
{
  const size_t n_t = realiz.XSize();
  const int    n_r = realiz.YSize();

  // dimension check
  if (n_t != mean_realiz.XSize()) {
    printf("meanRealiz: dimension error"); exit(1);
  }

  for(size_t i=0; i < n_t; i++){
    mean_realiz(i) = 0.0 ;
    for(int k=0; k < n_r; k++){
      mean_realiz(i) += realiz(i,k)/n_r;
    }
  }

  return;

}

void KLDecompUni::truncRealiz(const Array1D<double>& meanrea,const Array2D<double>& xi,
                              const int& nKL, Array2D<double>& trunc_realiz)
{

  const size_t n_t = meanrea.XSize();
  const int nsample=xi.XSize();

  // dimension check
  if ( ( nsample != (int) trunc_realiz.XSize() ) ||
       ( n_t     !=       trunc_realiz.YSize() ) ) {
    printf("truncRealiz: dimension error"); exit(1);
  }
  if ( nKL > (int) xi.YSize() ) {
    printf("truncRealiz: there are not enough xi variables."); exit(1);
  }

  for(int isample=0; isample < nsample; isample++){
    for(size_t i=0; i < n_t; i++){
      trunc_realiz(isample,i) = meanrea(i);
      for(int ikl=0; ikl < nKL; ikl++){
        trunc_realiz(isample,i) += xi(isample,ikl)*sqrt(eig_values_(ikl))*KL_modes_(i,ikl);
      }
    }
  }

  return;

}
