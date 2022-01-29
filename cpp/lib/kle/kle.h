/* =====================================================================================

                      The UQ Toolkit (UQTk) version 3.1.2
                          Copyright (2022) NTESS
                        https://www.sandia.gov/UQToolkit/
                        https://github.com/sandialabs/UQTk

     Copyright 2022 National Technology & Engineering Solutions of Sandia, LLC (NTESS).
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

     Questions? Contact the UQTk Developers at <uqtk-developers@software.sandia.gov>
     Sandia National Laboratories, Livermore, CA, USA
===================================================================================== */
/// \file kldecompuni.h
/// \author B. Debusschere, K. Sargsyan, C. Safta, 2007 -
/// \brief Header for Karhunen-Loeve decomposition class


#ifndef KLDECOMPUNI_H_SEEN
#define KLDECOMPUNI_H_SEEN

#include <iostream>
#include "ftndefs.h"
#include "Array1D.h"
#include "Array2D.h"
// #include "Array3D.h"

/*!
  \class  KLDecompUni
  \brief  Computes the Karhunen-Loeve decomposition of a univariate stochastic process
  \details
  @f[
  F(t,\theta) = \left < F(t,\theta) \right >_{\theta}
                + \sum_{k=1}^{\infty} \sqrt{\lambda_k} f_k(t) \xi_k@f]
*/
class KLDecompUni {

  public:

  /*!
  \brief
  Constructor that takes the autocorrelation matrix "corr" (\f$C\f$) of the process
  we are studying as well as the array "tsamples" (\f$t\f$) with the points in time where
  snapshots of the system were taken.

  \details
  Constructs weights (\f$w\f$) needed for the Nystrom method to solve the Fredholm integral equation
  @f[ \int C(s,t)f(t)dt=\lambda f(s) \rightarrow \sum w_j C(s_i,t_j) f_k(t_j) = \lambda_k f_k(s_i)@f]
  */
  KLDecompUni(const Array1D<double>& tSamples);

  KLDecompUni();

  /// \brief Destructor
  ~KLDecompUni() {};

  void Init() ;

  ///  \brief Set weights for computing the integral needed for
  ///  Nystrom's method for solving the Fredholm integral equation
  void SetWeights(const Array1D<double>& weights) ;

  ///  \brief Set weights for computing the integral needed for
  ///  Nystrom's method for solving the Fredholm integral equation
  void SetWeights(const double *weights, const int npts) ;

  /*!
  \brief Perform KL decomposition into nKL modes and return actual number of
   modes that were obtained
  \details Further manipulation of the discretized Fredholm equation leads to the eigenvalue problem
    @f[A g=\lambda g @f]
  where \f$A=W K W\f$ and \f$g=Wf\f$, with \f$W\f$ being the diagonal matrix, \f$W_{ii}=\sqrt{w_i}\f$ and
  \f$K_{ij}=Cov(t_i,t_j)\f$. Solutions consist of pairs of eigenvalues \f$\lambda_k\f$ and KL modes \f$f_k=W^{-1}g_k\f$.
   */
  int decompose(const Array2D<double>& corr, const int& nKL);

  /*!
  \brief Perform KL decomposition into nKL modes and return actual number of
   modes that were obtained
  \details Further manipulation of the discretized Fredholm equation leads to the eigenvalue problem
    @f[A g=\lambda g @f]
  where \f$A=W K W\f$ and \f$g=Wf\f$, with \f$W\f$ being the diagonal matrix, \f$W_{ii}=\sqrt{w_i}\f$ and
  \f$K_{ij}=Cov(t_i,t_j)\f$. Solutions consist of pairs of eigenvalues \f$\lambda_k\f$ and KL modes \f$f_k=W^{-1}g_k\f$.
   */
  int decompose(const double *corr, const int &nKL);

  /*!
  \brief Project realizations \f$F(t,\theta_l)\f$ to the KL modes and store them in xi (\f$\xi_k\f$)

  \details Samples of random variables \f$\xi_k\f$ are obtained by projecting
  realizations of the random process \f$F\f$ on the eigenmodes \f$f_k\f$
  @f[ \left.\xi_k\right\vert_{\theta_l}=\left <F(t,\theta_l)-\left <
    F(t,\theta) \right >_{\theta}, f_k(t) \right >_t/\sqrt{\lambda_k} @f]
  ... or numerically
  @f[
  \left.\xi_k\right\vert_{\theta_l}=\sum_{i=1}^{N_p} w_i\left(F(t_i,\theta_l)-\left <
    F(t_i,\theta) \right >_{\theta} \right) f_k(t_i)/\sqrt{\lambda_k} @f]
  */
  void KLproject(const Array2D<double>& realiz, Array2D<double>& xi);

  /// \brief Get eigenvalues in descending order
  const Array1D<double>& eigenvalues() const;
  void eigenvalues(const int nEIG, double *eigs) const;

  /// \brief Get associated KL modes
  const Array2D<double>& KLmodes() const;

  /// \brief Get associated KL modes
  void KLmodes(const int npts, const int nKL, double *klModes ) const;

  /// \brief Calculate (in meanRealiz) the mean realizations
  void meanRealiz(const Array2D<double>& realiz, Array1D<double>& mean_realiz);

  /*!
  \brief
  Returns the truncated KL sum
  \details
  @f[
    F(t_i,\theta_l) = \left < F(t_i,\theta) \right >_{\theta}
                          + \sum_{k=1}^{nKL} \sqrt{\lambda_k} f_k(t_i) \left. \xi_k\right\vert_{\theta_l}

  @f]
  */
  void truncRealiz(const Array1D<double>& meanrea,const Array2D<double>& xi,const int& nKL, Array2D<double>& trunc_realiz);

  private:
  /// \brief Dummy default constructor, which should not be used as it is not well defined
  //KLDecompUni() {};

  /// \brief Dummy copy constructor, which should not be used as it is currently not well defined
  KLDecompUni(const KLDecompUni &) {};

  /// \brief Flag to determine whether KL decomposition has taken place (and consequently
  /// that the interal data structures contain meaningful eigenvalues and vectors ... )
 bool decomposed_;

 /// \brief Matrix to hold the upper triangular part of the matrix to get eigenvalues of
  Array2D<double> whcwh_;

  /// \brief Array to hold weights for Nystrom's method for Fredholm integral equation solution
  Array1D<double> w_;

  /// \brief Array to hold square roots of weights
  Array1D<double> wh_;

  /// \brief Option to determine what to compute (eigenvalues and eigenvectors)
  char jobz_ ;
  /// \brief Option to set the type of range for eigenvalues
  char eigRange_;
  /// \brief Option to indicate how matrix is stored
  char uplo_;
  /// \brief Lower bound for range of eigenvalues
  double vl_;
  /// \brief Upper bound for range of eigenvalues
  double vu_;
  /// \brief Lower index of range of eigenvalues requested
  int il_;
  /// \brief Upper index of range of eigenvalues requested
  int iu_;
  /// \brief Absolute tolerance for convergence
  double absTol_;

  /// \brief Array to store eigenvalues
  Array1D<double> eig_values_;
  /// \brief Matrix to store KL modes
  Array2D<double> KL_modes_;

  /// \brief info on success of the eigenvector solutions
  int eig_info_;
  /// \brief Array to store indices of eigenvectors that failed to converge
  Array1D<int> ifail_;

};

#endif /* KLDECOMPUNI_H_SEEN */
