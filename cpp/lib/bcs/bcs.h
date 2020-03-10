/* =====================================================================================

                      The UQ Toolkit (UQTk) version 3.1.0
                          Copyright (2020) NTESS
                        https://www.sandia.gov/UQToolkit/
                        https://github.com/sandialabs/UQTk

     Copyright 2020 National Technology & Engineering Solutions of Sandia, LLC (NTESS).
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
/// \file bcs.h
/// \brief Header for the implementations of Bayesian compressive sensing algorithm.

#ifndef BCS_H
#define BCS_H

#include "Array1D.h"
#include "Array2D.h"


#define MAX_IT 1000



/// \brief Implements weighted version of the original Bayesian Compressive Sensing algorithm
/// \note This function has been written relying on the algorithm and MATLAB code presented in
/// http://ivpl.eecs.northwestern.edu/research/projects/bayesian-compressive-sensing-using-laplace-priors
/// and references therein
/// \todo The array manipulations are not optimized - perhaps they need to be reconsidered using,
/// say, fortran matrix-vector manipulation routines
/// \param[in] PHI         : design matrix
/// \param[in] y           : data vector
/// \param[in,out] sigma2  : initial noise variance (usually var(y)/1e2)
///                        : re-estimated on output
/// \param[in] eta         : stopping criterion (usually 1e-5)
/// \param[in] lambda_init : regularization weight vector, if empty array, it automatically computes the optimal, uniform weights
/// \param[in] adaptive    : generate basis for adaptive CS (usually 0)
/// \param[in] optimal     : use the rigorous implementation of adaptive CS (usually 1)
/// \param[in] scale       : diagonal loading parameter (usually 0.1)
/// \param[in] verbose     : verbosity flag
/// \param[out] weights    : sparse weights
/// \param[out] used       : the positions of sparse weights
/// \param[out] errbars    : one standard deviation around the sparse weights
/// \param[out] basis:     : if adaptive==1, then this is the next projection vector, see \cite Ji:2008
/// \param[out] alpha:     : estimated sparse hyperparameters (1/gamma), see \cite Babacan:2010
/// \param[out] Sig        : covariance matrix of the weights
void WBCS(Array2D<double> &PHI, Array1D<double> &y, double &sigma2,
                 double eta, Array1D<double> &lambda_init,
                 int adaptive, int optimal, double scale, int verbose,
                 Array1D<double> &weights, Array1D<int> &used,
                 Array1D<double> &errbars, Array1D<double> &basis,
                 Array1D<double> &alpha, Array2D<double> &Sig);



/// \brief Essentially same functionality as WBCS, but slightly altered I/O.
/// \note Kept for backward compatibility with PyUQTk and BCS tests
void BCS(Array2D<double> &PHI, Array1D<double> &y, double &sigma2,
                 double eta, Array1D<double> &lambda_init,
                 int adaptive, int optimal, double scale, int verbose,
                 Array1D<double> &weights, Array1D<int> &used,
                 Array1D<double> &errbars, Array1D<double> &basis,
                 Array1D<double> &alpha, double &lambda) ;



#endif // BCS_H
