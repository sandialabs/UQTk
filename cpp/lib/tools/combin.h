/* =====================================================================================

                      The UQ Toolkit (UQTk) version 3.1.3
                          Copyright (2023) NTESS
                        https://www.sandia.gov/UQToolkit/
                        https://github.com/sandialabs/UQTk

     Copyright 2023 National Technology & Engineering Solutions of Sandia, LLC (NTESS).
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
/// \file combin.h
/// \brief Header for combinatorial tools
/// \note Some functions are likely not optimal and could have been computed more efficiently.

#ifndef COMBIN_H
#define COMBIN_H

#include "Array2D.h"



/// \brief Calculates binomial coefficient C(n,k): n-choose-k
int choose(int n,int k);

/// \brief Calculates the factorial of a number
int factorial(int number);

/// \brief Calculates the logfactorial of a number
double logfactorial(int number);

/// \brief Computes all possible k-combinations of the first n non-negative integers 
/// and returns them in fullInd
void chooseComb(int n, int k,Array2D<int>& fullInd);


/// \brief Computes a random permutation of the first n non-negative integers
/// and returns is in perm 
void get_perm(int n, int* perm,int seed);

/// \brief Computes a random permutation of the first n non-negative integers
/// and returns is in perm 
/// \note n is the size of the array argument perm
void get_perm(Array1D<int>& perm, int seed);

 
/// \brief Compute the incomplete Gamma function with parameter a at point x
/// \note This is a slightly modified version of a code distributed by John Burkardt
/// \note see http://people.sc.fsu.edu/~jburkardt/cpp_src/asa147/asa147.html
/// \note see comments under the function definition
double gammai(const double a, const double x);

/// \brief Compute the Beta function at the point pair (z,w) 
double beta(const double z, const double w);

/// \brief Compute the incomplete Beta function with parameters a and b at point x
/// \note This is a slightly modified version of a code distributed by John Burkardt
/// \note see http://people.sc.fsu.edu/~jburkardt/cpp_src/asa063/asa063.html
/// \note see comments under the function file 
double betai(const double p, const double q, const double x);

/// \brief Computes the digamma, or psi, function, i.e. derivative of the logarithm of gamma function
/// \note This is a slightly modified version of a code distributed by John Burkardt
/// \note see http://people.sc.fsu.edu/~jburkardt/cpp_src/asa103/asa103.cpp
double digama ( double x );

/// \brief K-center clustering of data
/// \param[in]  data_in       : Nxd matrix of data
/// \param[in]  w             : Array of size d; dimension-wise scaling weights
/// \param[in]  ncl           : Number of clusters
/// \param[out] numData       : Array of size ncl; stores the number of elements for each cluster
/// \param[out] pClusterIndex : Array of size N indicating the cluster index for each data point
void clust(Array2D<double>& data_in, Array1D<double>& w,int ncl, Array1D<int>& numData,int *pClusterIndex);

/// \brief Multiple trials of K-center clustering and picking the best one according to explained variance criterion
double clust_best(Array2D<double>& data_in, Array1D<double>& w,int ncl, Array1D<int>& bestnumData,int *bestClusterIndex,int ntry);
/// \brief Find the best number of clusters in a dataset according to one of three (hardcoded) criteria
int findNumCl(Array2D<double>& data_in,Array1D<double>& w,int ntry);
//---------------------------------------------------------------------------------------
#endif // COMBIN_H
