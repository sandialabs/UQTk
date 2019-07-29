/* =====================================================================================
                     The UQ Toolkit (UQTk) version @UQTKVERSION@
                     Copyright (@UQTKYEAR@) Sandia Corporation
                     http://www.sandia.gov/UQToolkit/

     Copyright (@UQTKYEAR@) Sandia Corporation. Under the terms of Contract DE-AC04-94AL85000
     with Sandia Corporation, the U.S. Government retains certain rights in this software.

     This file is part of The UQ Toolkit (UQTk)

     UQTk is free software: you can redistribute it and/or modify
     it under the terms of the GNU Lesser General Public License as published by
     the Free Software Foundation, either version 3 of the License, or
     (at your option) any later version.

     UQTk is distributed in the hope that it will be useful,
     but WITHOUT ANY WARRANTY; without even the implied warranty of
     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
     GNU Lesser General Public License for more details.

     You should have received a copy of the GNU Lesser General Public License
     along with UQTk.  If not, see <http://www.gnu.org/licenses/>.

     Questions? Contact Bert Debusschere <bjdebus@sandia.gov>
     Sandia National Laboratories, Livermore, CA, USA
===================================================================================== */
/** \file probability.h
 * \brief Header for probability and random number generation- related tools.
 * \todo There shuold be a RNG class as a part of core UQTk - most of these functions will fit there.
 */

#ifndef PROBABILITY_H
#define PROBABILITY_H

#include "Array1D.h"
#include "Array2D.h"

#define DSFMT_DO_NOT_USE_OLD_NAMES
#include "dsfmt_add.h"

#include "KCenterClustering.h"
#include "figtree_internal.h"


/// \brief An implementation of error function using incomplete gamma function
double erff(const double x);

/// \brief Inverse error function, input scaled to [-1,1]
/// \note Cephes Math Library Release 2.8:  June, 2000.
/// Copyright 1984, 1987, 1988, 1992, 2000 by Stephen L. Moshier
// The website (http://www.boutell.com/lsm/lsmbyid.cgi/000626) states copying policy=freely distributable as of July 2012
//  \note modified by Sandia UQTk group to scale the input to [-1,1]
double inverf(double y0);

/// \brief Inverse of the CDF of the normal random variable, uses inverf
double invnormcdf(double y);

/// \brief Normal random variable CDF
double normcdf(double y);

/// \brief Complementary function for normcdf
double normcdfc(double y);

/// \brief Generates a vector of i.i.d. uniform(0,1) random variable samples of size ns*nd, given integer seed
void generate_uniform(double* rvar,int ns, int nd, int zSeed);

/// \brief Generates a matrix of i.i.d. uniform(0,1) random variable samples, given integer seed
void generate_uniform(Array2D<double>& rvar,int zSeed);

/// \brief Generates a vector of i.i.d. uniform(0,1) random variable samples of size ns*nd, given pointer to the state of current random number generator
void generate_uniform(double *rvar, int ns, int nd, dsfmt_t *rnstate);

/// \brief Generates a matrix of i.i.d. uniform(0,1) random variable samples, given pointer to the state of current random number generator
void generate_uniform(Array2D<double> &rvar, dsfmt_t *rnstate);

/// \brief Generates a vector of i.i.d. uniform(0,1) random variable LHS samples of size ns*nd, given integer seed
void generate_uniform_lhs(double* rvar,int ns, int nd, int zSeed);

/// \brief Generates a matrix of i.i.d. uniform(0,1) random variable LHS samples, given integer seed
void generate_uniform_lhs(Array2D<double>& rvar,int zSeed);

/// \brief Generates a vector of i.i.d. uniform(0,1) random variable LHS samples of size ns*nd, given pointer to the state of current random number generator
void generate_uniform_lhs(double *rvar, int ns, int nd, dsfmt_t *rnstate);

/// \brief Generates a matrix of i.i.d. uniform(0,1) random variable LHS samples, given pointer to the state of current random number generator
void generate_uniform_lhs(Array2D<double> &rvar, dsfmt_t *rnstate);

/// \brief Generates a matrix of i.i.d. normal(0,1) random variable samples
void generate_normal(Array2D<double>& rvar,int zSeed);

/// \brief Generates a matrix of i.i.d. normal(0,1) random variable LHS samples
/// \todo LHS generation is far from optimal, it is quite slow
void generate_normal_lhs(Array2D<double>& rvar,int zSeed);

/// \brief Returns the median of a data array
double get_median(const Array1D<double>& data);

/// \brief Returns the mean of a 1D data array
double get_mean(const Array1D<double>& data);

/// \brief Returns the mean of a 2D data array
double get_mean(const Array2D<double>& data);

/// \brief Returns the std of a data array
double get_std(const Array1D<double>& data);

/// \brief Returns the std of a data array
double get_var(const Array1D<double>& data);

/// \brief Vector-mean and weighted variance
double getMean_Variance(Array2D<double>& data_c, Array1D<double>& w, Array1D<double>& mean);

/// \brief Vector mean, column by column
void getMean(Array2D<double>& data_c, Array1D<double>& mean);

/// \brief Vector mean, either column by column for RC="C" or
/// row by row for RC="R"
void getMean(Array2D<double>& data_c, Array1D<double>& mean, char *RC);

/// \brief  Random permutation of 0..n-1
void rperm(int n, int *a, dsfmt_t *rnstate);

/// \brief KDE estimation of a PDF
void getPdf_figtree(Array2D<double>& source,Array2D<double>& target,Array1D<double>& sig, Array1D<double>& density, Array1D<double>& weight);

/// \brief Compute the PDF of data at the given points using given
/// number of clusters (if ncl=0, then find the optimal cluster
/// number) and a scale factor for the optimal bandwidth
void getPdf_cl(Array2D<double>& data, Array2D<double>& points, Array1D<double>& dens,int ncl, double sfac);

/// \brief Compute a few standard covariance functions C(x_1,x_2)
double covariance(Array1D<double>& x1, Array1D<double>& x2,Array1D<double>& param, string covtype);


void ihsU(Array2D<double> &rndnos, int dfac, dsfmt_t *rnstate) ;

void ihsU(int ns, int np, double *rndnos, int dfac, dsfmt_t *rnstate) ;

void ihsP(int ns, int np, int *rpos, int dfac, dsfmt_t *rnstate);

/// \brief Compute distance correlation factors given a set of samples
/// (no. of rows in spl) from a collection of random variables (no. of
/// columns in spl). dCor(i,j), with i>j stores the distance
/// correlation values between random variables i and j
void distCorr(const Array2D<double> &spl, Array2D<double> &dCor) ;

//---------------------------------------------------------------------------------------
#endif // PROBABILITY_H
