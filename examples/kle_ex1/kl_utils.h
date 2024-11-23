/* =====================================================================================

                      The UQ Toolkit (UQTk) version 3.1.5
                          Copyright (2024) NTESS
                        https://www.sandia.gov/UQToolkit/
                        https://github.com/sandialabs/UQTk

     Copyright 2024 National Technology & Engineering Solutions of Sandia, LLC (NTESS).
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
#ifndef KLUTILS_H_Seen
#define KLUTILS_H_Seen

#include <cstdio>
#include <stddef.h>
#include <fstream>
#include <string>
#include <math.h>
#include <iostream>
#include "assert.h"
#include <getopt.h>

using namespace std;

#include "Array1D.h"
#include "Array2D.h"
#include "kle.h"
#include "dsfmt_add.h"
#include "tools.h"
#include "arrayio.h"
#include "arraytools.h"
#include "deplapack.h"
#include "error_handlers.h"

/// \brief Computes analytical covariance values for a 1D configuration
///
/// Arguments:
///   \li const double x,y: coordinates for points x and y
///   \li const double clen: correlation length
///   \li const double sigma: noise level
///   \li const string type: covariance type; "SqExp" for
///   square-exponential, "Exp" for absolute distance
double genCovAnl1D(const double x1, const double x2, const double clen,
                 const double sigma, const string type) ;

/// \brief Computes analytical covariance values for a 2D structured grid
///
/// Arguments:
///   \li const double Array2D<double> &xy: Array with 2D coordinates
///   \li const int i and j: indices of points in the grid
///   \li const double clen: correlation length
///   \li const double sigma: noise level
///   \li const string type: covariance type; "SqExp" for
///   square-exponential, "Exp" for absolute distance
double genCovAnl2D( const Array2D<double> &xy, const int &i, const int &j,
const double &clen, const double &sigma, const string &covtype);

/// \brief Computes mean of a row/column in a 2D array
///
/// Arguments:
///   \li Array2D<double> y: 2D array with sample values
///   \li const string idir: either "L" for row mean or "C" for column mean
///   \li const int ij: row/column index
double getMean(const Array2D<double> &y, const string &idir, const int &ij) ;

/// \brief Generates non-uniform grid in [0,L] clustered towards the
/// ends of the interval
///
/// Arguments:
///   \li const double xi: coordinate in [0,1]; this grid should be
///   uniformly spaced
///   \li const double a: define where clustering takes place, range
///   in [0,1]; $a=1/2$ for same clustering at 0 and L; higher/lower
///   values bias the clustering towards L or 0.
///   \li const double b: compression factor, $b>1$; higher
///   compression as b gets closer to 1
///   \li const double L: length of interval in the computational space
/// \note is returns a double, the location in [0,L] corresponding to
///   xi in [0,1]
double dcomp(const double xi, const double a, const double b, const double L) ;

/// \brief Generates non-uniform grid in [0,L] clustered towards 0
///
/// Arguments:
///   \li const double xi: coordinate in [0,1]; this grid should be
///   uniformly spaced
///   \li const double b: compression factor, $b>1$; higher
///   compression as b gets closer to 1
///   \li const double L: length of interval in the computational space
/// \note is returns a double, the location in [0,L] corresponding to
///   xi in [0,1]
double scomp(const double xi, const double b, const double L) ;

/// \brief Generates array with 1D grid locations
///
/// Arguments:
///   \li Array1D<double> xgrid: array holding the grid points in
///   ascending order
///   \li const int npts: no. of grid points
///   \li const double L: length of 1D range [0,L]
///   \li const char *type: grid type: "unif" (uniform), "cl0"
///   (clustered near 0), "cl0L" (clustered both near 0 and L)
///   \li const double a and b: parameters for the clustering algorithms
void genGrid1D(Array1D<double> &xgrid, const int npts, const double L, const char *type,
	       const double a, const double b);


/// \brief Generates arrays with coordinates for a 2D grid
///
/// Arguments:
///   \li Array1D<double> xgrid: array holding the x-coordinates
///   \li Array1D<double> ygrid: array holding the y-coordinates
///   \li const int nx and ny: no. of grid points in x and y directions
///   \li const double L: length of 2D range [0,L]x[0,L]
///   \li const char *type: grid type: "unif" (uniform), "cl0"
///   (clustered near 0), "cl0L" (clustered both near 0 and L)
///   \li const double a and b: parameters for the clustering algorithms
void genGrid2D(Array1D<double> &xgrid,Array1D<double> &ygrid,
	       const int nx, const int ny, const double L, const char *type,
	       const double a, const double b) ;

/// \brief Create a 1D equivalent grid so that the weights computed in
/// the KLE class are consistent with the 2D grid
/// Arguments:
///   \li Array1D<double> xgrid: array holding the x-coordinates
///   \li Array1D<double> ygrid: array holding the y-coordinates
///   \li Array1D<double> xg1d:  array holding the equivalent 1D coordinates
void getGrid1dEquiv(const Array1D<double> &xgrid, const Array1D<double> &ygrid, Array1D<double> &xg1d) ;

/// \brief Compute the weights for a two-dimensional unstructured
/// triangular grid discretization
/// Arguments:
///   \li const Array2D<double> xgrid: array holding the (x,y)-coordinates
/// of vertices
///   \li const Array2D<int> tgrid:    array holding triangle connectivities
/// of vertices
///   \li Array1D<double> w:           array holding the nodal weights
/// of vertices
void getWeights2DU(const Array2D<double> &xgrid, const Array2D<int> &tgrid, Array1D<double> &w) ;

/// \brief Compute the area of a triangle given the coordinates of its vertices
/// Arguments:
///   \li const double x1,x2,x3: x-coordinates of the three vertices
///   \li const double y1,y2,y3: y-coordinates of the three vertices
double trArea(const double x1,const double x2,const double x3,
              const double y1,const double y2,const double y3) ;
#endif
