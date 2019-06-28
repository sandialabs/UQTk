/* =====================================================================================
                     The UQ Toolkit (UQTk) version 3.0.4
                     Copyright (2017) Sandia Corporation
                     http://www.sandia.gov/UQToolkit/

     Copyright (2017) Sandia Corporation. Under the terms of Contract DE-AC04-94AL85000
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
#ifndef GKPLIB
#define GKPLIB

/** \file gkplib.h
 * Functions related to Gauss-Kronrod-Patterson sparse quadrature construction
 */

/// \brief retrieve pointers to 1D Clenshaw-Curtis rules
void getCC ( int n, int *nq, double **x, double **w );

/// \brief get order of Clenshaw-Curtis rules based on level
int getOrderCC ( int lev ) ;

/// \brief retrieve pointers to 1D Gauss-Kronrod-Patterson rules for
/// uniform pdf based on the quadrature level
void getGKPunif ( int n, int *nq, double **x, double **w );

/// \brief retrieve pointers to 1D Kronrod-Patterson rules for
/// normal pdf based on the quadrature level
void getGKPnorm ( int n, int *nq, double **x, double **w );

/// \brief get order of uniform Gauss-Kronrod-Patterson rules based on level
int getOrderGKPunif ( int lev ) ;

/// \brief get order of normal Gauss-Kronrod-Patterson rules based on level
int getOrderGKPnorm ( int lev ) ;

/// \brief List of decompositions of 'n' into 'dim' parts. The
/// implementation is based on Algorithm 5 of Combinatorial Algorithms
/// by Albert Nijenhuis, Herbert Wilf
void getCompNintoDim(int n, int dim, int *nelem, int **plist) ;

/// \brief Initial estimate for sparse grid size
int getSpgSize ( int getOrder ( int ), int dim, int lev );

/// \brief Sort sparse grid in lexicographical order
void sortSpg ( int dim, int spgSize, double *qpts, double *w );

/// \brief compute dim-dimensional tensor grid based a series of 1D rules
void getTensorProd(int dim, double *qpts, double *w, int *spgSize, int *n1D,
                   double **x1D, double **w1D, double qfac);

/// \brief Main function that connects the user setup for pdftype,
/// dimensionality, and quadrature level and various pieces of the
/// sparse quadrature algorithm employing Gauss-Kronrod-Patterson rules
void getSpgQW ( void get1DQW ( int , int *, double **, double** ), int getOrder ( int ),
		int dim, int lev, int *nqpts, double **qpts, double
		**w );

void getSpgAnisQW ( void get1DQW ( int , int *, double **, double** ), int getOrder ( int ),
		    int dim, int *levList, int *nqpts, double **qpts, double **w ) ;

void getCC ( int n, int *nq, double **x, double **w );
int getOrderCC ( int lev );

/// brief Fortran function for sorting an array of items. The array
/// operations happen outside this function, based on a series of
/// flags passed between the user code and this function. This
/// implementation is based on Algorithm 15 of Combinatorial Algorithms
/// by Albert Nijenhuis, Herbert Wilf
extern "C" void heap_ext_(const int *,const int *, int *, int *, int *);

#endif
