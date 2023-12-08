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
#ifndef MULTIINDEX_H
#define MULTIINDEX_H

#include "Array1D.h"
#include "Array2D.h"


/** \file multiindex.h
 * \brief Header for tools that deal with integer multiindices.
 * \todo Multiindex could be a separate class and a part of core UQTk.
 */

extern "C" void heap_ext_(const int *,const int *, int *, int *, int *);

/// \brief Computes the number of PC basis terms 
/// for Total-Order truncation with a given dimensionality and order
/// \note The formula is (ndim+norder)!/(ndim!norder!)  
int computeNPCTerms(int ndim,int norder);

/// \brief Computes the multiindex set of a PC basis 
/// for Total-Order truncation with a given dimensionality and order
/// Also, returns the number of terms.
int computeMultiIndex(int ndim,int norder, Array2D<int> &mi);

/// \brief Computes the multiindex set of a PC basis 
/// for Total-Order truncation with a given dimensionality and order
/// Also, returns the number of terms. Note that here, the
/// multiindex array pointer stores indices in column-major format,
/// i.e. mi[j*ndim+i] holds the j-th index for dimension i
int computeMultiIndexT(int ndim,int norder, int *mi);

/// \brief Computes the multiindex set of a PC basis 
/// for Total-Order truncation with a given dimensionality and order
/// Also, returns the number of terms.
int computeMultiIndex(int ndim,int norder, Array2D<int> &mi, string ordtype);


/// \brief Computes the multiindex set of a PC basis 
/// for Tensor-Product truncation with a given maximum order per dimensionality 
/// Also, returns the number of terms.
int computeMultiIndexTP(Array1D<int>& maxorders, Array2D<int>& mindex);

/// \brief Computes the number of PC basis terms 
/// for HDMR truncation with a given dimensionality and maxorders array 
/// that contains maximal orders per interaction dimensionalities.
int computeNPCTermsHDMR(int ndim, Array1D<int>& maxorders);

/// \brief Computes  the multiindex set of a PC basis 
/// for HDMR truncation with a given dimensionality and maxorders array 
/// that contains maximal orders per interaction dimensionalities.
int computeMultiIndexHDMR(int ndim, Array1D<int>& maxorders,Array2D<int>& mindex);

/// \brief Decode a multiindex set from a sparse format to a regular format
/// \note For encoding and for more details on the format, see encodeMindex function of PCSet class
/// \sa PCSet.h
void decodeMindex(Array1D< Array2D<int> >& sp_mindex, int ndim, Array2D<int>& mindex);


/// \brief Given a multiindex set it computes a new multiindex set where only 'admissible' bases are added
/// \note A new basis is admissible, if by subtracting one order from any of the dimensions with 
/// non-zero order, one never leaves the set of old multiindices
void upOrder(Array2D<int>& mindex,Array2D<int>& new_mindex);

/// \brief A boolean check to see if a new basis term is admissible or not
bool is_admis(Array1D<int>& mindex_try,Array2D<int>& mindex);

/// \brief Given a multiindex set, it returns the orders of each basis term
/// \note Essentially, this function performs sums of each rows
void getOrders(Array2D<int>& mindex,Array1D<int>& orders);

/// \brief Given a single multiindex, this returns its relative position in the total-order multiindex set
int get_invmindex(Array1D<int> mi);

/// \brief Given a single multiindex, this returns its relative position in the total-order multiindex set among the bases of the same order
int get_invmindex_ord(Array1D<int> mi);


#endif // MULTIINDEX_H
