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
/// \file arraytools.h
/// \brief Header file for array tools
/// \todo Some functions are not optimal in terms of array access. 
/// \todo Some functions should be templated and or moved to array class

#ifndef ARRAYTOOLS_H
#define ARRAYTOOLS_H

#include <stdlib.h>
#include "Array1D.h"
#include "Array2D.h"

/// \brief Store a given 1d array in a 2d array with a single second dimension
template <typename T> void array1Dto2D(Array1D<T>& arr_1d,Array2D<T>& arr);

/// \brief Store a given 2d array with a single second dimension in a 1d array
template <typename T> void array2Dto1D(Array2D<T>& arr_2d,Array1D<T>& arr);

/// \brief Paste two 1d arrays of same size into a single 2d array with second dimension equal to two
template <typename T> void paste(Array1D<T>& arr1,Array1D<T>& arr2,Array2D<T>& arr);

/// \brief Generates multigrid as a cartesian product of each column of grid
/// \todo Should ideally be written in a recursive manner, similar to computeMultiIndexTP() in tools/multiindex.cpp
template <typename T> void generate_multigrid(Array2D<T>& multigrid,Array2D<T>& grid);

/// \brief Paste two 2D arrays next to each other (horizontal stack)
void paste(Array2D<double>& x, Array2D<double>& y, Array2D<double>& xy);

/// \brief Merge 2d double arrays (vertical stack)
void merge(Array2D<double>& x, Array2D<double>& y, Array2D<double>& xy);
/// \brief Merge 1d double arrays
void merge(Array1D<double>& x, Array1D<double>& y, Array1D<double>& xy);
/// \brief Merge 1d int arrays
void merge(Array1D<int>&    x, Array1D<int>&    y, Array1D<int>&    xy);

/// \brief Append array y to array x in place (double format)
void append(Array1D<double>& x, Array1D<double>& y);
/// \brief Append array y to array x in place (int format)
void append(Array1D<int>&    x, Array1D<int>&    y);

/// \brief Transpose a 2d double or int array x and return the result in xt
template <typename T> void transpose(Array2D<T> &x, Array2D<T> &xt);

/// \brief Unfold/flatten a 2d array into a 1d array (double format)
void flatten(Array2D<double>& arr_2, Array1D<double>& arr_1);

/// \brief Fold a 1d array into a 2d array (double format), row first
/// \note The dimension of the 1d array needs to be equal to
///  the product of the dimensions of the 2d array
void fold_1dto2d_rowfirst(Array1D<double>& x1, Array2D<double>& x2);

/// \brief Fold a 1d array into a 2d array (double format), column first
/// \note The dimension of the 1d array needs to be equal to
///  the product of the dimensions of the 2d array
void fold_1dto2d_colfirst(Array1D<double>& x1, Array2D<double>& x2);

/// \brief Swap i-th and j-th elements of the array arr
void swap(Array1D<double>& arr,int i,int j);

/// \brief Swap i-th and j-th rows of the 2d array arr
void swap(Array2D<double>& arr,int i,int j);

/// \brief Access element \f$j+i\times ny\f$ from 1D array 'arr_1'
double access(int nx, int ny, Array1D<double>& arr_1, int i, int j);

/// \brief Retrieves row 'k' from 2D array 'arr2d' and returns it in 1D array 'arr1d'
template <typename T> void getRow(Array2D<T> &arr2d, int k, Array1D<T> &arr1d);

/// \brief Retrieves column 'k' from 2D array 'arr2d' and returns it in 1D array 'arr1d'
template <typename T> void getCol(Array2D<T> &arr2d, int k, Array1D<T> &arr1d);

/// \brief Adds 'val' to the first n elements of an array pointer (double or int)
template <typename T> void addVal(int n, T *arr1d, T val) ;

/// \brief Adds 'val' to all elements of 1D array arr1d (double or int)
template <typename T> void addVal(Array1D<T> &arr1d, T val) ;

/// \brief Adds 'val' to all elements of 2D array arr2d (double or int)
template <typename T> void addVal(Array2D<T> &arr2d, T val) ;

/// \brief Extracts from 'vector', elements corresponding to indices 'ind' and returns them in 'subvector' (double or int)
template <typename T> void subVector(Array1D<T> &vector, Array1D<int> &ind, Array1D<T> &subvector);

/// \brief Extracts from 'matrix' rows corresponding to indices 'ind' and returns them in 'submatrix' (double or int)
template <typename T> void subMatrix_row(Array2D<T> &matrix, Array1D<int> &ind, Array2D<T> &submatrix);

/// \brief Extracts from 'matrix' columns corresponding to indices 'ind' and returns them in 'submatrix' (double or int)
template <typename T> void subMatrix_col(Array2D<T> &matrix, Array1D<int> &ind, Array2D<T> &submatrix);

/// \brief Adds scaled row or column to all rows / columns of a matrix (double or int)
/// \note RC is a character "R" or "C" for row or column, correspondingly
template <typename T> void matPvec(Array2D<T> &matrix, const Array1D<T> &rc, T alpha, char *RC);

/// \brief Returns maximum value in 'vector' and its location in *indx (double or int)
template <typename T> T maxVal(const Array1D<T>& vector, int *indx) ;

/// \brief Returns \f$ C=A\backslash B\f$ (C=Elements of A that are not in B); C is sorted in ascending order
void setdiff(Array1D<int> &A, Array1D<int> &B, Array1D<int> &C) ;

/// \brief Returns \f$ C=A\backslash B\f$ ( C=Elements of A that are not in B); C is sorted in ascending order
/// \note Assumes A is sorted and uses a faster algorithm than setdiff
/// \todo In future, this should sort A too and replace setdiff
/// \note B is sorted on output as well
void setdiff_s(Array1D<int> &A, Array1D<int> &B, Array1D<int> &C) ;

/// \brief Sorts integer array
void shell_sort (int *a, int n) ;
/// \brief Sorts integer array in ascending order
void shell_sort(Array1D<int>& array);
/// \brief Sorts double array in ascending order
void shell_sort(Array1D<double>& array);
/// \brief Sorts double array in ascending order according to a given column
void shell_sort_col(Array2D<double>& array,int col,Array1D<int>& newInd, Array1D<int>& oldInd);
/// \brief Sorts double array in ascending order according to first column, then second column breaks the tie, and so on
void shell_sort_all(Array2D<double>& array,Array1D<int>& newInd, Array1D<int>& oldInd);
/// \brief Quick-sort with 3-way partitioning of array between indices l and r
void quicksort3(Array1D<double>& arr, int l, int r);
/// \brief Quick-sort with 3-way partitioning of 2d array between indices l and r, according to column col
void quicksort3(Array2D<double>& arr,int left, int right,int col);
/// \brief Quick-sort with 3-way partitioning of 2d array between indices l and r, and sorting is done comparing rows (by first element, then by second, etc...)
void quicksort3(Array2D<double>& arr,int left, int right);

/// \brief Finds common entries in 1D arrays 'A' and 'B' and returns them in 'C', sorted in ascending order. It also
/// returns the original locations of these entries in 1D arrays 'iA' and 'iB', respectively
/// \note Currently, duplicated entries in either 'A' and 'B' will be duplicated in 'C'
void intersect(Array1D<int> &A, Array1D<int> &B, Array1D<int> &C, Array1D<int> &iA,Array1D<int> &iB) ;
/// \brief Find common entries in 1D arrays 'A' and 'B' and return them in 'C', sorted in ascending order
/// \note Currently, duplicated entries in either 'A' and 'B' will be duplicated in 'C'
void intersect(Array1D<int> &A, Array1D<int> &B, Array1D<int> &C) ;

/// \brief Return list of indices corresponding to elements of 1D array theta that are: larger ( type="gt" ),
/// larger or equal ( type="ge" ), smaller ( type="lt" ), smaller or equal ( type="le" ) than lmbda
template <typename T> void find(Array1D<T> &theta, T lmbda, string type, Array1D<int> &indx) ;

/// \brief Returns \f$y=\alpha Ax\f$, where 'A' is a \f$\left[n\times m\right]\f$ 2D array, 'x' is
/// 1D array of size \f$m\f$ and 'alpha' is a scalar. The 1D array 'y' has \f$n\f$ elements
void prodAlphaMatVec (Array2D<double>& A, Array1D<double>& x, double alpha, Array1D<double>& y) ;
/// \brief Returns \f$y=\alpha A^Tx\f$, where 'A' is a \f$\left[m\times n\right]\f$ 2D array, 'x' is
/// 1D array of size \f$m\f$ and 'alpha' is a scalar. The 1D array 'y' has \f$n\f$ elements
void prodAlphaMatTVec(Array2D<double>& A, Array1D<double>& x, double alpha, Array1D<double>& y) ;
/// \brief Returns \f$C=\alpha AB\f$, where 'A' and 'B' are \f$\left[m\times n\right]\f$ 2D arrays
/// and 'alpha' is a scalar. The 2D array 'C' has \f$m\times m\f$ elements
void prodAlphaMatMat(Array2D<double>& A, Array2D<double>& B, double alpha, Array2D<double>& C);
/// \brief Returns \f$C=\alpha A^TB\f$, where 'A' and 'B' are \f$\left[m\times n\right]\f$ 2D arrays
/// and 'alpha' is a scalar. The 2D array 'C' has \f$m\times m\f$ elements
void prodAlphaMatTMat(Array2D<double>& A, Array2D<double>& B, double alpha, Array2D<double>& C) ;
/// \brief Implements \f$x_i=x_i+\alpha y_i^ip\f$, where 'x' and 'y' are 1D arrays with \f$n\f$ elements
void addVecAlphaVecPow(Array1D<double>& x, double alpha, Array1D<double>& y, int ip) ;
/// \brief Returns \f$a^T B c\f$
double prod_vecTmatvec(Array1D<double>& a, Array2D<double>& B, Array1D<double>& c);
/// \brief Returns \f$A^T A\f$, where 'A' is a \f$\left[n\times k\right]\f$ 2D array
Array2D<double> MatTMat(Array2D<double>& A) ;


/// \brief Deletes row 'irow' from 2D array 'A'
/// \todo This should move to Array2D class
template <typename T> void delRow(Array2D<T>& A, int irow) ;

/// \brief Deletes column 'icol' from 2D array 'A'
/// \todo This should move to Array2D class
template <typename T> void delCol(Array2D<T> &A, int icol) ;

/// \brief Deletes element 'icol' from 1D array 'A'
/// \todo This should move to Array1D class
template <typename T> void delCol(Array1D<T> &x, int icol) ;

/// \brief Padds 2D array 'A' with the row 'x'
/// \note the number of elements in 'x' should be the same as the number of columns of 'A'
void paddMatRow(Array2D<double>& A, Array1D<double>& x) ;
/// \brief Padds 2D array 'A' with the column 'x'
/// \note the number of elements in 'x' should be the same as the number of rows in 'A'
void paddMatCol(Array2D<double>& A, Array1D<double>& x) ;
/// \brief Padds 2D array 'A' with the row 'x'
/// \note the number of elements in 'x' should be the same as the number of columns of 'A'
void paddMatRow(Array2D<int>& A, Array1D<int>& x) ;
/// \brief Padds 2D array 'A' with the column 'x'
/// \note the number of elements in 'x' should be the same as the number of rows in 'A'
void paddMatCol(Array2D<int>& A, Array1D<int>& x) ;
/// \brief Padds square 2D array 'A' \f$\left[n\times n\right]\f$ with the elements of 'x' and 'scal' as follows:
/// \f$A_{n+1,i}=A_{i,n+1}=x_i\f$ and \f$A_{n+1,n+1}=scal\f$
void paddMatColScal(Array2D<double>& A, Array1D<double>& x, double scal) ;

/// \brief Checks if two 1d int arrays are equal
bool is_equal(Array1D<int>& a, Array1D<int>& b);
/// \brief Checks if two 1d double arrays are equal
bool is_equal(Array1D<double>& a, Array1D<double>& b);
/// \brief Checks if one 1d int array is less than another (by first element, then by second, etc...)
bool is_less(Array1D<int>& a, Array1D<int>& b);
/// \brief Checks if one 1d double array is less than another (by first element, then by second, etc...)
bool is_less(Array1D<double>& a, Array1D<double>& b);

/// \brief Checks if vec matches with any of the rows of array
/// Returns the row number, or -1 if vec is not equal to any of the rows of array
int vecIsInArray(Array1D<int>& vec, Array2D<int>& array);

/// \brief Select the k-th smallest element of an array arr
double select_kth(int k, Array1D<double>& arr);

/// \brief Log-determinant of a real symmetric positive-definite matrix
/// \todo Check and catch the symmetric and positiv-definite conditions.
double logdeterm(Array2D<double>& mat);

/// \brief Trace of a matrix
double trace(Array2D<double>& mat);

/// \brief Evaluates the natural logarithm of a multivariate normal distribution
double evalLogMVN(Array1D<double>& x,Array1D<double>& mu,Array2D<double>& Sigma);

/// \brief Returns a diagonal matrix with a given diagonal
Array2D<double> diag(Array1D<double>& diagonal_array);
/**********************************************************
NEW ROUTINES - Kenny
***********************************************************/

/// \brief Returns a copy of 1D array
Array1D<double> copy(Array1D<double>&);

/// \brief Return a copy of 2D Array
Array2D<double> copy(Array2D<double>&);

/// \brief Deletes matrix columns or rows. Index specifies which column or row and dim = 1 deletes column, dim = 0 deletes the row
Array2D<double> mtxdel(Array2D<double>&, int index, int dim);

/// \brief Add two 1D Arrays and returns sum (must be of the same shape)
Array1D<double> add(Array1D<double>&, Array1D<double>&);

/// \brief Add two 2D Arrays and returns sum (must be of same shape)
Array2D<double> add(Array2D<double>&, Array2D<double>&);

/// \brief Add two 2D Arrays in place. Summation is returned as x.
void addinplace(Array2D<double>& x, Array2D<double>& y);

/// \brief Add two 1D Arrays in place. Summation is returned as x.
void addinplace(Array1D<double>& x, Array1D<double>& y);

/// \brief Returns subtraction of two 1D Arrays (must be of the same shape)
Array1D<double> subtract(Array1D<double>&, Array1D<double>&);

/// \brief Returns subtraction of two 2D Arrays (must be of the same shape)
Array2D<double> subtract(Array2D<double>&, Array2D<double>&);

/// \brief Subtract two 2D Arrays in place. Difference is returned as x.
void subtractinplace(Array2D<double>& x, Array2D<double>& y);

/// \brief Subtract two 1D Arrays in place. Difference is returned as x.
void subtractinplace(Array1D<double>& x, Array1D<double>& y);

/// \brief Returns 1D Arrays scaled by a double
Array1D<double> scale(Array1D<double>&, double);

/// \brief Returns 2D Array scaled by a double
Array2D<double> scale(Array2D<double>&, double);

/// \brief Multiply Array1D by double in place
void scaleinplace(Array1D<double>&, double);

/// \brief Multiply Array1D by int in place
void scaleinplace(Array1D<int>&, int);

/// \brief Multiply Array2D by double in place
void scaleinplace(Array2D<double>&, double);

/// \brief Multiply Array2D by int in place
void scaleinplace(Array2D<int>&, int);

/// \brief Returns the elementwise multiplication of two 2D Arrays
Array2D<double> dotmult(Array2D<double>&A,Array2D<double>&B );

/// \brief Returns the elementwise multiplication of two 1D Arrays
Array1D<double> dotmult(Array1D<double>&A,Array1D<double>&B );

/// \brief Returns the elementwise division of two 2D Arrays
Array2D<double> dotdivide(Array2D<double>&A,Array2D<double>&B );

/// \brief Returns the elementwise division of two 1D Arrays
Array1D<double> dotdivide(Array1D<double>&A,Array1D<double>&B );

/// \brief Returns norm of 1D Array (Euclidean)
double norm(Array1D<double>&);

/// \brief Weighted vector distance-squared
double dist_sq(Array1D<double>& x, Array1D<double>& y, Array1D<double>& w);

/// \brief Returns the transpose of a 2D Array
Array2D<double> Trans(Array2D<double>&);

/// \brief Returns the dot product of two 1D Arrays (must be of the same length)
double dot(Array1D<double>&, Array1D<double>&);

/// \brief Returns the matrix vector product
Array1D<double> dot(Array2D<double>&, Array1D<double>&);

/// \brief Returns the matrix matrix product
Array2D<double> dot(Array2D<double>&, Array2D<double>&);

/// \brief Returns the matrix matrix^T product
Array2D<double> dotT(Array2D<double>&, Array2D<double>&);

/// \brief Returns the inverse of a square 2D Array
Array2D<double> INV(Array2D<double> &A);

/// \brief Solves linear system AX=H, i.e. returns  A^(-1)*H, where A is real, symmetric and positive definite
Array2D<double> AinvH(Array2D<double> &A,Array2D<double> &H);

/// \brief Solves linear system Ax=b, i.e. return A^(-1)*b where A is real, symmetric and positive definite
Array1D<double> Ainvb(Array2D<double> &A,Array1D<double> &b);

/// \brief Least squares solution for overdetermined system. Note that A must be "taller than wide". Solution is returned in x.
void LSTSQ(Array2D<double> &A, Array1D<double> &b, Array1D<double> &x);

/// \brief Computes the QR factorization of a 2D Array (need not be square)
void QR(Array2D<double>&B,Array2D<double>&Q,Array2D<double>&R);

/// \brief Computes the SVD calculation of a 2D Array (need not be square)
void SVD(Array2D<double>&A,Array2D<double>&U,Array1D<double>&S,Array2D<double>&VT);

/// \brief Prints 1D double Array to screen (alternative to for loop using cout)
void printarray(Array1D<double>&);

/// \brief Prints 1D int Array to screen (alternative to for loop using cout)
void printarray(Array1D<int>&);

/// \brief Prints 2D double Array to screen (alternative to for loop using cout)
void printarray(Array2D<double>&);

/// \brief Prints 2D int Array to screen (alternative to for loop using cout)
void printarray(Array2D<int>&);


//---------------------------------------------------------------------------------------
#endif // ARRAYTOOLS_H
