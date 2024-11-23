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
/// \file arraytools.cpp
/// \brief Tools to manipulate Array 1D and 2D objects. Some tools mimick MATLAB functionalities.


#include "stdlib.h"
#include "stdio.h"
#include "math.h"
#include "assert.h"
#include <sstream>
#include <fstream>
#include <iomanip>

#include "arraytools.h"
#include "ftndefs.h"
#include "gen_defs.h"
#include "depblas.h"
#include "deplapack.h"

using namespace std;

// Store a given 1d array in a 2d array with a single second dimension
template <typename T>
void array1Dto2D(Array1D<T> &arr_1d, Array2D<T> &arr)
{
  int nd=arr_1d.XSize();
  arr.Resize(nd,1);

  for(int i=0;i<nd;i++)
    arr(i,0)=arr_1d(i);

  return;
}
template void array1Dto2D(Array1D<double> &arr_1d, Array2D<double> &arr);
template void array1Dto2D(Array1D<int>    &arr_1d, Array2D<int>    &arr);

// Store a given 2d array with a single second dimension in a 1d array
template <typename T>
void array2Dto1D(Array2D<T> &arr_2d, Array1D<T> &arr)
{
  int nd  = arr_2d.XSize();
  int one = arr_2d.YSize();

  // Size check
  CHECKEQ(one,1);

  arr.Resize(nd);

  for(int i=0;i<nd;i++)
    arr(i)=arr_2d(i,0);


  return;
}
template void array2Dto1D(Array2D<double> &arr_2d,Array1D<double> &arr);
template void array2Dto1D(Array2D<int>    &arr_2d,Array1D<int>    &arr);

// Paste two 1d arrays of same size into a single 2d array with second dimension equal to two
template <typename T>
void paste(Array1D<T>& arr1,Array1D<T>& arr2,Array2D<T>& arr)
{
  int n=arr1.XSize();
  int m=arr2.XSize();

  // Size check
  CHECKEQ(n,m);

  arr.Resize(n,2);
  for(int i=0;i<n;i++){
    arr(i,0)=arr1(i);
    arr(i,1)=arr2(i);
  }
  return;
}
template void paste(Array1D<double> &arr1, Array1D<double> &arr2, Array2D<double> &arr);
template void paste(Array1D<int>    &arr1, Array1D<int>    &arr2, Array2D<int>    &arr);

// Generates multigrid as a cartesian product of each column of grid
template <typename T>
void generate_multigrid(Array2D<T>& multigrid,Array2D<T>& grid)
{
  int ngrid=grid.XSize();
  int ndim=grid.YSize();

  int totgrid= (int) pow(ngrid,ndim);
  multigrid.Resize(totgrid,ndim);

  // Work arrays
  Array2D<int> pIndex(totgrid,ndim,0);

  // Indexing
  for (int it=0;it<totgrid;it++){
      pIndex(it,ndim-1)=(int) it/pow(ngrid,ndim-1);
      double dnum=0.0;
      for(int jdim=ndim-2; jdim>-1;jdim--){
        dnum = dnum + pIndex(it,jdim+1)*pow(ngrid,jdim+1);
        pIndex(it,jdim) =(int) (it-dnum)/pow(ngrid,jdim);
      }
    }

  // Fill in the array
  for (int it=0;it<totgrid;it++){
    for(int idim=0; idim<ndim;idim++){
      multigrid(it,idim)=grid(pIndex(it,idim),idim);
    }
  }

  return;
}
template void generate_multigrid(Array2D<double> &multigrid, Array2D<double> &grid);
template void generate_multigrid(Array2D<int>    &multigrid, Array2D<int>    &grid);

// Paste two 2D arrays next to each other (horizontal stack)
void paste(Array2D<double>& x, Array2D<double>& y, Array2D<double>& xy)
{

  int nsample = x.XSize();

  // Size check
  CHECKEQ(nsample, (int) y.XSize());


  int ndim1   = x.YSize();
  int ndim2   = y.YSize();
  int ndim    = ndim1+ndim2;

  int ntot1   = nsample*ndim1 ;
  int ntot2   = nsample*ndim2 ;
  int incr    = 1 ;

  xy.Resize(nsample,ndim);
  FTN_NAME(dcopy)(&ntot1,x.GetArrayPointer(),&incr,xy.GetArrayPointer(),      &incr);
  FTN_NAME(dcopy)(&ntot2,y.GetArrayPointer(),&incr,xy.GetArrayPointer()+ntot1,&incr);

  //// Older implementation
  // for (int i=0;i<nsample;i++){
  //   for(int idim=0;idim<ndim1;idim++){
  //     xy(i,idim)=x(i,idim);
  //   }
  // }

  // for (int i=0;i<nsample;i++){
  //   for(int idim=ndim1;idim<ndim;idim++){
  //     xy(i,idim)=y(i,idim-ndim1);
  //   }
  // }

  return;

}

// Merge 2d double arrays  (vertical stack)
void merge(Array2D<double>& x, Array2D<double>& y, Array2D<double>& xy)
{
  int ndim=x.YSize();

  // Size check
  CHECKEQ(ndim,(int) y.YSize());

  int nsample1=x.XSize();
  int nsample2=y.XSize();
  int  nsample=nsample1+nsample2;
  xy.Resize(nsample,ndim);


  for (int i=0;i<nsample1;i++){
    for(int idim=0;idim<ndim;idim++){
      xy(i,idim)=x(i,idim);
    }
  }

  for (int i=nsample1;i<nsample;i++){
    for(int idim=0;idim<ndim;idim++){
       xy(i,idim)=y(i-nsample1,idim);
    }
  }

  return;
}

// Merge 1d double arrays
void merge(Array1D<double>& x, Array1D<double>& y, Array1D<double>& xy)
{
  int ns1 = x.XSize();
  int ns2 = y.XSize();

  int ns12 = ns1+ns2;
  int incr = 1;

  xy.Resize(ns12,0.e0);

  FTN_NAME(dcopy)(&ns1,x.GetArrayPointer(),&incr,xy.GetArrayPointer(),    &incr);
  FTN_NAME(dcopy)(&ns2,y.GetArrayPointer(),&incr,xy.GetArrayPointer()+ns1,&incr);

  // for (int i=0;i<ns1;i++){
  //     xy(i)=x(i);
  // }
  // for (int i=ns1;i<ns12;i++){
  //     xy(i)=y(i-ns1);
  // }

  return;
}

/// \brief Merge 1d int arrays
void merge(Array1D<int>& x, Array1D<int>& y, Array1D<int>& xy)
{
  int ns1=x.XSize();
  int ns2=y.XSize();

  int ns12=ns1+ns2;

  xy.Resize(ns12,0);

  for (int i=0;   i<ns1;  i++) xy(i)=x(i);
  for (int i=ns1; i<ns12; i++) xy(i)=y(i-ns1);

  return;
}

// Append array y to array x in place (double format)
void append(Array1D<double>& x, Array1D<double>& y)
{
  int ns1  = x.XSize();
  int ns2  = y.XSize();

  int ns12 = ns1+ns2;
  int incr = 1;

  x.Resize(ns12);
  //for (int i=ns1;i<ns12;i++) x(i)=y(i-ns1) ;
  FTN_NAME(dcopy)(&ns2,y.GetArrayPointer(),&incr,x.GetArrayPointer()+ns1,&incr);

  return;

}

// Append array y to array x in place (int format)
void append(Array1D<int>& x, Array1D<int>& y)
{
  int ns1 = x.XSize();
  int ns2 = y.XSize();
  int ns12= ns1+ns2;

  x.Resize(ns12);
  for (int i=ns1;i<ns12;i++) x(i)=y(i-ns1) ;

  return;
}

// Transpose a 2d double or int array x and return the result in xt
template <typename T>
void transpose(Array2D<T> &x, Array2D<T> &xt)
{
  int nx = x.XSize();
  int ny = x.YSize();

  xt.Resize(ny, nx);

  for (int ix=0; ix<nx; ix++){
    for (int iy=0; iy<ny; iy++){
      xt(iy,ix) = x(ix,iy);
    }
  }

  return;

}
template void transpose(Array2D<double> &x, Array2D<double> &xt);
template void transpose(Array2D<int>    &x, Array2D<int>    &xt);

// Unfold/flatten a 2d array into a 1d array (double format)
void flatten(Array2D<double>& arr_2, Array1D<double>& arr_1)
{
    int nx=arr_2.XSize();
    int ny=arr_2.YSize();

    int nxy=nx*ny;

    arr_1.Resize(nxy,0.e0);

    for(int i=0;i<nx;i++){
      for(int j=0;j<ny;j++){
        arr_1(j+i*ny)=arr_2(i,j);
      }
    }

    return;
}

// Fold a 1d array into a 2d array (double format), row first
void fold_1dto2d_rowfirst(Array1D<double>& x1, Array2D<double>& x2)
{
  int nx=x2.XSize();
  int ny=x2.YSize();
  int nxy=nx*ny;

  // Size check
  CHECKEQ(nxy,(int) x1.XSize());

  for(int i=0;i<nx;i++){
    for(int j=0;j<ny;j++){
      x2(i,j)=x1(i*ny+j);
    }
  }

  return;
}

// Fold a 1d array into a 2d array (double format), column first
void fold_1dto2d_colfirst(Array1D<double>& x1, Array2D<double>& x2)
{
  int nx=x2.XSize();
  int ny=x2.YSize();
  int nxy=nx*ny;

  // Size check
  CHECKEQ(nxy,(int) x1.XSize());

  for(int i=0;i<nx;i++){
    for(int j=0;j<ny;j++){
      x2(i,j)=x1(j*nx+i);
    }
  }

  return;
}


// Swap i-th and j-th elements of the array arr
void swap(Array1D<double>& arr,int i,int j)
{

  double h=arr(i);
  arr(i)=arr(j);
  arr(j)=h;

  return;
}

// Swap i-th and j-th rows of the 2d array arr
void swap(Array2D<double>& arr,int i,int j)
{
    int n=arr.YSize();
    double tmp;
    for (int d=0;d<n;d++){
        tmp=arr(i,d);
        arr(i,d)=arr(j,d);
        arr(j,d)=tmp;
    }

    return;
}

// Access element \f$j+i\times ny\f$ from 1D array 'arr_1'
double access(int nx, int ny,Array1D<double>& arr_1, int i, int j)
{
  // Size check
  CHECKEQ(nx*ny,(int) arr_1.XSize());

  return arr_1(j+i*ny);
}

double accessPythonHelper(int nx, int ny,Array1D<double>& arr_1, int i, int j)
{
  return access(nx,ny,arr_1,i,j);
}

// Retrieves row 'k' from 2D array 'arr2d' and returns it in 1D array 'arr1d'
template <typename T>
void getRow(Array2D<T> &arr2d, int k, Array1D<T> &arr1d)
{

  arr1d.Clear();
  for(int i=0; i<(int)arr2d.YSize(); i++)
    arr1d.PushBack(arr2d(k,i));

  return;

}
template void getRow(Array2D<double> &arr2d, int k, Array1D<double> &arr1d);
template void getRow(Array2D<int>    &arr2d, int k, Array1D<int>    &arr1d);

// Retrieves column 'k' from 2D array 'arr2d' and returns it in 1D array 'arr1d'
template <typename T>
void getCol(Array2D<T>& arr2d, int k, Array1D<T>& arr1d)
{
  arr1d.Clear();
  for(int i=0;i<(int)arr2d.XSize();i++)
    arr1d.PushBack(arr2d(i,k));

  return;
}
template void getCol(Array2D<double> &arr2d, int k, Array1D<double> &arr1d);
template void getCol(Array2D<int>    &arr2d, int k, Array1D<int>    &arr1d);

// Adds 'val' to the first n elements of an array pointer (double or int)
template <typename T>
void addVal(int n, T *arr1d, T val)
{
  for(int i=0; i< n; i++) arr1d[i] += val;

  return;
}
template void addVal(int n, double *arr1d, double val);
template void addVal(int n, int    *arr1d, int    val);

// Adds 'val' to all elements of 1D array arr1d (double or int)
template <typename T>
void addVal(Array1D<T> &arr1d, T val) {
  for(int i=0; i< (int)arr1d.XSize(); i++) arr1d(i) += val;

  return;
}
template void addVal(Array1D<double> &arr1d, double val);
template void addVal(Array1D<int>    &arr1d, int    val);

// Adds 'val' to all elements of 2D array arr2d (double or int)
template <typename T>
void addVal(Array2D<T> &arr2d, T val) {
  for(int j=0; j< (int)arr2d.YSize(); j++)
    for(int i=0; i< (int)arr2d.XSize(); i++)
      arr2d(i,j) += val;
  return;
}
template void addVal(Array2D<double> &arr2d, double val);
template void addVal(Array2D<int>    &arr2d, int    val);

// Extracts from 'vector', elements corresponding to indices 'ind' and returns them in 'subvector' (double or int)
template <typename T>
void subVector(Array1D<T> &vector, Array1D<int> &ind, Array1D<T> &subvector)
{

  int n = vector.XSize();
  int k = ind.XSize();

  subvector.Resize(k, (T) 0);
  for (int ik=0; ik<k; ik++){
    if ( ind(ik)<0 || ind(ik)>=n ){
      printf("subVector()::Index ind(%d)=%d is not allowed. Exiting.\n",ik,ind(ik));
      exit(1);
    }
    subvector(ik)=vector(ind(ik));

  }

  return;

}
template void subVector(Array1D<double> &vector, Array1D<int> &ind, Array1D<double> &subvector);
template void subVector(Array1D<int>    &vector, Array1D<int> &ind, Array1D<int>    &subvector);

// Extracts from 'matrix' rows corresponding to indices 'ind' and returns them in 'submatrix' (double or int)
template <typename T>
void subMatrix_row(Array2D<T> &matrix, Array1D<int> &ind, Array2D<T> &submatrix) {
  int n=matrix.XSize();
  int m=matrix.YSize();
  int k=ind.XSize();

  submatrix.Resize(k,m, (T) 0);
  for (int ik=0;ik<k;ik++){
    if (ind(ik)<0 || ind(ik)>=n){
      printf("subMatrix()::Index ind(%d)=%d is not allowed. Exiting.\n",ik,ind(ik));
      exit(1);
    }
    for (int im=0;im<m;im++){
      submatrix(ik,im)=matrix(ind(ik),im);
    }
  }

  return;

}
template void subMatrix_row(Array2D<double> &matrix, Array1D<int> &ind, Array2D<double> &submatrix);
template void subMatrix_row(Array2D<int>    &matrix, Array1D<int> &ind, Array2D<int>    &submatrix);

// Extracts from 'matrix' columns corresponding to indices 'ind' and returns them in 'submatrix' (double or int)
template <typename T>
void subMatrix_col(Array2D<T> &matrix, Array1D<int> &ind, Array2D<T> &submatrix) {

  int n=matrix.XSize();
  int m=matrix.YSize();
  int k=ind.XSize();

  submatrix.Resize(n,k, (T) 0);
  for (int ik=0;ik<k;ik++){
    if (ind(ik)<0 || ind(ik)>=m){
      printf("subMatrix()::Index ind(%d)=%d is not allowed. Exiting.\n",ik,ind(ik));
      exit(1);
    }
    for (int in=0;in<n;in++){
      submatrix(in,ik)=matrix(in,ind(ik));
    }
  }

  return;

}
template void subMatrix_col(Array2D<double> &matrix, Array1D<int> &ind, Array2D<double> &submatrix);
template void subMatrix_col(Array2D<int>    &matrix, Array1D<int> &ind, Array2D<int>    &submatrix);

// Adds scaled row or column to all rows / columns of a matrix (double or int)
template <typename T>
void matPvec(Array2D<T> &matrix, const Array1D<T> &rc, T alpha, char *RC) {

  int nrows = matrix.XSize();
  int ncols = matrix.YSize();

  if ( std::string(RC) == std::string("C")) {
    if (nrows != (int) rc.Length()) {
      cout << "arraytools.cpp::matPvec(): Mismatch in array sizes: "<<nrows<<" vs "<<rc.XSize()
           << endl << flush;
      exit(1);
    }
  } else if ( std::string(RC) == std::string("R")) {
    if (ncols != (int) rc.Length()) {
      cout << "arraytools.cpp::matPvec(): Mismatch in array sizes: "<<ncols<<" vs "<<rc.XSize()
           << endl << flush;
      exit(1);
    }
  } else {
    cout << "arraytools.cpp:matPvec(): unknown flag. Please use R or C"
         << RC << endl << flush;
    exit(1);

  }

  Array1D<T> rcLoc = rc;
  if ( (double) alpha-1.0 != 0.0 ) {
    scaleinplace(rcLoc,alpha);
  }

  if ( std::string(RC) == std::string("C")) {
    for (int i2=0; i2<ncols; i2++)
      for (int i1=0; i1<nrows; i1++)
        matrix(i1,i2) += rcLoc(i2);
  } else {
    for (int i2=0; i2<ncols; i2++)
      for (int i1=0; i1<nrows; i1++)
        matrix(i1,i2) += rcLoc(i1);
  }
  return;

}
template void matPvec(Array2D<double> &matrix, const Array1D<double> &rc, double alpha, char *RC);
template void matPvec(Array2D<int>    &matrix, const Array1D<int>    &rc, int    alpha, char *RC);

// Returns maximum value in 'vector' and its location in *indx (double or int)
template <typename T>
T maxVal(const Array1D<T> &vector, int *indx)
{

  T maxVal_ = vector(0);
  (*indx) = 0 ;
  for(int i=1; i < (int) vector.XSize(); i++)
    if (vector(i) > maxVal_) {
      maxVal_ = vector(i);
      (*indx) = i ;
    }
  return maxVal_;

}
template double maxVal(const Array1D<double> &vector, int *indx);
template int    maxVal(const Array1D<int>    &vector, int *indx);

// Returns in C elements of A that are not in B; C is sorted in ascending order
void setdiff(Array1D<int> &A, Array1D<int> &B, Array1D<int> &C)
{
    C.Clear() ;
    bool fnd;
    for ( int i = 0; i < (int) A.XSize() ; i++ )
    {
        fnd = false;
        for ( int j = 0; j < (int) B.XSize() ; j++ )
            if ( A(i) == B(j) ) fnd = true ;
        if ( !fnd) C.PushBack(A(i));
    }

    /* order C in ascending order */
    shell_sort(C);
    return ;

}

// Returns in C elements of A that are not in B; C is sorted in ascending order
// Assumes A is sorted
void setdiff_s(Array1D<int> &A, Array1D<int> &B, Array1D<int> &C)
{
  shell_sort(B);

  C.Clear() ;
  int j=0;

  for ( int i = 0; i < (int) A.XSize() ; i++ )
  {
    while(A(i)>B(j)){
      j++;
    }
    if ( A(i) < B(j) )
      C.PushBack(A(i));
  }

  return ;

}

// Sorts integer array
void shell_sort (int *a, int n) {

  int j ;
  for (int h = n/2; h>0; h = h/2) {
    for ( int i = h; i < n; i++) {

      int k = a[i];
      for ( j = i; j >= h && k < a[j - h]; j -= h)
        a[j] = a[j - h];

      a[j] = k;

    }
  }

  return ;

}

// Sorts integer array in ascending order
void shell_sort(Array1D<int>& array)
{
   int flag = 1, length = array.XSize(), i;
   int temp;
   int d=length;
   while( flag || (d>1)){      // boolean flag (true when not equal to 0)
      flag = 0;           // reset flag to 0 to check for future swaps
      d = (d+1) / 2;
      for (i = 0; i < (length - d); i++){
        if (array(i + d) < array(i)){
          temp = array(i + d);      // swap items at positions i+d and d
          array(i + d)= array(i);
          array(i) = temp;
          flag = 1;                  // indicate that a swap has occurred
        }
      }
   }
   return;
}

// Sorts double array in ascending order
void shell_sort(Array1D<double>& array)
{
   int flag = 1, length = array.XSize(), i;
   double temp;
   int d=length;
   while( flag || (d>1))      // boolean flag (true when not equal to 0)
   {
      flag = 0;           // reset flag to 0 to check for future swaps
      d = (d+1) / 2;
      for (i = 0; i < (length - d); i++){
        if (array(i + d) < array(i))
        {
          temp = array(i + d);      // swap items at positions i+d and d
          array(i + d)= array(i);
          array(i) = temp;
          flag = 1;                  // indicate that a swap has occurred
        }
      }
   }

   return;
}

// Sorts double array in ascending order according to a given column
void shell_sort_col(Array2D<double>& array,int col, Array1D<int>& newInd, Array1D<int>& oldInd)
{

   int flag = 1, length = array.XSize(), ncol=array.YSize(),i,j;
   double temp;
   int d=length, tmp;


   newInd.Resize(length,0);

  while( flag || (d>1)){     // boolean flag (true when not equal to 0)
    flag = 0;           // reset flag to 0 to check for future swaps
    d = (d+1) / 2;
    for (i = 0; i < (length - d); i++){
      if (array(i + d,col) < array(i,col))
      {
        for (j=0;j<ncol;j++){
          temp = array(i + d,j);      // swap items at positions i+d and d
          array(i + d,j)= array(i,j);
          array(i,j) = temp;
  		  }
    		newInd(oldInd(i+d))=i;
    		newInd(oldInd(i))=i+d;

    		tmp=oldInd(i);
    		oldInd(i)=oldInd(i+d);
    		oldInd(i+d)=tmp;
        flag = 1;                  // indicate that a swap has occurred
      }
    }
  }
     return;
}

// Sorts double array in ascending order according to first column, then second column breaks the tie, and so on
void shell_sort_all(Array2D<double>& array,Array1D<int>& newInd, Array1D<int>& oldInd)
{
  int flag = 1, length = array.XSize(), ncol=array.YSize(),i,j;
  double temp;
  int d=length, tmp;

  newInd.Resize(length,0);

  while( flag || (d>1)){      // boolean flag (true when not equal to 0)
    flag = 0;           // reset flag to 0 to check for future swaps
    d = (d+1) / 2;
    for (i = 0; i < (length - d); i++){
  	  bool swflag=false;
      for(int col=0;col<ncol;col++){
	      //if ( fabs(array(i + d,col) - array(i,col) ) > 1e-10 ){
        if (fabs(array(i + d,col) - array(i,col) ) > 1e-10 && array(i + d,col) < array(i,col)){
          swflag=true; break;
	      }
        else if (fabs(array(i + d,col) - array(i,col) ) > 1e-10 && array(i + d,col) > array(i,col)){
		      swflag=false; break;
	      }
  	    else {}
	    }

      if (swflag){
		    for (j=0;j<ncol;j++){
          temp = array(i + d,j);      // swap items at positions i+d and d
          array(i + d,j)= array(i,j);
          array(i,j) = temp;
		    }
    		newInd(oldInd(i+d))=i;
    		newInd(oldInd(i))=i+d;

    		tmp=oldInd(i);
    		oldInd(i)=oldInd(i+d);
    		oldInd(i+d)=tmp;
        flag = 1;                  // indicate that a swap has occurred
      }
    }
  }

  return;

}

// Quick-sort with 3-way partitioning of array between indices l and r
void quicksort3(Array1D<double>& arr, int l, int r)
{
    if (l>=r) return;
    int i = l-1, j = r, p = l-1, q = r;
    double v = arr(r);
    if (r <= l) return;
    for (;;)
    {
      while (arr(++i) < v) ;
      while (v < arr(--j)) if (j == l) break;
      if (i >= j) break;
      swap(arr,i,j);
      if (arr(i) == v) { p++; swap(arr,p,i); }
      if (v == arr(j)) { q--; swap(arr,j,q); }
    }
    swap(arr,i,r); j = i-1; i = i+1;
    for (int k = l; k < p; k++, j--) swap(arr,k,j);
    for (int k = r-1; k > q; k--, i++) swap(arr,i,k);
    quicksort3(arr, l, j);
    quicksort3(arr, i, r);

    return;
}

// Quick-sort with 3-way partitioning of 2d array between indices l and r, according to column col
void quicksort3(Array2D<double>& arr,int l, int r,int col)
{
    if (l>=r) return;
    int i = l-1, j = r, p = l-1, q = r;
    double v = arr(r,col);
    if (r <= l) return;
    for (;;)
    {
      while (arr(++i,col) < v) ;
      while (v < arr(--j,col)) if (j == l) break;
      if (i >= j) break;
      swap(arr,i,j);
      if (arr(i,col) == v) { p++; swap(arr,p,i); }
      if (v == arr(j,col)) { q--; swap(arr,j,q); }
    }
    swap(arr,i,r); j = i-1; i = i+1;
    for (int k = l; k < p; k++, j--) swap(arr,k,j);
    for (int k = r-1; k > q; k--, i++) swap(arr,i,k);

    quicksort3(arr, l, j,col);
    quicksort3(arr, i, r,col);
    return;

}

// Quick-sort with 3-way partitioning of 2d array between indices l and r, and sorting is done comparing rows (by first element, then by second, etc...)
void quicksort3(Array2D<double>& arr, int l, int r)
{
    if (l>=r) return;
    int i = l-1, j = r, p = l-1, q = r;
    Array1D<double> v;
    getRow(arr,r,v);
    if (r <= l) return;
    Array1D<double> arri,arrj;
    for (;;)
    {
        do{getRow(arr,++i,arri);} while (is_less(arri,v));
        do{getRow(arr,--j,arrj);  if (j == l) break;}while (is_less(v,arrj));
        if (i >= j) break;
        swap(arr,i,j);
        getRow(arr,i,arri);
        if (is_equal(arri,v)) { p++; swap(arr,p,i); }
        getRow(arr,j,arrj);
        if (is_equal(v,arrj)) { q--; swap(arr,j,q); }
    }
    swap(arr,i,r); j = i-1; i = i+1;
    for (int k = l; k < p; k++, j--) swap(arr,k,j);
    for (int k = r-1; k > q; k--, i++) swap(arr,i,k);

    quicksort3(arr, l, j);
    quicksort3(arr, i, r);

    return;
}

// Finds common entries in 1D arrays 'A' and 'B' and returns them in 'C', sorted in ascending order
// It also returns the original locations of these entries in 1D arrays 'iA' and 'iB', respectively
void intersect(Array1D<int> &A, Array1D<int> &B, Array1D<int> &C,Array1D<int> &iA,Array1D<int> &iB)
{
    C.Clear() ;
    iA.Clear() ;
    iB.Clear() ;
    for ( int i = 0; i < (int) A.XSize() ; i++ )
        for ( int j = 0; j < (int) B.XSize() ; j++ )
            if ( A(i) == B(j) )
            {
              C.PushBack(A(i));
              iA.PushBack(i);
              iB.PushBack(j);
            }

    /* order C in ascending order */
    bool chgOrd=true;
    while (chgOrd)
    {
        chgOrd=false;
        for ( int i = 0; i < (int) C.XSize()-1 ; i++ )
        {
            if (C(i)>C(i+1))
            {
                chgOrd=true;
                int itmp ;
                itmp = C(i);  C(i)  = C(i+1);  C(i+1) =itmp ;
                itmp = iA(i); iA(i) = iA(i+1); iA(i+1)=itmp ;
                itmp = iB(i); iB(i) = iB(i+1); iB(i+1)=itmp ;
            }
        }
    }
    return ;

}

// Finds common entries in 1D arrays 'A' and 'B' and returns them in 'C', sorted in ascending order
void intersect(Array1D<int> &A, Array1D<int> &B, Array1D<int> &C)
{
    C.Clear() ;
    for ( int i = 0; i < (int) A.XSize() ; i++ )
        for ( int j = 0; j < (int) B.XSize() ; j++ )
            if ( A(i) == B(j) )
                C.PushBack(A(i));

    /* order C in ascending order */
    bool chgOrd=true;
    while (chgOrd)
    {
        chgOrd=false;
        for ( int i = 0; i < (int) C.XSize()-1 ; i++ )
        {
          if (C(i)>C(i+1))
          {
            chgOrd=true;
            int itmp ;
            itmp = C(i);  C(i)  = C(i+1);  C(i+1) =itmp ;
          }
        }
    }
    return ;

}

// Return list of indices corresponding to elements of 1D array theta that are: larger ( type="gt" ),
// larger or equal ( type="ge" ), smaller ( type="lt" ), smaller or equal ( type="le" ) than lmbda
template <typename T>
void find(Array1D<T> &theta, T lmbda, string type, Array1D<int> &indx)
{
  indx.Clear();
  if ( type == "gt" ) {
    for ( int i = 0; i<(int) theta.XSize(); i++)
      if ( theta(i) > lmbda ) indx.PushBack(i) ;
    return ;
  }
  if ( type == "ge" ) {
    for ( int i = 0; i<(int) theta.XSize(); i++)
      if ( theta(i) >= lmbda ) indx.PushBack(i) ;
    return ;
  }
  if ( type == "lt" ) {
    for ( int i = 0; i<(int) theta.XSize(); i++)
      if ( theta(i) < lmbda ) indx.PushBack(i) ;
    return ;
  }
  if ( type == "le" ) {
    for ( int i = 0; i<(int) theta.XSize(); i++)
      if ( theta(i) <= lmbda ) indx.PushBack(i) ;
    return ;
  }
  if ( type == "eq" ) {
    for ( int i = 0; i<(int) theta.XSize(); i++)
      if ( theta(i) == lmbda ) indx.PushBack(i) ;
    return ;
  }
  return ;
}
template void find(Array1D<double> &theta, double lmbda, string type, Array1D<int> &indx);
template void find(Array1D<int>    &theta, int    lmbda, string type, Array1D<int> &indx);

// Implements y = a A x
void prodAlphaMatVec(Array2D<double>& A, Array1D<double>& x, double alpha, Array1D<double>& y)
{

    int n=A.XSize();
    int m=A.YSize();

    // Size check
    CHECKEQ(m, (int) x.XSize());

    y.Resize(n,0.e0);

    char trans='n';
    double beta=0.e0;
    int xinc=1;
    int yinc=1;
    FTN_NAME(dgemv)(&trans, &n, &m, &alpha, A.GetArrayPointer(), &n, x.GetArrayPointer(), &xinc,  &beta, y.GetArrayPointer(), &yinc );

   /* older implementation
    for (int i=0;i<n;i++)
    {
        for (int j=0;j<m;j++)
            y(i) += A(i,j)*x(j) ;
        y(i) *= alpha ;
    }
    */


    return;

}

// Implements y = a A^T x
void prodAlphaMatTVec(Array2D<double>& A, Array1D<double>& x, double alpha, Array1D<double>& y)
{

    int n=A.YSize();
    int m=A.XSize();

    // Size check
    CHECKEQ(m, (int) x.XSize());


    y.Resize(n,0.e0);

    char trans='t';
    double beta=0.e0;
    int xinc=1;
    int yinc=1;
    FTN_NAME(dgemv)(&trans, &m, &n, &alpha, A.GetArrayPointer(), &m, x.GetArrayPointer(), &xinc,  &beta, y.GetArrayPointer(), &yinc );

      /*
    for (int i=0;i<n;i++)
    {
        for (int j=0;j<m;j++)
            y(i) += A(j,i)*x(j) ;
        y(i) *= alpha ;
    }
      */
    return;

}

// Implements C = a A B
void prodAlphaMatMat(Array2D<double>& A, Array2D<double>& B, double alpha, Array2D<double>& C)
{

  int n=A.YSize();
  int m=A.XSize();
  int k=B.YSize();

  // Size check
  CHECKEQ(n, (int) B.XSize());


  C.Resize(m,k,0.e0);

  char transa='n';
  char transb='n';

  double beta=0.e0;


  FTN_NAME(dgemm)(&transa, &transb, &m, &k, &n, &alpha, A.GetArrayPointer(), &m, B.GetArrayPointer(), &n,  &beta, C.GetArrayPointer(), &m );

  /*
   for (int j=0;j<m;j++)
   for (int i=0;i<m;i++) {
   for ( int k=0; k<n; k++)
   C(i,j) += A(k,i)*B(k,j) ; // reversed indices for A
   C(i,j) *= alpha ;
   }
   */
  return;

}

// Implements C = a A^T B
void prodAlphaMatTMat(Array2D<double>& A, Array2D<double>& B, double alpha, Array2D<double>& C)
{

    int n=A.XSize();
    int m=A.YSize();
    int k=B.YSize();

    // Size check
    CHECKEQ(n, (int) B.XSize());


    C.Resize(m,k,0.e0);

    char transa='t';
    char transb='n';

    double beta=0.e0;


     FTN_NAME(dgemm)(&transa, &transb, &m, &k, &n, &alpha, A.GetArrayPointer(), &n, B.GetArrayPointer(), &n,  &beta, C.GetArrayPointer(), &m );

     /*
    for (int j=0;j<m;j++)
      for (int i=0;i<m;i++) {
        for ( int k=0; k<n; k++)
          C(i,j) += A(k,i)*B(k,j) ; // reversed indices for A
        C(i,j) *= alpha ;
      }
     */
    return;

}

// Implements x = x + a y**ip
void addVecAlphaVecPow(Array1D<double>& x, double alpha, Array1D<double>& y, int ip)
{
    // Size check
    CHECKEQ( (int) x.XSize(), (int) y.XSize() );

    for (int i = 0; i < (int) x.XSize(); i++)
        x(i) += alpha*pow(y(i),ip) ;

    return ;

}

// Returns a^T B c
double prod_vecTmatvec(Array1D<double>& a, Array2D<double>& B, Array1D<double>& c)
{
  double prod=0.e0;

  Array1D<double> tmp;
  prodAlphaMatTVec(B,a,1.0,tmp);
  prod=dot(tmp,c);

  return prod;
}

// Returns A^T A
Array2D<double> MatTMat(Array2D<double>& A)
{
  int n=A.XSize();
  int k=A.YSize();

  Array2D<double> B(k,k,0.e0);
  for (int i=0;i<k;i++){
    for (int j=0;j<=i;j++){
      for (int in=0;in<n;in++){
        B(i,j)+= A(in,i)*A(in,j);
      }
      B(j,i)=B(i,j);
    }
  }

  return B;
}

// Deletes a row from a matrix
template <typename T>
void delRow(Array2D<T>& A, int irow)
{

    int n = A.XSize() ;
    int m = A.YSize() ;

    if ( n <= 1 || m == 0 ) return ;

    Array2D<T> B(n-1,m) ;
    for ( int i = 0; i < irow; i++ )
        for ( int j = 0; j < m; j++)
            B(i,j) = A(i,j) ;
    for ( int i = irow+1; i < n; i++ )
        for ( int j = 0; j < m; j++)
            B(i-1,j) = A(i,j) ;

    A.Resize(n-1,m);
    for ( int i = 0; i < n-1; i++ )
        for ( int j = 0; j < m; j++)
            A(i,j) = B(i,j) ;

    return ;

}
template void delRow(Array2D<double> &A, int irow);
template void delRow(Array2D<int>    &A, int irow);

// Deletes a column from a matrix
template <typename T>
void delCol(Array2D<T> &A, int icol)
{

    int n = A.XSize() ;
    int m = A.YSize() ;

    if ( n == 0 || m <= 1 ) return ;

    Array2D<T> B(n,m-1) ;
    for ( int i = 0; i < n; i++ )
        for ( int j = 0; j < icol; j++)
            B(i,j) = A(i,j) ;
    for ( int i = 0; i < n; i++ )
        for ( int j = icol+1; j < m; j++)
            B(i,j-1) = A(i,j) ;

    A.Resize(n,m-1);
    for ( int i = 0; i < n; i++ )
        for ( int j = 0; j < m-1; j++)
            A(i,j) = B(i,j) ;

    return ;

}
template void delCol(Array2D<double> &A, int icol);
template void delCol(Array2D<int>    &A, int icol);

// Deletes an element from an array
template <typename T>
void delCol(Array1D<T> &x, int icol)
{

    int n = x.XSize() ;

    if ( n == 0 ) return ;

    Array1D<T> y(n-1) ;
    for ( int i = 0; i<icol; i++)   y(i ) = x(i) ;
    for ( int i = icol+1; i<n; i++) y(i-1) = x(i) ;

    x.Resize(n-1);
    for ( int i = 0; i<n-1; i++) x(i) = y(i) ;

    return ;

}
template void delCol(Array1D<double> &x, int icol) ;
template void delCol(Array1D<int>    &x, int icol) ;

// Padds 2D double array 'A' with the row 'x'
void paddMatRow(Array2D<double>& A, Array1D<double>& x)
{

    int n = A.XSize() ;
    int m = A.YSize() ;

    // Size check
    CHECKEQ(m, (int) x.XSize());

    Array2D<double> B ;
    B=A;
    A.Resize(n+1,m);
    for ( int i = 0; i < m; i++ )
    {
        for ( int j = 0; j < n; j++)
            A(j,i) = B(j,i) ;
        A(n,i) = x(i) ;
    }

    return ;

}

// Padds 2D array 'A' with the column 'x'
void paddMatCol(Array2D<double>& A, Array1D<double>& x)
{

    int n = A.XSize() ;
    int m = A.YSize() ;

    // Size check
    CHECKEQ(n, (int) x.XSize());


    Array2D<double> B ;
    B=A;
    A.Resize(n,m+1);
    for ( int i = 0; i < n; i++ )
    {
        for ( int j = 0; j < m; j++)
            A(i,j) = B(i,j) ;
        A(i,m) = x(i) ;
    }

    return ;

}

// Padds 2D int array 'A' with the row 'x'
void paddMatRow(Array2D<int>& A, Array1D<int>& x)
{

    int n = A.XSize() ;
    int m = A.YSize() ;

    // Size check
    CHECKEQ(m, (int) x.XSize());

    Array2D<int> B ;
    B=A;
    A.Resize(n+1,m);
    for ( int i = 0; i < m; i++ )
    {
        for ( int j = 0; j < n; j++)
            A(j,i) = B(j,i) ;
        A(n,i) = x(i) ;
    }

    return ;

}

// Padds 2D int array 'A' with the column 'x'
void paddMatCol(Array2D<int>& A, Array1D<int>& x)
{

    int n = A.XSize() ;
    int m = A.YSize() ;

    // Size check
    CHECKEQ(n, (int) x.XSize());



    Array2D<int> B ;
    B=A;
    A.Resize(n,m+1);
    for ( int i = 0; i < n; i++ )
    {
        for ( int j = 0; j < m; j++)
            A(i,j) = B(i,j) ;
        A(i,m) = x(i) ;
    }

    return ;

}

// Padds square matrix A with a row and column x, and adds an element A_{n+1,n+1} to obtain a larger square matrix
void paddMatColScal(Array2D<double>& A, Array1D<double>& x, double scal)
{

  int n = x.XSize() ;

  // Size check
  CHECKEQ(n, (int) A.YSize() );
  CHECKEQ(n, (int) A.XSize());


  Array2D<double> B ;
  B=A;
  A.Resize(n+1,n+1);
  for ( int i = 0; i < n; i++ )
  {
      for ( int j = 0; j < n; j++)
          A(i,j) = B(i,j) ;
      A(i,n) = A(n,i)=x(i) ;
  }
  A(n,n) = scal ;

  return ;

}

// Checks if two 1d int arrays are equal
bool is_equal(Array1D<int>& a, Array1D<int>& b){

    int n=a.XSize();
    int m=b.XSize();
    if (n!=m)
        return false;

    for (int i=0;i<n;i++){

        if (a(i)!=b(i))
            return false;

    }

    return true;
}

// Checks if two 1d double arrays are equal
bool is_equal(Array1D<double>& a, Array1D<double>& b)
{

  int n=a.XSize();
  int m=b.XSize();
  if (n!=m)
      return false;

  for (int i=0;i<n;i++){

      if (a(i)!=b(i))
 //     if (fabs(a(i)-b(i))>1.e-10)
        return false;

  }

  return true;
}

// Checks if one 1d int array is less than another (by first element, then by second, etc...)
bool is_less(Array1D<int>& a, Array1D<int>& b)
{

  int n=a.XSize();
  int m=b.XSize();

  // Size check
  CHECKEQ(n,m);

  for (int i=0;i<n;i++){

      if (a(i)>b(i))
          return false;
      else if (a(i)<b(i))
          return true;

  }

    return false;
}

// Checks if one 1d double array is less than another (by first element, then by second, etc...)
bool is_less(Array1D<double>& a, Array1D<double>& b)
{

    int n=a.XSize();
    int m=b.XSize();

    // Size check
    CHECKEQ(n,m);

    for (int i=0;i<n;i++){

      //  if (fabs(a(i)-b(i))>1.e-10)
        //{
        if (a(i)>b(i))
            return false;
        else if (a(i)<b(i))
            return true;

    //    }
    }


    return false;
}

// Checks if vec matches with any of the rows of a 2d array
int vecIsInArray(Array1D<int>& vec, Array2D<int>& array)
{
  int dd=array.XSize();
  int dim=vec.XSize();

  // Size check
  CHECKEQ(dim,(int) array.YSize());

  for(int i=0;i<dd;i++){
    int fl=0;
    for(int id=0;id<dim;id++){
      if (vec(id) == array(i,id))
        fl+=1;
      else break;
    }
    if (fl==dim)
      return i;
  }

  return -1;

}

// Select the k-th smallest element of an array arr
double select_kth(int k, Array1D<double>& arr)
{
  int i,ir,j,l,mid;
  double a;

  int n=arr.XSize();
  l=0;
  ir=n-1;
  for(;;){
    if (ir<=l+1){
      if (ir==l+1 && arr(ir)<arr(l))
	      swap(arr,l,ir);
      return arr(k);
    }
    else {
      mid=(l+ir) >> 1;
      swap(arr,mid,l+1);
      if (arr(l)>arr(ir))
	      swap(arr,l,ir);
      if (arr(l+1)>arr(ir))
	      swap(arr,l+1,ir);
      if (arr(l)>arr(l+1))
	      swap(arr,l,l+1);
      i=l+1;
      j=ir;
      a=arr(l+1);
      for(;;){
	      do i++; while (arr(i)<a);
	      do j--; while (arr(j)>a);
	      if (j<i) break;
	      swap(arr,i,j);
      }
      arr(l+1)=arr(j);
      arr(j)=a;
      if (j>=k) ir=j-1;
      if (j<=k) l=i;
    }
  }
}


// Log-determinant of a real symmetric positive-definite matrix
double logdeterm(Array2D<double>& mat)
{
    Array2D<double> A;
    A=mat;

    int nd=A.XSize();
    int chol_info=0;
    char lu='L';
    double logDet=0.0;

    // Cholesky factorization, done in-place

    //for(int i=0;i<nd;i++){
    //   for(int j=0;j<nd;j++)
    //       printf("%lg ",A(i,j));
    //    printf("\n");
    //}

    FTN_NAME(dpotrf)(&lu,&nd, A.GetArrayPointer(),&nd,&chol_info);

    // Catch the error in Cholesky factorization
    if (chol_info != 0 ) {
        printf("logdeterm():Error in Cholesky factorization, info=%d, printing the matrix below:\n", chol_info);

        for(int i=0;i<nd;i++){
            for(int j=0;j<nd;j++)
            printf("%lg ",A(i,j));
            printf("\n");
        }


        exit(1);
    }


    //get the log determinant
    for(int i = 0; i < nd; i++) logDet += 2*log(A(i,i));

    return logDet;
}

// Computes trace of a matrix
double trace(Array2D<double>& mat)
{

    int nd=mat.XSize();
    double trace = 0.0;

    for(int i = 0; i < nd; i++) trace += mat(i,i);

    return trace;
}

// Evaluates the natural logarithm of a multivariate normal distribution
double evalLogMVN(Array1D<double>& x,Array1D<double>& mu,Array2D<double>& Sigma)
{
  // Check that the dimesnions match
  if(Sigma.XSize() != Sigma.YSize())
    throw Tantrum((string) "Error in evalMVN: passed matrix is not square");

  // Check that the dimesnions match
  if(Sigma.XSize() != x.XSize())
    throw Tantrum((string) "Error in evalMVN: dimension mismatch in passed matrix and vector");

  // Check that the dimesnions match
  if(Sigma.XSize() != mu.XSize())
    throw Tantrum((string) "Error in evalMVN: dimension mismatch in passed matrix and vector");

  int dim = Sigma.XSize();

  // Compute pi
  double pi=4.0*atan(1.0);

  // Compute the inverse of the covariance matrix
  Array2D<double> invSigma(dim,dim);
  invSigma = INV(Sigma);

  // Compute the argument of the exponential
  Array1D<double> diff(dim,0.0);
  Array1D<double> matvec;
  for (int i=0; i<dim; i++)
    diff(i)= x(i)- mu(i);
  matvec= dot(invSigma, diff);
  double exparg= dot(diff, matvec);

  // Compute the logarithm of the pre-exponential factor
  //double detSigma = DET(Sigma);
    double logdetSigma=logdeterm(Sigma);
    double logPreexp= (-dim/2.0)*log(2.0*pi)- 0.5*logdetSigma; //0.5*log(detSigma);

    double logMVN= -0.5*exparg+ logPreexp;

  return logMVN;
}

// Returns a diagonal matrix with a given diagonal
Array2D<double> diag(Array1D<double>& diagonal_array)
{

  int nx=diagonal_array.Length();
  Array2D<double> diagonal_matrix(nx,nx,0.e0);
  for (int i=0;i<nx;i++)
    diagonal_matrix(i,i)=diagonal_array(i);

  return diagonal_matrix;
}


/**********************************************************
New LA Op - Kenny
***********************************************************/

// Returns a copy of an array
Array1D<double> copy(Array1D<double>& in_array){
  int n = in_array.Length();
  Array1D<double> out(n,0e0);
  out = in_array;
  return out;
}

// Returns a copy of an array
Array2D<double> copy(Array2D<double>& in_array){
  int m = in_array.XSize();
  int n = in_array.YSize();
  Array2D<double> out(m,n,0e0);
  out = in_array;
  return out;
}

// Deletes matrix columns or rows.
Array2D<double> mtxdel(Array2D<double>& A, int index, int dim){
  // deletes column when dim = 1
  // deletes row when dim = 0
  int n = A.XSize() ;
  int m = A.YSize() ;

  if ( n == 0 || m <= 1 ) return A;

  if (dim == 1){
    Array2D<double> B(n,m-1) ;
    for ( int i = 0; i < n; i++ )
      for ( int j = 0; j < index; j++)
        B(i,j) = A(i,j) ;
    for ( int i = 0; i < n; i++ )
      for ( int j = index+1; j < m; j++)
        B(i,j-1) = A(i,j) ;
    return B;
  }

  if (dim == 0){
    Array2D<double> B(n-1,m) ;
    for ( int i = 0; i < index; i++ )
      for ( int j = 0; j < m; j++)
        B(i,j) = A(i,j) ;
    for ( int i = index+1; i < n; i++ )
      for ( int j = 0; j < m; j++)
        B(i-1,j) = A(i,j) ;
    return B;
  }
}

// add two vectors
Array1D<double> add(Array1D<double>& x, Array1D<double>& y){

  int nx = x.Length();
  int ny = y.Length();

  if ( nx != ny )
  {
    printf("add() : Error : no. of elements in x and size of y are not the same : %d %d\n",
       nx,ny);
    exit(1);
  }

  Array1D<double> ytemp = copy(y);
  double alpha = 1;
  int incr = 1;
  FTN_NAME(daxpy)(&ny, &alpha, x.GetArrayPointer(), &incr, ytemp.GetArrayPointer(), &incr);

  return ytemp;
}

// add two matrices of the same size
Array2D<double> add(Array2D<double>& x, Array2D<double>& y){

  int nx = x.XSize();
  int ny = y.XSize();

  // Size check
  CHECKEQ(nx,ny);

  int nx2 = x.YSize();
  int ny2 = y.YSize();

  // Size check
  CHECKEQ(nx2,ny2);

  Array2D<double> temp(nx,nx2,0.0);
  for (int i = 0; i < nx; i++){
    for (int j = 0; j < nx2; j++){
      temp(i,j) = x(i,j) + y(i,j);
    }
  }

  return temp;
}

// add two matrices of the same size
void addinplace(Array2D<double>& x, Array2D<double>& y){

  int nx = x.XSize();
  int ny = y.XSize();

  // Size check
  CHECKEQ(nx,ny);

  int nx2 = x.YSize();
  int ny2 = y.YSize();

   // Size check
  CHECKEQ(nx2,ny2);


  // Array2D<double> temp(nx,ny,0.0);
  for (int i = 0; i < nx; i++){
    for (int j = 0; j < nx2; j++){
      x(i,j) = x(i,j) + y(i,j);
    }
  }

  return;
}

// add two vectors of the same size
void addinplace(Array1D<double>& x, Array1D<double>& y){

  int nx = x.XSize();
  int ny = y.XSize();

  // Size check
  CHECKEQ(nx,ny);

  for (int i = 0; i < nx; i++){
    x(i) = x(i) + y(i);
  }

  return;
}


// subtract two vectors
Array1D<double> subtract(Array1D<double>& x, Array1D<double>& y){

  int nx = x.Length();
  int ny = y.Length();

  // Size check
  CHECKEQ(nx,ny);

  Array1D<double> xtemp = copy(x);
  double alpha = -1.0;
  int incr = 1;
  FTN_NAME(daxpy)(&ny, &alpha, y.GetArrayPointer(), &incr, xtemp.GetArrayPointer(), &incr);


  return xtemp;
}

// subtract two matrices of the same size
Array2D<double> subtract(Array2D<double>& x, Array2D<double>& y){

  int nx = x.XSize();
  int ny = y.XSize();

  // Size check
  CHECKEQ(nx,ny);

  int nx2 = x.YSize();
  int ny2 = y.YSize();

  // Size check
  CHECKEQ(nx2,ny2);

  Array2D<double> temp(nx,nx2,0.0);
  for (int i = 0; i < nx; i++){
    for (int j = 0; j < nx2; j++){
      temp(i,j) = x(i,j) - y(i,j);
    }
  }

  return temp;
}


// add two matrices of the same size
void subtractinplace(Array2D<double>& x, Array2D<double>& y){

  int nx = x.XSize();
  int ny = y.XSize();

  // Size check
  CHECKEQ(nx,ny);

  int nx2 = x.YSize();
  int ny2 = y.YSize();

   // Size check
  CHECKEQ(nx2,ny2);


  // Array2D<double> temp(nx,ny,0.0);
  for (int i = 0; i < nx; i++){
    for (int j = 0; j < nx2; j++){
      x(i,j) = x(i,j) - y(i,j);
    }
  }

  return;
}

// add two vectors of the same size
void subtractinplace(Array1D<double>& x, Array1D<double>& y){

  int nx = x.XSize();
  int ny = y.XSize();

  // Size check
  CHECKEQ(nx,ny);

  for (int i = 0; i < nx; i++){
    x(i) = x(i) - y(i);
  }

  return;
}

// multiply Array1D by double
Array1D<double> scale(Array1D<double>& x, double alpha){

  int incr = 1;
  int n    = x.XSize();

  Array1D<double> temp = x;
  FTN_NAME(dscal)(&n, &alpha, temp.GetArrayPointer(), &incr);

  return temp;

}

// multiply Array2D by double
Array2D<double> scale(Array2D<double>& x, double alpha){

  int incr = 1;
  int n    = x.XSize();
  int m    = x.YSize();
  int nm    = n*m;

  Array2D<double> temp = x;
  FTN_NAME(dscal)(&nm, &alpha, temp.GetArrayPointer(), &incr);

  return temp;

}

// multiply Array1D by double, in place
void scaleinplace(Array1D<double>& x, double alpha){

  int n = x.Length();
  int incr = 1;
  FTN_NAME(dscal)(&n, &alpha, x.GetArrayPointer(), &incr);

  return ;

}

// multiply Array1D by int, in place
void scaleinplace(Array1D<int>& x, int alpha){

  int n = x.Length();
  for ( int i=0; i<n; i++ ) x(i) *= alpha;

  return ;

}

// multiply Array2D by double, in place
void scaleinplace(Array2D<double>& x, double alpha){

  int incr = 1;
  int n    = x.XSize();
  int m    = x.YSize();
  int nm   = n*m;

  FTN_NAME(dscal)(&nm, &alpha, x.GetArrayPointer(), &incr);
  return ;

}

// multiply Array2D by int, in place
void scaleinplace(Array2D<int>& x, int alpha){

  int incr = 1;
  int n    = x.XSize();
  int m    = x.YSize();
  for ( int j=0; j<m; j++ )
    for ( int i=0; i<n; i++ )
      x(i,j) *= alpha;
  return ;

}

// Returns the elementwise multiplication of two 2D Arrays
Array2D<double> dotmult(Array2D<double>& A, Array2D<double>& B){

  int n = A.XSize();
  int m = B.XSize();

  // Size check
  CHECKEQ(n,m);

  int n1 = A.YSize();
  int m1 = B.YSize();

  // Size check
  CHECKEQ(n1,m1);

  Array2D<double> C(n,n1,0.0);
  for (int i = 0; i < n; i++){
    for (int j = 0; j < n1; j++){
      C(i,j) = A(i,j) * B(i,j);
    }
  }

  return C;
}

// Returns the elementwise multiplication of two 1D Arrays
Array1D<double> dotmult(Array1D<double>& A, Array1D<double>& B){

  int n = A.XSize();
  int m = B.XSize();

  // Size check
  CHECKEQ(n,m);

  Array1D<double> C(n,0.0);
  for (int i = 0; i < n; i++){
      C(i) = A(i) * B(i);
  }

  return C;
}

// Returns the elementwise division of two 2D Arrays
Array2D<double> dotdivide(Array2D<double>& A, Array2D<double>& B){

  int n = A.XSize();
  int m = B.XSize();

  // Size check
  CHECKEQ(n,m);

  int n1 = A.YSize();
  int m1 = B.YSize();

  // Size check
  CHECKEQ(n1,m1);

  Array2D<double> C(n,n1,0.0);
  for (int i = 0; i < n; i++){
    for (int j = 0; j < n1; j++){
      C(i,j) = A(i,j) / B(i,j);
    }
  }

  return C;
}

// Returns the elementwise division of two 1D Arrays
Array1D<double> dotdivide(Array1D<double>& A, Array1D<double>& B){

  int n = A.XSize();
  int m = B.XSize();

  // Size check
  CHECKEQ(n,m);

  Array1D<double> C(n,0.0);
  for (int i = 0; i < n; i++){
      C(i) = A(i) / B(i);
  }

  return C;
}

// Get norm of array 1d
double norm(Array1D<double>& x){
  int n = x.Length();
  int incr = 1;
  return dnrm2_(&n, x.GetArrayPointer(), &incr);
}

// Returns weighted vector distance-squared
double dist_sq(Array1D<double>& x, Array1D<double>& y, Array1D<double>& w)
{
  double dist=0., t;
  int ndim=x.XSize();//=y.XSize();//=w.XSize();

  for (int idim=0;idim<ndim;idim++){
    t=(x(idim)-y(idim))/w(idim);
    dist += t * t;
  }

  return dist;
}

// Get transpose of 2d array
Array2D<double> Trans(Array2D<double> &A){

  int n=A.XSize();
  int m=A.YSize();

  Array2D<double> B(m,n,0e0);

  for (int i = 0; i < m; i++){
  for (int j = 0; j < n; j++){
    B(i,j) = A(j,i);
  }
  }
  return B;
}

// get dot product between two vectors
double dot(Array1D<double>& v1, Array1D<double>& v2){
  int n1 = v1.Length();
  int n2 = v2.Length();

  // Size check
  CHECKEQ(n1,n2);

  int incr = 1;
  return FTN_NAME(ddot)(&n1, v1.GetArrayPointer(), &incr, v2.GetArrayPointer(), &incr);
}

// get matrix vector product
Array1D<double> dot(Array2D<double>& A, Array1D<double>& x){
  int n=A.XSize();
  int m=A.YSize();

  // Size check
  CHECKEQ(m, (int) x.XSize());


  Array1D<double> y(n,0e0);

  char trans='n';
  double beta=0.e0;
  int xinc=1;
  int yinc=1;
  double alpha = 1;
  FTN_NAME(dgemv)(&trans, &n, &m, &alpha, A.GetArrayPointer(), &n, x.GetArrayPointer(), &xinc,  &beta, y.GetArrayPointer(), &yinc );

  return y;
}

// get matrix vector product
Array1D<double> dot(Array1D<double>& x, Array2D<double>& A){
  int n=A.XSize();
  int m=A.YSize();

  // Size check
  CHECKEQ(n, (int) x.XSize());


  Array1D<double> y(n,0e0);

  char trans='n';
  double beta=0.e0;
  int xinc=1;
  int yinc=1;
  double alpha = 1;
  FTN_NAME(dgemv)(&trans, &n, &m, &alpha, A.GetArrayPointer(), &n, x.GetArrayPointer(), &xinc,  &beta, y.GetArrayPointer(), &yinc );

  return y;
}

// get matrix matrix product
Array2D<double> dot(Array2D<double>& A, Array2D<double>& B){

  int m=A.XSize();
  int k=A.YSize();
  int n=B.YSize();

  // Size check
  CHECKEQ(k, (int) B.XSize());



  Array2D<double> C(m,n,0.e0);

  char transa='n';
  char transb='n';
  double beta = 0.e0;
  double alpha = 1.e0;


  FTN_NAME(dgemm)(&transa, &transb, &m, &n, &k, &alpha, A.GetArrayPointer(), &m, B.GetArrayPointer(), &k,  &beta, C.GetArrayPointer(), &m );


  return C;
}

// get matrix^T matrix product
Array2D<double> dotT(Array2D<double>& A, Array2D<double>& B){

  // transose A before doing matrix matrix operation
  int n=A.XSize();
  int m=A.YSize();
  int k=B.YSize();

  // Size check
  CHECKEQ(n, (int) B.XSize());

  Array2D<double> C(m,k,0.e0);

  char transa='t';
  char transb='n';
  double beta=0.e0;
  double alpha = 1;


  FTN_NAME(dgemm)(&transa, &transb, &m, &k, &n, &alpha, A.GetArrayPointer(), &n, B.GetArrayPointer(), &n,  &beta, C.GetArrayPointer(), &m );


  return C;
}

// inverse of a real square matrix
Array2D<double> INV(Array2D<double> &A){

  int m = A.XSize();
  int n = A.YSize();
  int LUinfo;
  Array1D<int> IPIV(m,0);

  Array2D<double> B = A;

  // Get pivot indices from LU factorization
  FTN_NAME(dgetrf)(&m, &n, B.GetArrayPointer(), &n, IPIV.GetArrayPointer(), &LUinfo);

  Array1D<double> work(min(m,n),0e0);
  int lwork = -1;
  int info;

  // get work space first
  FTN_NAME(dgetri)(&m, B.GetArrayPointer(), &m, IPIV.GetArrayPointer(), work.GetArrayPointer(), &lwork, &info);

  int nwork = work(0);
  lwork = -min(-1,-nwork); // find maximum
  work.Resize(lwork,0e0);

  // Get inverse
  FTN_NAME(dgetri)(&m, B.GetArrayPointer(), &m, IPIV.GetArrayPointer(), work.GetArrayPointer(), &lwork, &info);

  return B;
}

// Solving AX=H where A is real, symmetric and positive definite
Array2D<double> AinvH(Array2D<double> &A,Array2D<double> &H){

  int n = A.XSize();
  CHECKEQ(n,A.YSize());
  CHECKEQ(n,H.XSize());
  int k=H.YSize();
  int info;

  char uplo='L';
    Array2D<double> L = A;
  Array2D<double> X = H;

  FTN_NAME(dposv)( &uplo, &n,&k, L.GetArrayPointer(), &n, X.GetArrayPointer(), &n, &info  );

  return X;
}

// Solving Ax=b where A is real, symmetric and positive definite
Array1D<double> Ainvb(Array2D<double> &A,Array1D<double> &b){

  int n = A.XSize();
  CHECKEQ(n,A.YSize());
  CHECKEQ(n,b.XSize());
  int k=1;
  int info;

  char uplo='L';
    Array2D<double> L = A;
  Array1D<double> x = b;

  FTN_NAME(dposv)( &uplo, &n,&k, L.GetArrayPointer(), &n, x.GetArrayPointer(), &n, &info  );

  return x;
}
// Least squares solution for overdetermined system
// A must be "taller than wide"
void LSTSQ(Array2D<double> &A, Array1D<double> &b, Array1D<double> &x){

  int m = A.XSize();
  int n = A.YSize();

  x.Resize(n,0);

  Array2D<double> B = A;

  Array1D<double> work(min(m,n),0e0);
  int lwork = -1;
  int info;
  char trans = 'N';
  int nrhs = 1;

  // get work space
  FTN_NAME(dgels)(&trans, &m, &n, &nrhs, B.GetArrayPointer(), &m, b.GetArrayPointer(), &m, work.GetArrayPointer(), &lwork, &info);

  int nwork = work(0);
  lwork = -min(-1,-nwork); // find maximum
  work.Resize(lwork,0e0);

  // now get least squares solution
  FTN_NAME(dgels)(&trans, &m, &n, &nrhs, B.GetArrayPointer(), &m, b.GetArrayPointer(), &m, work.GetArrayPointer(), &lwork, &info);

  for (int i = 0; i < n; i++){
    x(i) = b(i);
  }
  return;
}

// QR factorization
void QR(Array2D<double>& B, Array2D<double>& Q, Array2D<double>& R ){


  Q = B;

  int m = Q.XSize();
  int n = Q.YSize();
  int r = m - n;

  // Q.Resize(m,m,0.0);
  R.Resize(m,n,0.0);

  // resize Q
  Array1D<double> z(m,0.0);
  if (r >= 1){
    for (int i = 0; i < r; i++){
      Q.insertCol(z,n+i);
    }
  }
  int nnew = Q.YSize();

  //*** GET R MATRIX

  Array1D<double> tau(min(m,nnew),1);
  Array1D<double> work(min(m,nnew),0e0);
  int lwork = -1;
  int info;

  // get work space first
  FTN_NAME(dgeqrf)(&m, &nnew, Q.GetArrayPointer(), &m, tau.GetArrayPointer(), work.GetArrayPointer(), &lwork, &info);

  int nwork = work(0);
  lwork = -min(-1,-nwork); // find maximum
  work.Resize(lwork,0e0);

  // now run real QR factorization
  FTN_NAME(dgeqrf)(&m, &nnew, Q.GetArrayPointer(), &m, tau.GetArrayPointer(), work.GetArrayPointer(), &lwork, &info);

  // copy upper diagonal contents of A to R
  for (int i = 0; i < n; i++){
    for (int j = i; j < n; j++){
      R(i,j) = Q(i,j);
    }
  }

  /**********************/
  int k = tau.Length();
  lwork = -1;

  FTN_NAME(dorgqr)(&m, &nnew, &k, Q.GetArrayPointer(), &m, tau.GetArrayPointer(), work.GetArrayPointer(), &lwork, &info);

  nwork = work(0);
  lwork = -min(-1,-nwork); // find maximum
  work.Resize(lwork,0e0);


  FTN_NAME(dorgqr)(&m, &nnew, &k, Q.GetArrayPointer(), &m, tau.GetArrayPointer(), work.GetArrayPointer(), &lwork, &info);

  return;
}

// SVD calculation
void SVD(Array2D<double>& A,Array2D<double>& U,Array1D<double>& S,Array2D<double>& VT){

    int m = A.XSize();
  int n = A.YSize();

  Array2D<double> B = A;

  U.Resize(m,m,0.0);
  S.Resize(min(m,n),0.0);
  VT.Resize(n,n,0.0);

  char jobu = 'A';
  char jobvt = 'A';
  int lwork = -1;
  int info;
  Array1D<double> work(min(m,n),0e0);

  // determine optimal length for work first
  // DGESVD( JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, INFO )
  FTN_NAME(dgesvd)(&jobu, &jobvt, &m, &n, B.GetArrayPointer(), &m, S.GetArrayPointer(), U.GetArrayPointer(), &m, VT.GetArrayPointer(), &n, work.GetArrayPointer(), &lwork, &info);

  int nwork = work(0);
  lwork = -min(-1,-nwork); // find maximum
  work.Resize(lwork,0e0);

  FTN_NAME(dgesvd)( &jobu, &jobvt, &m, &n, B.GetArrayPointer(), &m, S.GetArrayPointer(), U.GetArrayPointer(), &m, VT.GetArrayPointer(), &n, work.GetArrayPointer(), &lwork, &info);

  return;
}

// Print array to screen
void printarray(Array1D<double>& x){
    int n = x.Length();
    cout << endl;
    cout << "====================================================" << endl;
    cout << "n-dim double Array " << endl;
    cout << "with n = " << n << endl;
    cout << "----------------------------------------------------" << endl;
    cout << "[( ";
    int cutoff = 25;
    if (n < cutoff){
        for (int i = 0; i < x.Length(); i++){
            cout << x(i) << ", ";
        }
        cout << ")]" << endl;
    }
    else{
        for (int i = 0; i < cutoff; i++){
            cout << x(i) << ", ";
        }
        cout << " ... ";
        cout << ")]" << endl;
    }
}

// Print array to screen
void printarray(Array1D<int>& x){
  int n = x.Length();
  cout << endl;
  cout << "====================================================" << endl;
  cout << "n-dim integer Array " << endl;
  cout << "with n = " << n << endl;
  cout << "----------------------------------------------------" << endl;
  cout << "[( ";
  for (int i = 0; i < x.Length(); i++){
    cout << x(i) << ", ";
  }
  cout << ")]" << endl;
  return;
}

// Print array to screen
void printarray(Array2D<double>& x){
  int m = x.XSize();
  int n = x.YSize();
  cout << "====================================================" << endl;
  cout << "mxn double Array " << endl;
  cout << "with m = " << m << ", and n = " << n << endl;
  cout << "----------------------------------------------------" << endl;
  for(int ip=0; ip < m; ip++){
    cout << setw(5) << ip+1 << "   | ";
    for(int idim=0; idim < n; idim++){
      cout << setw(10) << x(ip,idim) << "  ";
    }
    cout << endl;
  }
  return;
}

// Print array to screen
void printarray(Array2D<int>& x){
  int m = x.XSize();
  int n = x.YSize();
    cout << "====================================================" << endl;
  cout << "mxn integer Array " << endl;
  cout << "with m = " << m << ", and n = " << n << endl;
  cout << "----------------------------------------------------" << endl;
  for(int ip=0; ip < m; ip++){
    cout << setw(5) << ip+1 << "   | ";
    for(int idim=0; idim < n; idim++){
      cout << setw(10) << x(ip,idim) << "  ";
    }
    cout << endl;
  }
  return;
}
