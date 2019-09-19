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

     Questions? Contact the UQTk Developers at <uqtk-developers@software.sandia.gov>
     Sandia National Laboratories, Livermore, CA, USA
===================================================================================== */
/** \file minmax.cpp
 * \brief Tools to find min/max values of arrays.
 */

#include "Array1D.h"
#include "Array2D.h"
#include <math.h>
#include "minmax.h"

// Get the domain bounds given a multidimensional data
void getDomain(Array2D<double>& data_in,Array1D<double>& a, Array1D<double>& b)
{
  // Get the sizes of data
  int nsample=data_in.XSize();
  int ndim=data_in.YSize();

  // Set the sizes of domain bounds
  a.Resize(ndim);
  b.Resize(ndim);
    
  // Work array
  Array1D<double> data_1d(nsample,0.e0);

  // For each dimension, get the associated data and get the bounds
  for (int idim=0;idim<ndim;idim++){
    for (int is=0;is<nsample;is++){
      data_1d(is)=data_in(is,idim);
    }
    // Add 'cushions'
    a(idim)=data_1d(minIndex(data_1d))-1.e-6;;
    b(idim)=data_1d(maxIndex(data_1d))+1.e-6;;
  }

  return;
}


double maxVal(Array1D<double>& vector)
{

  double maxVal_ = vector(0);
  for(int i=1;i< (int) vector.XSize();i++)
    if (vector(i) > maxVal_) maxVal_ = vector(i);
    
  return maxVal_;

}


int maxVal(const Array1D<int>& vector)
{

  int maxVal_ = vector(0);
  for(int i=1;i< (int) vector.XSize();i++)
    if (vector(i) > maxVal_) maxVal_ = vector(i);
    
  return maxVal_;

}

double maxVal(const Array2D<double>& vector)
{

  double maxVal_ = vector(0,0);
  for(int j=0;j< (int) vector.YSize();j++)
    for(int i=0;i< (int) vector.XSize();i++)
      if (vector(i,j) > maxVal_) maxVal_ = vector(i,j);

  return maxVal_;

}

int maxVal(const Array2D<int>& vector)
{

  int maxVal_ = vector(0,0);
  for(int j=0;j< (int) vector.YSize();j++)
    for(int i=0;i< (int) vector.XSize();i++)
      if (vector(i,j) > maxVal_) maxVal_ = vector(i,j);

  return maxVal_;

}

double minVal(const Array1D<double>& vector)
{
  double minVal_ = vector(0);
  for(int i=1;i<(int) vector.XSize();i++){
    if (vector(i) < minVal_){
      minVal_ = vector(i);
    }
  }
  return minVal_;
}

int minVal(const Array1D<int>& vector)
{
  int minVal_ = vector(0);
  for(int i=1;i<(int) vector.XSize();i++){
    if (vector(i) < minVal_){
      minVal_ = vector(i);
    }
  }
  return minVal_;
}

double minVal(const Array2D<double>& vector)
{

  double minVal_ = vector(0,0);
  for(int j=0;j< (int) vector.YSize();j++)
    for(int i=0;i< (int) vector.XSize();i++)
      if (vector(i,j) < minVal_) minVal_ = vector(i,j);

  return minVal_;

}

int minVal(const Array2D<int>& vector)
{

  int minVal_ = vector(0,0);
  for(int j=0;j< (int) vector.YSize();j++)
    for(int i=0;i< (int) vector.XSize();i++)
      if (vector(i,j) < minVal_) minVal_ = vector(i,j);

  return minVal_;

}

int maxIndex(Array1D<double>& vector)
{
  int maxInd_ = 0;
  for(int i=1;i< (int) vector.XSize();i++){
    if (vector(i) > vector(maxInd_)){
      maxInd_ = i;
    }
  }
  return maxInd_;
}

int minIndex(Array1D<double>& vector)
{
   int minInd_ = 0;
   for(int i=1;i< (int) vector.XSize();i++){
     if (vector(i) < vector(minInd_)){
       minInd_ = i;
     }
   }
   return minInd_;
}

int maxIndex(Array1D<int>& vector)
{
  int maxInd_ = 0;
  for(int i=1;i< (int) vector.XSize();i++){
    if (vector(i) > vector(maxInd_)){
      maxInd_ = i;
    }
  }
  return maxInd_;
}

 int minIndex(Array1D<int>& vector)
 {
   int minInd_ = 0;
   for(int i=1;i< (int) vector.XSize();i++){
     if (vector(i) < vector(minInd_)){
       minInd_ = i;
     }
   }
   return minInd_;
 }


// int maxIndexR_2D(const Array2D<double> a2d, const int irow)
// {
//   if ( ( irow < 0 ) ||( irow >= (int) a2d.XSize() ) ) {
//     printf("Error in maxIndexR_2D() : illegal row index %d\n",irow) ;
//     exit(1) ;
//   }
 
//   int maxInd_ = 0;
//   for( int j = 1; j < (int) a2d.YSize(); j++){
//     if (a2d(irow,j) > a2d(irow,maxInd_)){
//       maxInd_ = j;
//     }
//   }

//   return ( maxInd_ ) ;
// }


// int minIndexR_2D(const Array2D<double> a2d, const int irow)
// {
//   if ( ( irow < 0 ) ||( irow >= (int) a2d.XSize() ) ) {
//     printf("Error in minIndexR_2D() : illegal row index %d\n",irow) ;
//     exit(1) ;
//   }
 
//   int minInd_ = 0;
//   for( int j = 1; j < (int) a2d.YSize(); j++){
//     if (a2d(irow,j) < a2d(irow,minInd_)){
//       minInd_ = j;
//     }
//   }

//   return ( minInd_ ) ;
// }



int maxIndexC_2D(const Array2D<double>& a2d, const int icol)
{
  if ( ( icol < 0 ) ||( icol >= (int) a2d.YSize() ) ) {
    printf("Error in maxIndexC_2D() : illegal column index %d\n",icol) ;
    exit(1) ;
  }
 
  int maxInd_ = 0;
  for( int i = 1; i < (int) a2d.XSize(); i++){
    if (a2d(i,icol) > a2d(maxInd_,icol)){
      maxInd_ = i;
    }
  }

  return ( maxInd_ ) ;
}



int minIndexC_2D(const Array2D<double>& a2d, const int icol)
{
  if ( ( icol < 0 ) ||( icol >= (int) a2d.YSize() ) ) {
    printf("Error in minIndexC_2D() : illegal column index %d\n",icol) ;
    exit(1) ;
  }
 
  int minInd_ = 0;
  for( int i = 1; i < (int) a2d.XSize(); i++){
    if (a2d(i,icol) < a2d(minInd_,icol)){
      minInd_ = i;
    }
  }

  return ( minInd_ ) ;
}
