// File: figtree.cpp
// Created:  11-03-06 by Vlad Morariu
//
// Modified:  6-22-07 by Vlad Morariu 
//   Initial changes from previous version of the IFGT code (written by Vikas C.
//   Raykar and Changjiang Yang) and FIGTree code (written by Vikas C. Raykar).  
//
//   Modifications include:
//   1) Code can compile into a dynamic library that provides C-style interface
//      without requiring Matlab.
//
//   2) Added an improved parameter selection method that removes assumption that 
//      sources are uniformly distributed (observed large speedup in cases where 
//      sources were not uniformly distributed, and often little slowdown from 
//      overhead when the sources were actually uniformly distributed).
//
//   3) Changed the IFGT code to take multiple sets of weights for same set of 
//      sources and targets instead of having to call IFGT code multiple times.  
//      By computing a set of coefficients for each weight set, much overhead is 
//      saved (eg. computing monomials, and so on), resulting in significant 
//      speedup.
//
//   4) Added function (figtree()) that performs all parameter selection/clustering 
//      using any choice of parameter selection and evaluation algorithms.
// 
//   5) Some bugs/problem cases were fixed (some bugs caused seg faults, others 
//      were certain problem cases that could result in bad parameter selection 
//      and, as a result, memory allocation errors or seg faults).
//
//   6) In the original implementation, most code resided in the constructor and 
//      Evaluate() functions of a class, and was actually called in sequential 
//      order as if it were a C function (thus not using any real advantages of 
//      C++ classes).  Thus, all code except for that of KCenterClustering, which
//      seems to fit better in a class, has been put in C-style functions inside 
//      of figtree.cpp.  The original location of the original source is indicated 
//      in figtree.cpp before each function.
//  
//   7) Stylistic changes (eg. variable naming conventions, function names, ...)
//    
// Modified:  9-23-07 by Vlad Morariu 
//   Change code to compile on linux and solaris.
//
// Modified: 10-03-07 by Vlad Morariu 
//   Remove requirement that data is in unit hypercube by adding 
//   maxRange parameter to figtreeChoose* functions.
//
// Modified: 01-22-08 by Vlad Morariu 
//   Rename library to FIGTree (and some
//   other function remanimg)
//
// Modified: 02-20-08 by Vlad Morariu
//   Added nchoosek_double function to use 'double' instead of 'int' to prevent
//   overflow issues.  The overflow would cause incorrect parameter estimation
//   which then resulted in out of memory errors.
//
// Modified: 02-21-08 by Vlad Morariu
//   Allow rx to be zero (each pt has a cluster center on it), and allow
//   figtreeChooseParametersNonUniform to choose a value of K that gives rx=0. 
//   In some cases in higher dimensions, it is significantly cheaper to have
//   a center at each pt (i.e. rx=0) than having even one cluster with nonzero
//   radius (since the radius might require excessively high pMax).
//   Also added FIGTREE_CHECK_POS_DOUBLE macro to allow rx to be zero when 
//   checking input parameters.
//
// Modified: 05-03-08 by Vlad Morariu 
//   Add method selection code.  This set of functions uses a tree data structure
//   and k-center clustering to estimate number of source neighbors, target neighbors, 
//   and ifgt parameters, which then allows us to estimate how much it would cost
//   to evaluate using any of direct, direct+tree, ifgt, or ifgt+tree.
//
// Modified: 05-27-08 by Vlad Morariu 
//   Change figtreeChooseParameters* and figtreeChooseTruncationNumber functions
//   to return the predicted errorBound.  This can then be used to check if
//   the parameters chosen will satisfy the desired error bound (they may not 
//   since we enforce a limit on pMax (the truncation number).
//
// Modified: 05-29-08 by Vlad Morariu
//   Fixed a small computational error in the error bound (when choosing truncation
//   number).  Because the bounds are loose, the desired error was still met even
//   with the computation error.
//
// Modified: 05-29-08 to 06-10-08 by Vlad Morariu
//   A few changes were made:
//   1) Added code to choose individual truncation numbers for both targets and sources
//       using pointwise error bounds. 
//   2) Added code to choose individual truncation numbers for targets and sources
//       using clusterwise error bounds.
//   3) Reuse K-center clustering computed during method selection if 
//       FIGTREE_EVAL_AUTO is chosen.
//   4) Changed ANN code to compute unordered nearest neighbors, saving time
//       by not using a priority queue and also because now the fixed radius 
//       nearest neighbor computation and retrieval is done in one step, 
//       instead of first finding # of nn's and then doing the search again to
//       retrieve the nn's.  This really speeds up direct+tree since the ANN
//       priority queue was implemented using insertion sort.
//
// Modified: 11-02-08, 12-01-08 to 12-05-08 by Vlad Morariu
//   Made some revisions before posting new version online.  
//   1) Changed interface of figtree so users can choose truncation method
//      (also changed figtree() to automatically revert to the simplest 
//      truncation method in cases where the two more complex methods 
//      cannot give a speedup).
//   2) Revised some comments
//   3) Found all parameters used throughout code and defined constants for them
//      so that users can change them and recompile.  In future releases, 
//      these should only be defaults, and users should be able to modify them
//      at runtime.
//   4) Removed most helper functions from figtree.h (all but figtree(), 
//      figtreeChooseEvaluationMethod(), and figtreeKCenterClustering() ) and 
//      placed them in figtree_internal.h.
//   5) Changed floating op estimation functions to reflect revised versions 
//      of code
//
// Modified: 2010/05/12 by Vlad Morariu
//   Added stdlib.h and string.h includes to for exit(), strcmp(), memset(), etc.
//   (they used to be implicitly included by gcc headers, but they are not anymore).
//
//
//------------------------------------------------------------------------------
// The code was written by Vlad Morariu, Vikas Raykar, and Changjiang Yang 
// and is copyrighted under the Lesser GPL: 
//
// Copyright (C) 2008 Vlad Morariu and Vikas Raykar and Changjiang Yang 
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 or later.
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
// See the GNU Lesser General Public License for more details. 
// You should have received a copy of the GNU Lesser General Public
// License along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place - Suite 330, Boston, 
// MA 02111-1307, USA.  
//
// The author may be contacted via email at:
// morariu(at)umd(.)edu, vikas(at)umiacs(.)umd(.)edu, cyang(at)sarnoff(.)com
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// Parameters used by the algorithm
//
// Currently, the code must be recompiled when these parameters are changed.
// In future releases, users should be able to dynamically set these at runtime.
//
// FLOPS_EXP is especially important since each processor might have different
//    values for this, which will affect prediction performance
//------------------------------------------------------------------------------

#define P_UPPER_LIMIT 100        // upper limit on truncation numbers
//#define C_UPPER_LIMIT 536870912  // upper limit for total amt of memory 
                                 // coefficients can use... NOT USED FOR NOW

// these parameters are all for the method selection portion of the code
#define M_SAMPLE                  50  // number of queries when querying tree
#define N_SS_MIN                 100  // min number of source samples
#define N_SS_POW                 .75  // the exponent for determining size of 
                                      //   subsampled set Nss = N^N_SS_POW
#define K_LIMIT_TO_AVG_NBR_RATIO   2  // the max K we allow as ratio of avg source neighbors
#define FLOPS_EXP                 28  // how many floating point ops does exp() 
                                      //   take (machine dependent)

//------------------------------------------------------------------------------
// headers
//------------------------------------------------------------------------------
#define _SECURE_SCL 0

#include "figtree.h"
#include "figtree_internal.h"

#include <stddef.h>            // for definition of NULL
#include <math.h>              // for rounding (floor)
#include <stdio.h>             // for printf
#include <string.h>            // for memset()

#include <algorithm>           // for lower_bound and random_sample
#include <functional>          // for greater<double>

#include "KCenterClustering.h" // provides class for KCenterClustering

#ifndef FIGTREE_NO_ANN
#include "ANN.h"           // ANN library used for kd-tree in FIGTree code
#endif

#ifndef INT_MAX
#include <limits.h>
#endif

#ifndef DBL_MAX
#include <float.h>
#endif 

// define MAX and MIN if not yet defined
#ifndef MAX
#define MAX(a,b) ((a) > (b) ? (a) : (b) )
#endif

#ifndef MIN
#define MIN(a,b) ((a) < (b) ? (a) : (b) )
#endif

//------------------------------------------------------------------------------
// change some functions (printf, new, delete) to use Matlab functions
//   instead if we are compiling library for use in matlab mex file.
//------------------------------------------------------------------------------
#ifdef FIGTREE_USE_MATLAB_MEX
#include "mex.h"

// use mexPrintf instead of printf for messages
#undef printf
#define printf mexPrintf

/*  
//Commented this out because when ANN is compiled as part of the same DLL
//it uses the 'new()' definition below, which significantly slows it down.
//I am not sure why the code in the ANN source files uses the 'new()' 
//definition in figtree.cpp.
inline
void * operator new (size_t size)
{
  void *p=mxMalloc(size);
  return p;
}

inline
void operator delete (void *p)
{
  mxFree(p); 
}
*/
#endif

//------------------------------------------------------------------------------
// define some macros to use for checking values of input args
// without the macros, all the if statements take up a lot of space
//------------------------------------------------------------------------------
#define FIGTREE_CHECK_POS_NONZERO_DOUBLE( VAR, FCN )                 \
  if( (VAR) <= 0.0 )                                                 \
  {                                                                  \
    printf( #FCN ": Input '" #VAR "' must be a positive number.\n"); \
    return -1;                                                       \
  }

#define FIGTREE_CHECK_POS_DOUBLE( VAR, FCN )                         \
  if( (VAR) < 0.0 )                                                 \
  {                                                                  \
    printf( #FCN ": Input '" #VAR "' must be a positive number.\n"); \
    return -1;                                                       \
  }

#define FIGTREE_CHECK_POS_NONZERO_INT( VAR, FCN )                    \
  if( (VAR) <= 0.0 )                                                 \
  {                                                                  \
    printf( #FCN ": Input '" #VAR "' must be a positive number.\n"); \
    return -1;                                                       \
  }

#define FIGTREE_CHECK_NONNULL_PTR( VAR, FCN )            \
  if( (VAR) == NULL )                                    \
  {                                                      \
  printf( #FCN ": Input pointer '" #VAR "' is NULL.\n"); \
    return -1;                                           \
  }

#ifndef ANN_H
#define ANNpointArray void*
#define ANNkd_tree void*
#endif

typedef struct _FigtreeData
{
//  // general params - not used yet, but will be used when saving data-structured for multiple calls
//  int d;
//  int N;
//  int M;
//  int W;
//  double epsilon;
//  double * x;
//  double h;
//  double * q; 
//  double * y;
  
  // params for IFGT
  int pMax;
  int pMaxTotal;
  int K;
  int * clusterIndex;
  double * clusterCenters;
  double * clusterRadii;
  int * numPoints;
  double r;
  double rx;

  // params for IFGT + Tree
  ANNpointArray annClusters;
  ANNkd_tree * annClustersKdTree;
  
  // params for Direct + Tree
  ANNpointArray annSources;
  ANNkd_tree * annSourcesKdTree;

} FigtreeData;

FigtreeData figtreeCreateData()
{
  FigtreeData data;

  data.pMax = 0;
  data.pMaxTotal = 0;
  data.K = 0;
  data.clusterIndex = NULL;
  data.clusterCenters = NULL;
  data.clusterRadii = NULL;
  data.numPoints = NULL;
  data.r = 0;
  data.rx = 0;

  data.annClusters = NULL;
  data.annClustersKdTree = NULL;

  data.annSources = NULL;
  data.annSourcesKdTree = NULL;

  return data;
}

void figtreeReleaseData( FigtreeData * data )
{
  data->pMax = 0;
  data->pMaxTotal = 0;
  data->K = 0;
  if( data->clusterIndex != NULL )
  {
    delete [] data->clusterIndex;
    data->clusterIndex = NULL;
  }
  if( data->clusterCenters != NULL )
  {
    delete [] data->clusterCenters;
    data->clusterCenters = NULL;
  }
  if( data->clusterRadii != NULL )
  {
    delete [] data->clusterRadii;
    data->clusterRadii = NULL;
  }
  if( data->numPoints != NULL )
  {
    delete [] data->numPoints;
    data->numPoints = NULL;
  }
  data->r = 0;
  data->rx = 0;

#ifndef FIGTREE_NO_ANN
  if( data->annClusters != NULL )
  {
    annDeallocPts(data->annClusters);
    data->annClusters = NULL;
  }
  if( data->annClustersKdTree != NULL )
  {
    delete data->annClustersKdTree;  
    data->annClustersKdTree = NULL;
  }

  if( data->annSources != NULL )
  {
    annDeallocPts(data->annSources);
    data->annSources = NULL;
  }
  if( data->annSourcesKdTree != NULL )
  {
    delete data->annSourcesKdTree;  
    data->annSourcesKdTree = NULL;
  }
#endif
}

////////////////////////////////////////////////////////////////////////////////
// Helper functions (their prototpyes do not appear in the header file).
////////////////////////////////////////////////////////////////////////////////

//------------------------------------------------------------------------------
// Compute the combinatorial number nchoosek.
// Originally from ImprovedFastGaussTransform.cpp (IFGT source code)
//
// Modified by Vlad Morariu on 2007-06-19.
//------------------------------------------------------------------------------
int nchoosek(int n, int k)
{
  int n_k = n - k;
  
  if (k < n_k)
  {
    k = n_k;
    n_k = n - k;
  }

  int nchsk = 1; 
  for ( int i = 1; i <= n_k; i++)
  {
    nchsk *= (++k);
    nchsk /= i;
  }

  return nchsk;
}

//------------------------------------------------------------------------------
// Compute the combinatorial number nchoosek, using double precision.
// This prevents some overflow issues for large n.
//
// Created by Vlad Morariu on 2008-02-20.
//------------------------------------------------------------------------------
double nchoosek_double(int n, int k)
{
  int n_k = n - k;
  
  if (k < n_k)
  {
    k = n_k;
    n_k = n - k;
  }

  double nchsk = 1; 
  for ( int i = 1; i <= n_k; i++)
  {
    nchsk *= (++k);
    nchsk /= i;
  }

  return nchsk;
}

//------------------------------------------------------------------------------
// This function computes the constants  2^alpha/alpha!.
// Originally compute_constant_series from ImprovedFastGaussTransform.cpp (IFGT 
// source code).
//
// Modified by Vlad Morariu on 2007-06-19.
//------------------------------------------------------------------------------
void computeConstantSeries( int d, int pMaxTotal, int pMax, double * constantSeries )
{ 
  int *heads = new int[d+1];
  int *cinds = new int[pMaxTotal];
  
  for (int i = 0; i < d; i++)
    heads[i] = 0;
  heads[d] = INT_MAX;
  
  cinds[0] = 0;
  constantSeries[0] = 1.0;
  for (int k = 1, t = 1, tail = 1; k < pMax; k++, tail = t)
  {
    for (int i = 0; i < d; i++)
    {
      int head = heads[i];
      heads[i] = t;
      for ( int j = head; j < tail; j++, t++)
      {
        cinds[t] = (j < heads[i+1])? cinds[j] + 1 : 1;
        constantSeries[t] = 2.0 * constantSeries[j];
        constantSeries[t] /= (double) cinds[t];
      }
    }
  }
  
  delete [] cinds;
  delete [] heads; 
}

//------------------------------------------------------------------------------
// This function computes the monomials [(x_i-c_k)/h]^{alpha} and 
// norm([(x_i-c_k)/h])^2.
// Originally compute_source_center_monomials from 
// ImprovedFastGaussTransform.cpp (IFGT source code).
//
// Modified by Vlad Morariu on 2007-06-19.
//------------------------------------------------------------------------------
void computeSourceCenterMonomials( int d, double h, double * dx, 
                                   int p, double * sourceCenterMonomials )
{    
  int * heads = new int[d];

  for (int i = 0; i < d; i++)
  {
    dx[i]=dx[i]/h;
    heads[i] = 0;
  }
    
  sourceCenterMonomials[0] = 1.0;
  for (int k = 1, t = 1, tail = 1; k < p; k++, tail = t)
  {
    for (int i = 0; i < d; i++)
    {
      int head = heads[i];
      heads[i] = t;
      for ( int j = head; j < tail; j++, t++)
        sourceCenterMonomials[t] = dx[i] * sourceCenterMonomials[j];
    }            
  }          

  delete [] heads;
}

//------------------------------------------------------------------------------
// This function computes the monomials [(y_j-c_k)/h]^{alpha}
// Originally compute_target_center_monomials from 
// ImprovedFastGaussTransform.cpp (IFGT source code).
//
// Modified by Vlad Morariu on 2007-06-19.
//------------------------------------------------------------------------------
void computeTargetCenterMonomials( int d, double h, double * dy, 
                                   int pMax, double * targetCenterMonomials )
{    
  int *heads = new int[d];

  for (int i = 0; i < d; i++)
  {
    dy[i] = dy[i]/h;
    heads[i] = 0;
  }
    
  targetCenterMonomials[0] = 1.0;
  for (int k = 1, t = 1, tail = 1; k < pMax; k++, tail = t)
  {
    for (int i = 0; i < d; i++)
    {
      int head = heads[i];
      heads[i] = t;
      for ( int j = head; j < tail; j++, t++)
        targetCenterMonomials[t] = dy[i] * targetCenterMonomials[j];
    }            
  }          

  delete [] heads;
}
//------------------------------------------------------------------------------
// Given error(a,b,p) = (2^p/p!) * (a/h)^p * (b/h)^p * e^(-(a-b)^2/h^2), 
// calculates the maximum error given only a and the maximum possible b value
//   a - radius of source (target) point
//   b_max - maximum radius of target (source) point
//   c - constant term (2^p/p!)
//   h2 - bandwidth, squared
//   p - truncation number
//
// Created by Vlad Morariu on 2008-06-04.
//------------------------------------------------------------------------------
inline
double figtreeOneSidedErrorBound( double a, double b_max, double c, double h2, int p )
{
  double b = MIN( b_max, .5*(a + sqrt(a*a + 2*p*h2)) ); // this is the value of at which error(a,b,p) reaches maximum for a given 'a' and 'p'
  double d_ab = a-b;
  return c * pow(a*b/h2,p) * exp( -d_ab*d_ab/h2 );
}


//------------------------------------------------------------------------------
// Let error(a,b_max,p) = (2^p/p!) * (a/h)^p * (b/h)^p * e^(-(a-b)^2/h^2). This function
// finds radius values of lo_out and hi_out such that error(a_lo,b_max,p) <= epsilon < error(a_hi,b_max,p)
// We assume that the input arguments a_hi and a_lo initially satisfy this property.
//   a - radius of source (target) point
//   b_max - maximum radius of target (source) point
//   c - constant term (2^p/p!)
//   h2 - bandwidth, squared
//   p - truncation number 
//   epsilon - desired error bound
//   max_it - number of iterations of halving the interval [lo,hi]
//
// Created by Vlad Morariu on 2008-06-04
//------------------------------------------------------------------------------
void figtreeFindRadiusBounds( double a_lo, double a_hi, double b_max, 
                              double c, double h2, int p, double epsilon,
                              int max_it, double * lo_out, double * hi_out )
{
  // compute bounds at hi
  bool sat_hi = (figtreeOneSidedErrorBound( a_hi, b_max, c, h2 ,p ) <= epsilon);
  if( sat_hi )
  {
    // the bounds are already satisfied even at a_hi, where error is highest
    *hi_out = a_hi;
    *lo_out = a_hi;
  }
  else
  {
    bool sat_lo = (figtreeOneSidedErrorBound( a_lo, b_max, c, h2 ,p ) <= epsilon);
    if( !sat_lo )
    {
      // the bounds are not satisfied at a_lo (and since we assume error to increase 
      // monotonically from a_lo to a_hi, it is not satisfied in this range)
      *hi_out = a_hi;
      *lo_out = 2*a_lo - a_hi;  // go a little past a_lo to signal that not even a_lo satisfied error
    }
    else
    {
      for( int i = 0; i < max_it; i++ )
      {
        double a_mid = .5*(a_lo + a_hi);
        bool sat_mid = (figtreeOneSidedErrorBound( a_mid, b_max, c, h2 ,p ) <= epsilon);
        if( sat_mid ) // move the lo or hi value to keep property that error is not satisfied at a_hi, but is satisfied at a_lo
          a_lo = a_mid;
        else
          a_hi = a_mid;
      }
      *hi_out = a_hi;
      *lo_out = a_lo;
    }
  }
}


//------------------------------------------------------------------------------
// This function precomputes, for each truncation number, the range in the
//   distance of a source from the cluster center so that error is still satisfied.
//   This is used in the point-wise adaptive version of the IFGT. 
//
// Created by Vlad Morariu on 2008-06-04.
//------------------------------------------------------------------------------
void figtreeSourceTruncationRanges( double r, double rx, double h, double epsilon, int pMax, double * max_source_dists2 )
{
  double h2 = h*h;
  for( int i = 0; i < pMax-1; i++ )
    max_source_dists2[i] = -1; // negative numbers indicate no distance satisfies error bounds for the particular value of p
  max_source_dists2[pMax-1] = rx;

  double c = 1;
  for( int i = 0; i < pMax-1; i++ )
  {
    c *= (2.0/(i+1));
    double a_lo = 0, a_hi = rx;
    figtreeFindRadiusBounds( a_lo, a_hi, r + rx, c, h2, i+1, epsilon, 10, &a_lo, &a_hi );
    max_source_dists2[i] = a_lo*a_lo;
  }
}

//------------------------------------------------------------------------------
// This function precomputes, for each truncation number, the range in the
//   distance of a target from the cluster center so that error is still satisfied.
//   This is used in the point-wise adaptive version of the IFGT.
//
// Created by Vlad Morariu on 2008-06-04.
//------------------------------------------------------------------------------
void figtreeTargetTruncationRanges( double r, double rx, double h, double epsilon, int pMax, double * max_target_dists2, double * min_target_dists2 )
{
  double h2 = h*h;
  double ry = r + rx;
  for( int i = 0; i < pMax-1; i++ )
  {
    max_target_dists2[i] = -1;   
    min_target_dists2[i] = ry*ry+1;
  }

  double c = 1;
  for( int i = 0; i < pMax-1; i++ )
  {
    c *= (2.0/(i+1));

    double peak_dist = .5*(rx + sqrt( rx*rx + 2*h2*(i+1) ));

    // here we calculate for each value of p the maximum distance from a cluster center 
    //   that a target can be to satisfy the error bounds, provided the distance is 
    //   in the portion of the error bound that monotonically increases with distance
    double a_lo = 0, a_hi = MIN(ry,peak_dist);
    figtreeFindRadiusBounds( a_lo, a_hi, rx, c, h2, i+1, epsilon, 10, &a_lo, &a_hi );
    max_target_dists2[i] = a_lo*a_lo;

    // here we calculate for each value of p the minimum distance from a cluster center 
    //   that a target can be to satisfy the error bounds, provided the radius is 
    //   in the portion of the error bound that monotonically increases with distance
    if( peak_dist <= ry )
    {
      a_lo = ry, a_hi = peak_dist;
      figtreeFindRadiusBounds( a_lo, a_hi, rx, c, h2, i+1, epsilon, 10, &a_lo, &a_hi );
      min_target_dists2[i] = a_lo*a_lo;
    } // otherwise we leave it at ry*ry+1 (it should have been initialized to this

    if( i > 0 && min_target_dists2[i] > min_target_dists2[i-1] )
    { 
      min_target_dists2[i] = min_target_dists2[i-1];
    }
  }
  if( pMax > 1 && min_target_dists2[pMax-1] > min_target_dists2[pMax-2] )
  { 
    min_target_dists2[pMax-1] = min_target_dists2[pMax-2];
  }

}


//------------------------------------------------------------------------------
// Given the precomputed distances from the cluster center for which 
//   the error bound is satisfied for each truncation number, and the actual
//   distance from the cluster center, this function finds the lowest truncation
//   number so that error is still satisfied.
//
// Created by Vlad Morariu on 2008-06-04.
//------------------------------------------------------------------------------
inline
int figtreeSourceTruncationNumber( double dx2, int pMax, double * max_source_dists2 )
{
  return (int)(std::lower_bound( max_source_dists2, max_source_dists2 + pMax - 1, dx2) - max_source_dists2) + 1;  
}

//------------------------------------------------------------------------------
// Given the precomputed distances from the cluster center for which 
//   the error bound is satisfied for each truncation number, and the actual
//   distance from the cluster center, this function finds the lowest truncation
//   number so that error is still satisfied.
//
// Created by Vlad Morariu on 2008-06-04.
//------------------------------------------------------------------------------
inline
int figtreeTargetTruncationNumber( double dy2, int pMax, double * max_target_dists2, double * min_target_dists2 )
{
  if( dy2 <= max_target_dists2[pMax-2] )
    return (int)(std::lower_bound( max_target_dists2, max_target_dists2 + pMax - 1, dy2) - max_target_dists2) + 1;
  else if( dy2 >= min_target_dists2[pMax-2] )
    return (int)(std::lower_bound( min_target_dists2, min_target_dists2 + pMax - 1, dy2, std::greater<double>() ) - min_target_dists2) + 1;
  else 
    return pMax;
}

//------------------------------------------------------------------------------
// This function computes the coefficients C_k for all clusters.
// Originally compute_C from ImprovedFastGaussTransform.cpp (IFGT source code)
//
// Modified by Vlad Morariu on 2007-06-19.
//------------------------------------------------------------------------------
void computeC( int d, int N, int W, int K, int pMaxTotal, int pMax, 
               double h, int * clusterIndex, double * x, double * q,
               double * clusterCenter, double * C )
{
  double * sourceCenterMonomials = new double[pMaxTotal];
  double * constantSeries = new double[pMaxTotal];
  double hSquare = h*h;
  double * dx = new double[d];

  for (int i = 0; i < W*K*pMaxTotal; i++)
  {
    C[i] = 0.0;
  }

  for(int i = 0; i < N; i++)
  {
    int k = clusterIndex[i];
    int sourceBase = i*d;
    int centerBase = k*d;
    double sourceCenterDistanceSquare = 0.0;

    for (int j = 0; j < d; j++)
    {
      dx[j] = (x[sourceBase+j] - clusterCenter[centerBase+j]);
      sourceCenterDistanceSquare += (dx[j]*dx[j]);
    }
  
    computeSourceCenterMonomials( d, h, dx, pMax, sourceCenterMonomials );    
    
    for(int w = 0; w < W; w++ )
    {
      double f = q[N*w + i]*exp(-sourceCenterDistanceSquare/hSquare);
      for(int alpha = 0; alpha < pMaxTotal; alpha++)
      {
          C[(K*w + k)*pMaxTotal + alpha] += (f*sourceCenterMonomials[alpha]);
      }
    }  
  }

  computeConstantSeries( d, pMaxTotal, pMax, constantSeries );

  for(int w = 0; w < W; w++)
  {   
    for(int k = 0; k < K; k++)
    {
      for(int alpha = 0; alpha < pMaxTotal; alpha++)
      {
        C[(K*w + k)*pMaxTotal + alpha] *= constantSeries[alpha];
      } 
    }
  }

  delete [] sourceCenterMonomials;
  delete [] constantSeries;
  delete [] dx;
}

//------------------------------------------------------------------------------
// This function computes a separate truncation number for each cluster that 
//   satisfies the total allowed cluster-wise error.  This means that some 
//   sources are allowed to contribute more error as long as others contribute 
//   less.  The speedup is observed mostly in higher dimensions since finding 
//   the cluster-wise truncation numbers adds to the overhead cost.  The 
//   truncations for each cluster are found by doing a binary search over 
//   the truncation number, p.
// 
// See 'Automatic online tuning for fast Gaussian summation,' by Morariu et al, 
//   NIPS 2008 for details.  
//
// NOTES: 1) This function currently only works when only one set of weights 'q' 
//           are used (i.e. W=1).  It can be extended to W>1, but I have not had 
//           time to clean up that part of the code for release yet.
//   
//        2) Also, this function does not currently compute multiple 'regions' 
//           for targets where different cluster-wise truncation numbers are used 
//           depending on what region targets are in.
//
// Created by Vlad Morariu on 2008-06-04.
//------------------------------------------------------------------------------
void figtreeFindClusterTruncations( int d, int N, double * x, double * q, double h, double epsilon, double r, int pMax, int K, int * clusterIndex, int * numPoints, double * clusterCenters, double * clusterRadii, int * clusterTruncations )
{
  double * clusterWeights = new double[K];
  double * pointClusterDists = new double[N];
  double h2 = h*h;

  double * q_reordered = new double[N];

  memset(clusterTruncations,0,sizeof(int)*K);
  memset(clusterWeights,0,sizeof(double)*K);
  memset(pointClusterDists,0,sizeof(double)*N);

  // find out total sum of weights for each cluster
  for( int i = 0; i < N; i++ )
    clusterWeights[ clusterIndex[i] ] += fabs(q[i]);

  // precompute a list of points for each cluster so that we can reorganize
  //   some frequently used data in memory to make memory is accessed
  //   sequentially.
  int * clusterMembers = new int[N];
  int * clusterStart = new int[K], * clusterEnd = new int[K];
  clusterStart[0] = 0;
  clusterEnd[0] = clusterStart[0];
  for( int i = 1; i < K; i++ )
  {
    clusterStart[i] = clusterStart[i-1] + numPoints[i-1];
    clusterEnd[i] = clusterStart[i];
  }
 
  for( int i = 0; i < N; i++ )
    clusterMembers[clusterEnd[clusterIndex[i]]++] = i;

  // compute the distance from points to clusters, and 
  //   reorder them so that the weights 'q_reordered' and cluster 
  //   distances 'pointClusterDists' for points that belong to the
  //   same cluster are contiguous in memory (reduces memory access overhead).
  for( int k = 0; k < K; k++ )
  {
    int start = clusterStart[k], end = clusterEnd[k];
    for( int i = start; i < end; i++ )
    {
      q_reordered[i] = fabs(q[clusterMembers[i]]) / clusterWeights[k];
      for( int j = 0; j < d; j++ )
      {
        double dx = clusterCenters[k*d+j] - x[clusterMembers[i]*d+j];
        pointClusterDists[i] += dx*dx;
      }
      pointClusterDists[i] = sqrt(pointClusterDists[i]);
    }
  }

  // precompute constant in front of error term that depends only on p
  double * constants = new double[pMax];
  constants[0] = 2;
  for( int p = 2; p <= pMax; p++ )
    constants[p-1] = constants[p-2]*2.0/p;

  // find the lowest p for each cluster such that the total error per cluster meets the error bound using
  // a binary search based algorithm.  If for some reason, there exists a p_i > p_j such that 
  // error is met for p_j but not for p_i, (i.e. the p's are not ordered), the resulting p is still guaranteed
  // to satisfy error, though it may not be the lowest p that satisfies error.
  for( int k = 0; k < K; k++ )
  {
    int start = clusterStart[k], end = clusterEnd[k];
  
    int p_lo = 1, p_hi = pMax;
    while( p_lo < p_hi )
    {
      int p_mid = (p_hi + p_lo)/2;
      double clusterBound = 0;
      for( int i = start; i < end && clusterBound <= epsilon; i++ )
      {
        double error = q_reordered[i]*figtreeOneSidedErrorBound( pointClusterDists[i], clusterRadii[k]+r, constants[p_mid-1], h2 , p_mid );
        clusterBound += error;
      }

      if( clusterBound > epsilon )
        p_lo = p_mid + 1;
      else
        p_hi = p_mid;
    }
    
    clusterTruncations[k] = p_hi;
  }

  delete [] constants;
  delete [] clusterWeights;
  delete [] pointClusterDists;
  delete [] clusterMembers;
  delete [] clusterStart;
  delete [] clusterEnd;
  delete [] q_reordered;
}

//------------------------------------------------------------------------------
// This function computes the cluster coefficients using a separate truncation
//   number for each cluster that satisfies the total allowed cluster-wise error.  
//   The clusterwise truncation numbers are returned by the function
//   figtreeFindClusterTruncations().  
// 
// See 'Automatic online tuning for fast Gaussian summation,' by Morariu et al, 
//   NIPS 2008 for details.  
//
// NOTES: 1) Because figtreeFindClusterTruncations() only takes one set of weights
//           into account when computing cluster-wise truncations, this function
//           will not give accurate results for multiple sets of weights (W>1).
//           When figtreeFindClusterTruncations is modified to consider all W sets
//           of weights, this function should work properly.  Thus, W=1 is assumed
//           for now even though it is part of the input arguments.
//   
//        2) Also, this function does not currently use multiple 'regions' 
//           for targets where different cluster-wise truncation numbers are used 
//           depending on what region targets are in.
//
// Created by Vlad Morariu on 2008-06-04.
//------------------------------------------------------------------------------
void computeCAdaptiveCluster( int d, int N, int W, int K, int pMaxTotal, int pMax, 
                       double h, int * clusterIndex, double * x, double * q,
                       double * clusterCenter, int * clusterTruncations, int * pMaxTotals, double * C )
{
  double * sourceCenterMonomials = new double[pMaxTotal];
  double * constantSeries = new double[pMaxTotal];
  double hSquare = h*h;
  double * dx = new double[d];

  memset( C, 0, sizeof(double)*W*K*pMaxTotal );

  for(int i = 0; i < N; i++)
  {
    int k = clusterIndex[i];
    int sourceBase = i*d;
    int centerBase = k*d;
    double sourceCenterDistanceSquare = 0.0;

    for (int j = 0; j < d; j++)
    {
      dx[j] = (x[sourceBase+j] - clusterCenter[centerBase+j]);
      sourceCenterDistanceSquare += (dx[j]*dx[j]);
    }
  
    int p = clusterTruncations[k];
    int pTotal = pMaxTotals[p-1];
    computeSourceCenterMonomials( d, h, dx, p, sourceCenterMonomials );    
    
    for(int w = 0; w < W; w++ )
    {
      double f = q[N*w + i]*exp(-sourceCenterDistanceSquare/hSquare);
      for(int alpha = 0; alpha < pTotal; alpha++)
      {
          C[(K*w + k)*pMaxTotal + alpha] += (f*sourceCenterMonomials[alpha]);
      }
    }  
  }

  computeConstantSeries( d, pMaxTotal, pMax, constantSeries );

  for(int w = 0; w < W; w++)
  {   
    for(int k = 0; k < K; k++)
    {
      for(int alpha = 0; alpha < pMaxTotal; alpha++)
      {
        C[(K*w + k)*pMaxTotal + alpha] *= constantSeries[alpha];
      } 
    }
  }

  delete [] sourceCenterMonomials;
  delete [] constantSeries;
  delete [] dx;
}

//------------------------------------------------------------------------------
// This function computes the cluster coefficients allowing each source point 
//   to have a variable truncation number that still satisfies the desired error
//   assuming the worst-case target point.
// 
// See 'Automatic online tuning for fast Gaussian summation,' by Morariu et al, 
//   NIPS 2008 for details.  
//
// NOTES: Unlike the cluster-wise adaptive version, this does work for W>1.
//
// Created by Vlad Morariu on 2008-06-04.
//------------------------------------------------------------------------------
void computeCAdaptivePoint( int d, int N, int W, int K, int pMaxTotal, int pMax, 
                       double h, int * clusterIndex, double * x, double * q,
                       double * clusterCenter, double * maxSourceDists2, int * pMaxTotals, double * C )
{
  double * sourceCenterMonomials = new double[pMaxTotal];
  double * constantSeries = new double[pMaxTotal];
  double hSquare = h*h;
  double * dx = new double[d];

  memset( C, 0, sizeof(double)*W*K*pMaxTotal );

  //int * pHistogram = new int[pMax];
  //memset(pHistogram,0,sizeof(int)*pMax);

  for(int i = 0; i < N; i++)
  {
    int k = clusterIndex[i];
    int sourceBase = i*d;
    int centerBase = k*d;
    double sourceCenterDistanceSquare = 0.0;

    for (int j = 0; j < d; j++)
    {
      dx[j] = (x[sourceBase+j] - clusterCenter[centerBase+j]);
      sourceCenterDistanceSquare += (dx[j]*dx[j]);
    }
  
    int p = figtreeSourceTruncationNumber( sourceCenterDistanceSquare, pMax, maxSourceDists2 );
    int pTotal = pMaxTotals[p-1];
    //pHistogram[p-1]++;    
    computeSourceCenterMonomials( d, h, dx, p, sourceCenterMonomials );    
    
    for(int w = 0; w < W; w++ )
    {
      double f = q[N*w + i]*exp(-sourceCenterDistanceSquare/hSquare);
      for(int alpha = 0; alpha < pTotal; alpha++)
      {
          C[(K*w + k)*pMaxTotal + alpha] += (f*sourceCenterMonomials[alpha]);
      }
    }  
  }

  computeConstantSeries( d, pMaxTotal, pMax, constantSeries );

  for(int w = 0; w < W; w++)
  {   
    for(int k = 0; k < K; k++)
    {
      for(int alpha = 0; alpha < pMaxTotal; alpha++)
      {
        C[(K*w + k)*pMaxTotal + alpha] *= constantSeries[alpha];
      } 
    }
  }

  //printf( "source p histogram: ");
  //for( int i = 0; i < pMax; i++ )
  //  printf( " %i", pHistogram[i] );
  //printf( "\n" );
  //delete [] pHistogram;

  delete [] sourceCenterMonomials;
  delete [] constantSeries;
  delete [] dx;
}

//------------------------------------------------------------------------------
// This function provides a simple interface for evaluation of gauss transforms
// in several ways (direct, direct with approximate nearest-neighbors
// structure on sources, ifgt, and ifgt with approximate nearest-neighbors) and
// using different parameter selection methods (can assume uniform distribution
// or use the actual distribution in estimating parameters).
//
// Created by Vlad Morariu on 2006-11-03.
//
// Modified by Vlad Morariu on 2007-06-21.
//
// Modified by Vlad Morariu on 2008-06-04. 
//
//    - Added changes described in 'Automatic online tuning for fast Gaussian
//      summation,' by Morariu et al, NIPS 2008 for details).
//    - added FIGTREE_EVAL_AUTO eval method which allows FIGTREE to pick 
//      the evaluation method that is predicted to be the fastest, making
//      FIGTREE a black box apprach
//    - added point-wise and cluster-wise adaptive versions of the IFGT (instead
//      of using the max truncation number for all sources, clusters, and targets,
//      the truncation number is varied either point-wise or cluster-wise to improve
//      performance while still satisfying the desired error bound).
//
// Modified by Vlad Morariu on 2008-12-01
//    - add a third parameter ifgtTruncMethod to reduce the number of evalMethod
//      possibilities, separating evaluation method from truncation number choices
//    - removed 'forceK' option, since it was used only for debugging
//------------------------------------------------------------------------------
int figtree( int d, int N, int M, int W, double * x, double h, 
             double * q, double * y, double epsilon, double * g,
             int evalMethod, int ifgtParamMethod, int ifgtTruncMethod, int verbose )
{
  int ret = 0;

  FigtreeData data = figtreeCreateData();

  // if the evalMethod is FIGTREE_EVAL_AUTO, choose the method that is estimated to be the fastest
  if( evalMethod == FIGTREE_EVAL_AUTO )
  {
    ret = figtreeChooseEvaluationMethod( d, N, M, W, x, h, y, epsilon, ifgtParamMethod, verbose, &evalMethod, NULL, &data );
  }

  // for FIGTREE_EVAL_DIRECT and FIGTREE_EVAL_DIRECT_TREE, we don't need to compute
  //   parameters, so we just run the fcns directly, once for each set of weights
  if( evalMethod == FIGTREE_EVAL_DIRECT )
  {
    verbose && printf("figtreeEvalMethod() chose the direct method.\n");
    for( int i = 0; i < W; i++ )
      ret = figtreeEvaluateDirect( d, N, M, x, h, q+i*N, y, g+i*M );
  }

  if( evalMethod == FIGTREE_EVAL_DIRECT_TREE )
  {
    verbose && printf("figtreeEvalMethod() chose direct+tree method.\n");
    for( int i = 0; i < W; i++ )
      ret = figtreeEvaluateDirectTree( d, N, M, x, h, q+i*N, y, epsilon, g+i*M );
  }

  // for FIGTREE_EVAL_IFGT and FIGTREE_EVAL_IFGT_TREE, we must first compute
  //   parameters
  if( evalMethod == FIGTREE_EVAL_IFGT || 
      evalMethod == FIGTREE_EVAL_IFGT_TREE )
  {
    if(verbose && evalMethod == FIGTREE_EVAL_IFGT)
      verbose && printf("figtreeEvalMethod() chose the IFGT method.\n");
    if(verbose && evalMethod == FIGTREE_EVAL_IFGT_TREE)
      verbose && printf("figtreeEvalMethod() chose the IFGT+tree method.\n");

    bool alreadyHaveClustering = (data.clusterCenters != NULL); // quick and dirty test    
    double maxRange = 0;
    if( !alreadyHaveClustering )
    {
      int kLimit = N, kMax;

      //
      // calculate R, the diameter of the hypercube that encompasses sources and targets
      //
      double * mins = new double[d];
      double * maxs = new double[d];
      figtreeCalcMinMax( d, N, x, mins, maxs, 0 );
      figtreeCalcMinMax( d, M, y, mins, maxs, 1 );
      figtreeCalcMaxRange( d, mins, maxs, &maxRange );
      delete [] mins;
      delete [] maxs;

      //
      // choose parameters for IFGT
      //
      if( ifgtParamMethod == FIGTREE_PARAM_NON_UNIFORM )
        ret = figtreeChooseParametersNonUniform( d, N, x, h, epsilon, kLimit, maxRange, &kMax, &data.pMax, &data.r, NULL );
      else
        ret = figtreeChooseParametersUniform( d, h, epsilon, kLimit, maxRange, &kMax, &data.pMax, &data.r, NULL );
      if( ret < 0 )
      {
        printf("figtree: figtreeChooseParameters%sUniform() failed.\n", 
               ((ifgtParamMethod == FIGTREE_PARAM_NON_UNIFORM) ? "Non" : ""));
        return ret;
      }

      verbose && printf("figtreeChooseParameters%sUniform() chose p=%i, k=%i.\n", 
                        ((ifgtParamMethod == FIGTREE_PARAM_NON_UNIFORM) ? "Non" : ""), data.pMax, kMax );

      //
      // do k-center clustering
      //
      data.clusterIndex = new int[N];
      data.numPoints    = new int[kMax]; 
      data.clusterCenters = new double[d*kMax];
      data.clusterRadii = new double[kMax];

      ret = figtreeKCenterClustering( d, N, x, kMax, &data.K, &data.rx, data.clusterIndex, 
                                   data.clusterCenters, data.numPoints, data.clusterRadii );
      if( ret < 0 )
        printf("figtree: figtreeKCenterClustering() failed.\n");
    }

    double errorBound = epsilon + 1;
    if( ret >= 0 && !alreadyHaveClustering )
    {
      // choose truncation number again now that clustering is done
      ret = figtreeChooseTruncationNumber( d, h, epsilon, data.rx, maxRange, &data.pMax, &errorBound );
      if( ret < 0 )
        printf("figtreeChooseTruncatoinNumber() failed.\n");
      else
      {
        if( verbose && errorBound > epsilon )
          printf("figtreeChooseTruncationNumber(): could not find p within limit that satisfies error bound!\n" );
      }
    } 

    if( ret >= 0 )
    {
      // evaluate IFGT
      verbose && printf( "Eval IFGT(h= %3.2e, pMax= %i, K= %i, r= %3.2e, rx= %3.2e, eps= %3.2e)\n", 
                         h, data.pMax, data.K, data.r, data.rx, epsilon);

      // if maximum truncation is 1, then nothing can be gained by doing individual truncations
      if( data.pMax == 1 && (ifgtTruncMethod != FIGTREE_TRUNC_MAX) )
      {
        evalMethod = FIGTREE_EVAL_IFGT;
        verbose && printf("figtree(): max truncation is 1, so adaptive truncations are unnecessary.\n  Switching to FIGTREE_TRUNC_MAX...\n");
      }

      // if W>1, and we want to use adaptive truncations, we cannot use cluster-wise truncations yet
      if( W > 1 && (ifgtTruncMethod == FIGTREE_TRUNC_CLUSTER ) )
      {
        ifgtTruncMethod = FIGTREE_TRUNC_POINT;
        verbose && printf( "figtree(): FIGTREE_TRUNC_CLUSTER is not yet implemented\n  to handle W > 1.  Switching to FIGTREE_TRUNC_POINT...\n");
      }

      if( evalMethod == FIGTREE_EVAL_IFGT && ifgtTruncMethod == FIGTREE_TRUNC_POINT )
      {
        ret = figtreeEvaluateIfgtAdaptivePoint( d, N, M, W, x, h, q, y, data.pMax, data.K, data.clusterIndex, 
                                                    data.clusterCenters, data.clusterRadii, data.r, epsilon, g );
      }

      if( evalMethod == FIGTREE_EVAL_IFGT_TREE && ifgtTruncMethod == FIGTREE_TRUNC_POINT )
      {
        ret = figtreeEvaluateIfgtTreeAdaptivePoint( d, N, M, W, x, h, q, y, data.pMax, data.K, data.clusterIndex, 
                                                    data.clusterCenters, data.clusterRadii, data.r, epsilon, g );
      }

      if( ifgtTruncMethod == FIGTREE_TRUNC_CLUSTER )
      {
        int * clusterTruncations = new int[data.K];
        figtreeFindClusterTruncations( d, N, x, q, h, epsilon, data.r, data.pMax, data.K, data.clusterIndex, data.numPoints, data.clusterCenters, data.clusterRadii, clusterTruncations );
        int pMaxNew = 0;
        for( int i = 0; i < data.K; i++ )
          pMaxNew = MAX(pMaxNew, clusterTruncations[i]);
        if( evalMethod == FIGTREE_EVAL_IFGT )
        {
          ret = figtreeEvaluateIfgtAdaptiveCluster( d, N, M, W, x, h, q, y, pMaxNew, data.K, data.clusterIndex, 
                                                      data.clusterCenters, data.clusterRadii, data.r, epsilon, clusterTruncations, g );
        }
        else // evalMethod == FIGTREE_EVAL_IFGT_TREE
        {
          ret = figtreeEvaluateIfgtTreeAdaptiveCluster( d, N, M, W, x, h, q, y, pMaxNew, data.K, data.clusterIndex, 
                                                      data.clusterCenters, data.clusterRadii, data.r, epsilon, clusterTruncations, g );
        }

        delete [] clusterTruncations;
      }

      if( evalMethod == FIGTREE_EVAL_IFGT && ifgtTruncMethod == FIGTREE_TRUNC_MAX )
      {
        ret = figtreeEvaluateIfgt( d, N, M, W, x, h, q, y, data.pMax, data.K, data.clusterIndex, 
                                     data.clusterCenters, data.clusterRadii, data.r, epsilon, g );
      }

      if( evalMethod == FIGTREE_EVAL_IFGT_TREE && ifgtTruncMethod == FIGTREE_TRUNC_MAX )
      {
        ret = figtreeEvaluateIfgtTree( d, N, M, W, x, h, q, y, data.pMax, data.K, data.clusterIndex, 
                                        data.clusterCenters, data.clusterRadii, data.r, epsilon, g );
      }

      if( ret < 0 )
      {
        printf("figtree: figtreeEvaluateIfgt%s*() failed.\n", 
               ((evalMethod == FIGTREE_EVAL_IFGT_TREE) ? "Tree" : ""));
      }
    }
  } // if we are doing IFGT

  // release data if any was allocated
  figtreeReleaseData( &data );

  return ret;
}

//------------------------------------------------------------------------------
// Chooses minimum truncation number that satisfies desired error, given
//   the maximum radius of any cluster (rx).
//
// Originally constructor from 
// ImprovedFastGaussTransformChooseTruncationNumber.cpp (IFGT source code) by 
// Vikas C. Raykar.
//
// Modified by Vlad Morariu on 2007-06-20
// Modified by Vlad Morariu on 2007-10-03 - pass R (the max distance between
//     any source and target) as argument instead of assuming that data fits
//     in unit hypercube
//------------------------------------------------------------------------------
int figtreeChooseTruncationNumber( int d, double h, double epsilon, 
                                double rx, double maxRange, int * pMax, double * errorBound )
{
  // check input arguments
  FIGTREE_CHECK_POS_NONZERO_INT( d, figtreeChooseTruncationNumber );
  FIGTREE_CHECK_POS_NONZERO_DOUBLE( h, figtreeChooseTruncationNumber );
  FIGTREE_CHECK_POS_NONZERO_DOUBLE( epsilon, figtreeChooseTruncationNumber );
  FIGTREE_CHECK_POS_DOUBLE( rx, figtreeChooseTruncationNumber );
  FIGTREE_CHECK_POS_NONZERO_DOUBLE( maxRange, figtreeChooseTruncationNumber );
  FIGTREE_CHECK_NONNULL_PTR( pMax, figtreeChooseTruncationNumber );

  double R = maxRange*sqrt((double)d);
  double hSquare = h*h;
  double r = MIN(R, h*sqrt(log(1/epsilon)));
  double rxSquare = rx*rx;
  
  double error = epsilon + 1;
  double temp = 1;
  int p = 0;
  /*
  while((error > epsilon) & (p <= P_UPPER_LIMIT)){
    p++;
    double b = MIN(((rx + sqrt((rxSquare) + (2*p*hSquare)))/2), rx + r);
    double c = rx - b;
    temp = temp*(((2*rx*b)/hSquare)/p);
    error = temp*(exp(-(c*c)/hSquare));      
  } */ 
  while((error > epsilon) & (p <= P_UPPER_LIMIT))
  {
    p++;
    double b = MIN(((rx + sqrt((rxSquare) + (2*p*hSquare)))/2), rx + r);
    double c = rx - b;
    temp = 1;
    for( int i = 1; i <= p; i++ )
      temp = temp*((2.0*rx*b/hSquare)/i);
    error = temp*(exp(-(c*c)/hSquare));      
  }  
  if( pMax != NULL )
    *pMax = p;
  if( errorBound != NULL )
    *errorBound = error;

  return 0;
}

//------------------------------------------------------------------------------
// Parameter selection for the Improved Fast Gauss Transform (IFGT).
//
// Implementation based on:
//
// Fast computation of sums of Gaussians in high dimensions. 
// Vikas C. Raykar, C. Yang, R. Duraiswami, and N. Gumerov,
// CS-TR-4767, Department of computer science,
// University of Maryland, Collegepark.
//
// Originally constructor from ImprovedFastGaussTransformChooseParameters.cpp 
//   by Vikas C. Raykar. (IFGT source code)
//
// Modified by Vlad Morariu on 2007-06-20
// Modified by Vlad Morariu on 2007-10-03 - pass R (the max distance between
//     any source and target) as argument instead of assuming that data fits
//     in unit hypercube
// Modified by Vlad Morariu on 2008-06-04 - change the way bound is computed to
//     be more precise
//------------------------------------------------------------------------------
int figtreeChooseParametersUniform( int d, double h, double epsilon, 
                                 int kLimit, double maxRange, int * K, int * pMax, double * r, double * errorBound )
{
  // check input arguments
  FIGTREE_CHECK_POS_NONZERO_INT( d, figtreeChooseParametersUniform );
  FIGTREE_CHECK_POS_NONZERO_DOUBLE( h, figtreeChooseParametersUniform );
  FIGTREE_CHECK_POS_NONZERO_DOUBLE( maxRange, figtreeChooseParametersUniform );
  FIGTREE_CHECK_POS_NONZERO_DOUBLE( epsilon, figtreeChooseParametersUniform );
  FIGTREE_CHECK_POS_NONZERO_INT( kLimit, figtreeChooseParametersUniform );

  double R = maxRange*sqrt((double)d);
  double hSquare = h*h;
  double complexityMin = DBL_MAX;

  // These variables will hold the values that will then be returned.
  // We use temporary variables in case caller does not care about a variable 
  // and passes a NULL pointer.
  int kTemp = 1;
  int pMaxTemp = P_UPPER_LIMIT + 1;
  double rTemp = MIN(R,h*sqrt(log(1/epsilon)));
  double errorTemp = epsilon + 1;

  for(int i = 0; i < kLimit; i++)
  {
    double rx = maxRange*pow((double)i + 1,-1.0/(double)d);
    double rxSquare = rx*rx;
    double n = MIN(i + 1, pow(rTemp/rx,(double)d));
    double error = epsilon + 1;
    double temp = 1;
    int p = 0;
    while((error > epsilon) & (p <= P_UPPER_LIMIT))
    {
      p++;
      double b = MIN(((rx + sqrt((rxSquare) + (2*p*hSquare)))/2), rx + rTemp);
      double c = rx - b;
      temp = 1;
      for( int j = 1; j <= p; j++ )
        temp = temp*((2.0*rx*b/hSquare)/j);
      error = temp*(exp(-(c*c)/hSquare));      
    }  
    double complexity = (i + 1) + log((double)i + 1) + ((1 + n)*nchoosek_double(p - 1 + d, d));
    if (complexity < complexityMin )
    {
      complexityMin = complexity;
      kTemp = i + 1;
      pMaxTemp = p;
      errorTemp = error;
    }
  }

  // added this to catch case where desired error is never reached.
  // The best thing is to have as many clusters and terms in the taylor
  // series as possible (which will give lowest error)
  if( errorTemp > epsilon )
  {
    kTemp = kLimit;
  }

  // set output variables to computed values
  if( K != NULL )
    *K = kTemp;
  if( pMax != NULL )
    *pMax = pMaxTemp;
  if( r != NULL )
    *r = rTemp;
  if( errorBound != NULL )
    *errorBound = errorTemp;

  return 0;
}

//------------------------------------------------------------------------------
// Parameter selection scheme that does not assume uniform distribution.
// In cases where sources are not uniformly distribution, this can lead to 
// very large performance increases
// because as the number of clusters increases, the max radius of any cluster
// decreases MUCH faster than it would if the sources were uniformly distributed.
// This function is based on ImprovedFastGaussTransformChooseParameters.cpp from
// the IFGT source code, by Vikas C. Raykar.  
//
// Initially created by Vlad Morariu on 2007-01-24.
// Modified by Vlad Morariu on 2007-06-20.
// Modified by Vlad Morariu on 2007-10-03 - pass R (the max distance between
//     any source and target) as argument instead of assuming that data fits
//     in unit hypercube
// Modified by Vlad Morariu on 2008-02-21 - allow rx to be zero.  In some cases,
//     if there isn't a center at each source pt to give rx=0, an excessively
//     large pMax is needed, and it is faster to just have a center at each pt.
// Modified by Vlad Morariu on 2008-06-04 - change the way bound is computed to
//     be more precise     
// Modified by Vlad Morariu on 2008-12-05 - began changing function to incorporate
//     memory limit for coefficient storage... not done yet
//------------------------------------------------------------------------------
int figtreeChooseParametersNonUniform( int d, int N, double * x, 
                                    double h, double epsilon, int kLimit, double maxRange,
                                    int * K, int * pMax, double * r, double * errorBound )
{
  // check input arguments
  FIGTREE_CHECK_POS_NONZERO_INT( d, figtreeChooseParametersNonUniform );
  FIGTREE_CHECK_POS_NONZERO_INT( N, figtreeChooseParametersNonUniform );
  FIGTREE_CHECK_NONNULL_PTR( x, figtreeChooseParametersNonUniform );
  FIGTREE_CHECK_POS_NONZERO_DOUBLE( h, figtreeChooseParametersNonUniform );
  FIGTREE_CHECK_POS_NONZERO_DOUBLE( epsilon, figtreeChooseParametersNonUniform );
  FIGTREE_CHECK_POS_NONZERO_INT( kLimit, figtreeChooseParametersNonUniform );
  FIGTREE_CHECK_POS_NONZERO_DOUBLE( maxRange, figtreeChooseParametersNonUniform );

  // allocate temporary memory, and set some variables
  int * pClusterIndex = new int[N];
  KCenterClustering * kcc = new KCenterClustering( d, N, x, pClusterIndex, kLimit );

  double R = maxRange*sqrt((double)d);
  double hSquare = h*h;
  double complexityMin = DBL_MAX;
  double complexityLast = DBL_MAX;

  double rTemp = MIN(R,h*sqrt(log(1/epsilon)));
  int kTemp = 1;
  int pMaxTemp = P_UPPER_LIMIT + 1;
  double errorTemp = epsilon + 1;

  int numClusters;
  double rx;
  
  // Vlad 01/24/07 - add first cluster and get rx
  kcc->ClusterIncrement( &numClusters, &rx ); 

  // evaluate complexity for increasing values of K
  for(int i = 0; i < kLimit; i++)
  {
    double rxSquare = rx*rx;
    double n = MIN(i + 1, pow(rTemp/rx,(double)d));
    double error = epsilon + 1;
    double temp = 1;
    int p = 0;
    /*
    while((error > epsilon) & (p <= P_UPPER_LIMIT))
    {
      p++;
      double b = MIN(((rx + sqrt((rxSquare) + (2*p*hSquare)))/2), rx + rTemp);
      double c = rx - b;
      temp = temp*(((2*rx*b)/hSquare)/p);
      error = temp*(exp(-(c*c)/hSquare));      
    }*/  
    //double memTotalC = 8*(i+1);  // 8 is number of bytes, i+1 is K
    while((error > epsilon) && (p <= P_UPPER_LIMIT) ) //&& (memTotalC <= C_UPPER_LIMIT) )
    {
      p++;
      double b = MIN(((rx + sqrt((rxSquare) + (2*p*hSquare)))/2), rx + rTemp);
      double c = rx - b;
      temp = 1;
      for( int j = 1; j <= p; j++ )
        temp = temp*((2.0*rx*b/hSquare)/j);
      error = temp*(exp(-(c*c)/hSquare));

      //memTotalC *= (d+p);
      //memTotalC /= p;
    }
    double complexity = d*(i + 1) + d*log((double)i + 1) + ((1 + n)*nchoosek_double(p - 1 + d, d));
    if ( (complexity < complexityMin) && (error <= epsilon))
    {
      complexityMin = complexity;
      kTemp = i + 1;
      pMaxTemp = p;  
      errorTemp = error;
    }
    
    // try to guess if we have gone past the minimum (the complexity function
    // zigzags as we increase number of clusters, but if it goes up enough,
    // we'll assume we've passed the global minimum).
    // Also stop if truncation number is only 1 or if the max number of unique 
    // clusters are reached (rx = 0).
    double nextComplexityEstimate = d*(i + 1) + d*log((double)i + 1) + ((1 + n)*nchoosek_double(p - 2 + d, d));   
    if( (p == 1) || (rx <= 0) || ( nextComplexityEstimate > 2*complexityMin || complexity > 2*complexityMin ) )
    {
      break;    
    }

    // add another cluster center, and get new max cluster radius
    kcc->ClusterIncrement( &numClusters, &rx );
    complexityLast = complexity;
  }  

  // added this to catch case where desired error is never reached.
  // The best thing is to have as many clusters and terms in the taylor
  // series as possible (which will give lowest error)
  if( errorTemp > epsilon )
  {
    kTemp = kLimit;
  }

  //printf("memLimit = %e, memTotalC = %e\n", (double)C_UPPER_LIMIT, kTemp*8*nchoosek_double(pMaxTemp-1+d,d));

  // copy results 
  if( K != NULL )
    *K = kTemp;
  if( pMax != NULL )
    *pMax = pMaxTemp;
  if( r != NULL )
    *r = rTemp;
  if( errorBound != NULL )
    *errorBound = errorTemp;

  delete [] pClusterIndex;
  delete kcc;

  return 0;
}

//------------------------------------------------------------------------------
// This function is an interface to the C++ KCenterClustering class from the
// original IFGT library.
//
// Created by Vlad Morariu 2007-06-19.
//------------------------------------------------------------------------------
int figtreeKCenterClustering( int d, int N, double * x, int kMax, int * K,
                           double * rx, int * clusterIndex, double * clusterCenters, 
                           int * numPoints, double * clusterRadii )
{
  // check input arguments
  FIGTREE_CHECK_POS_NONZERO_INT( d, figtreeKCenterClustering );
  FIGTREE_CHECK_POS_NONZERO_INT( N, figtreeKCenterClustering );
  FIGTREE_CHECK_NONNULL_PTR( x, figtreeKCenterClustering );
  FIGTREE_CHECK_POS_NONZERO_INT( kMax, figtreeKCenterClustering );
  FIGTREE_CHECK_NONNULL_PTR( K, figtreeKCenterClustering );
  FIGTREE_CHECK_NONNULL_PTR( rx, figtreeKCenterClustering );
  FIGTREE_CHECK_NONNULL_PTR( clusterIndex, figtreeKCenterClustering );
  FIGTREE_CHECK_NONNULL_PTR( clusterCenters, figtreeKCenterClustering );
  FIGTREE_CHECK_NONNULL_PTR( numPoints, figtreeKCenterClustering );
  FIGTREE_CHECK_NONNULL_PTR( clusterRadii, figtreeKCenterClustering );

  //k-center clustering
  KCenterClustering* pKCC = new KCenterClustering( d, N, x, clusterIndex, kMax );
  *K = pKCC->Cluster();
  if( rx != NULL )
    *rx = pKCC->MaxClusterRadius;
  pKCC->ComputeClusterCenters(*K, clusterCenters, numPoints, clusterRadii);

  delete pKCC;
  return 0;
}

//------------------------------------------------------------------------------
// Actual function to evaluate the exact Gauss Transform directly.
// Originally Evaluate() from GaussTransform.cpp, written by Vikas C. Raykar.
//
// Modified by Vlad Morariu on 2007-06-19.
//------------------------------------------------------------------------------
int figtreeEvaluateDirect( int d, int N, int M, double * x, double h, 
                        double * q, double * y, double * g )
{
  // check input arguments
  FIGTREE_CHECK_POS_NONZERO_INT( d, figtreeEvaluateDirect );
  FIGTREE_CHECK_POS_NONZERO_INT( N, figtreeEvaluateDirect );
  FIGTREE_CHECK_POS_NONZERO_INT( M, figtreeEvaluateDirect );
  FIGTREE_CHECK_NONNULL_PTR( x, figtreeEvaluateDirect );
  FIGTREE_CHECK_POS_NONZERO_DOUBLE( h, figtreeEvaluateDirect );
  FIGTREE_CHECK_NONNULL_PTR( q, figtreeEvaluateDirect );
  FIGTREE_CHECK_NONNULL_PTR( y, figtreeEvaluateDirect );
  FIGTREE_CHECK_NONNULL_PTR( g, figtreeEvaluateDirect );

  // evaluate
  double hSquare = h*h;
  for(int j = 0; j < M; j++)
  {
    g[j] = 0.0;
    for(int i = 0; i < N; i++)
    {
      double norm = 0.0;
      for (int k = 0; k < d; k++)
      {
        double temp = x[(d*i) + k] - y[(d*j) + k];
        norm = norm + (temp*temp);
      }
      g[j] = g[j] + (q[i]*exp(-norm/hSquare));
    }
  }

  return 0;
}

//------------------------------------------------------------------------------
// This function approximates Gauss Transform.
// Originally constructor, Evaluate(), and destructor from 
// ImprovedFastGaussTransform.cpp (IFGT source code).
//
// Modified by Vlad Morariu on 2007-06-19.
//------------------------------------------------------------------------------
int figtreeEvaluateIfgt( int d, int N, int M, int W, double * x, 
                           double h, double * q, double * y, 
                           int pMax, int K, int * clusterIndex, 
                           double * clusterCenter, double * clusterRadii,
                           double r, double epsilon, double * g )
{
  // check input arguments
  FIGTREE_CHECK_POS_NONZERO_INT( d, figtreeEvaluateIfgt );
  FIGTREE_CHECK_POS_NONZERO_INT( N, figtreeEvaluateIfgt );
  FIGTREE_CHECK_POS_NONZERO_INT( M, figtreeEvaluateIfgt );
  FIGTREE_CHECK_POS_NONZERO_INT( W, figtreeEvaluateIfgt );
  FIGTREE_CHECK_NONNULL_PTR( x, figtreeEvaluateIfgt );
  FIGTREE_CHECK_POS_NONZERO_DOUBLE( h, figtreeEvaluateIfgt );
  FIGTREE_CHECK_NONNULL_PTR( q, figtreeEvaluateIfgt );
  FIGTREE_CHECK_NONNULL_PTR( y, figtreeEvaluateIfgt );
  FIGTREE_CHECK_POS_NONZERO_INT( pMax, figtreeEvaluateIfgt );
  FIGTREE_CHECK_POS_NONZERO_INT( K, figtreeEvaluateIfgt );
  FIGTREE_CHECK_NONNULL_PTR( clusterIndex, figtreeEvaluateIfgt );
  FIGTREE_CHECK_NONNULL_PTR( clusterCenter, figtreeEvaluateIfgt );
  FIGTREE_CHECK_NONNULL_PTR( clusterRadii, figtreeEvaluateIfgt );
  FIGTREE_CHECK_POS_NONZERO_DOUBLE( r, figtreeEvaluateIfgt );
  FIGTREE_CHECK_POS_NONZERO_DOUBLE( epsilon, figtreeEvaluateIfgt );
  FIGTREE_CHECK_NONNULL_PTR( g, figtreeEvaluateIfgt );

  //Memory allocation
  int pMaxTotal = nchoosek(pMax - 1 + d, d);
  double hSquare=h*h;
  double * targetCenterMonomials = new double[pMaxTotal];
  double * dy = new double[d];
  double * C = new double[W*K*pMaxTotal];
  double * ry = new double[K];
  double * rySquare = new double[K];

  for(int i = 0; i < K; i++)
  {
    ry[i] = r + clusterRadii[i];
    rySquare[i] = ry[i]*ry[i];
  }   

  //////////////////////////////////////////////////////////////////////////////
  // Evaluate  
  //////////////////////////////////////////////////////////////////////////////
  computeC( d, N, W, K, pMaxTotal, pMax, h, clusterIndex, x, q, clusterCenter, C );  

  for(int j = 0; j < M; j++)
  {
    for( int w = 0; w < W; w++ )
    {
      g[M*w + j] = 0.0;
    }

    int targetBase = j*d;        
    for(int k = 0; k < K; k++)
    {
      int centerBase = k*d;
      double targetCenterDistanceSquare = 0.0;
      for(int i = 0; i < d; i++)
      {
        dy[i] = y[targetBase + i] - clusterCenter[centerBase + i];
        targetCenterDistanceSquare += dy[i]*dy[i];
        if(targetCenterDistanceSquare > rySquare[k]) break;
      }

      if(targetCenterDistanceSquare <= rySquare[k])
      {
        computeTargetCenterMonomials( d, h, dy, pMax, targetCenterMonomials );
        double f=exp(-targetCenterDistanceSquare/hSquare);
        for(int w = 0; w < W; w++ )
        {        
          for(int alpha = 0; alpha < pMaxTotal; alpha++)
          {
            g[M*w + j] += (C[(K*w + k)*pMaxTotal + alpha]*f*targetCenterMonomials[alpha]);
          }
        }                      
      }
    }
  }

  //////////////////////////////////////////////////////////////////////////////
  // Release memory
  //////////////////////////////////////////////////////////////////////////////
  delete [] rySquare;
  delete [] ry;
  delete [] C;
  delete [] dy;
  delete [] targetCenterMonomials;

  return 0;
}


//------------------------------------------------------------------------------
// This function evaluates the IFGT by using a different for each source assuming
//   the worst-case placement of any target and for each target assuming the 
//   worst-case placement of any source.  Because the error bound is guaranteed
//   to be satisfied for each source-target point pair, it is a point-wise adaptive
//   version of the IFGT.
//
// See 'Automatic online tuning for fast Gaussian summation,' by Morariu et al, 
//   NIPS 2008 for details.  
//
// NOTES: Unlike the cluster-wise adaptive version, this does work for W>1.
//
// Created by Vlad Morariu on 2008-06-04.
//------------------------------------------------------------------------------
int figtreeEvaluateIfgtAdaptivePoint( int d, int N, int M, int W, double * x, 
                                      double h, double * q, double * y, 
                                      int pMax, int K, int * clusterIndex, 
                                      double * clusterCenter, double * clusterRadii,
                                      double r, double epsilon,
                                      double * g )
{
  // check input arguments
  FIGTREE_CHECK_POS_NONZERO_INT( d, figtreeEvaluateIfgtAdaptive );
  FIGTREE_CHECK_POS_NONZERO_INT( N, figtreeEvaluateIfgtAdaptive );
  FIGTREE_CHECK_POS_NONZERO_INT( M, figtreeEvaluateIfgtAdaptive );
  FIGTREE_CHECK_POS_NONZERO_INT( W, figtreeEvaluateIfgtAdaptive );
  FIGTREE_CHECK_NONNULL_PTR( x, figtreeEvaluateIfgtIfgtAdaptive );
  FIGTREE_CHECK_POS_NONZERO_DOUBLE( h, figtreeEvaluateIfgtAdaptive );
  FIGTREE_CHECK_NONNULL_PTR( q, figtreeEvaluateIfgtAdaptive );
  FIGTREE_CHECK_NONNULL_PTR( y, figtreeEvaluateIfgtAdaptive );
  FIGTREE_CHECK_POS_NONZERO_INT( pMax, figtreeEvaluateIfgtAdaptive );
  FIGTREE_CHECK_POS_NONZERO_INT( K, figtreeEvaluateIfgtAdaptive );
  FIGTREE_CHECK_NONNULL_PTR( clusterIndex, figtreeEvaluateIfgtAdaptive );
  FIGTREE_CHECK_NONNULL_PTR( clusterCenter, figtreeEvaluateIfgtAdaptive );
  FIGTREE_CHECK_NONNULL_PTR( clusterRadii, figtreeEvaluateIfgtAdaptive );
  FIGTREE_CHECK_POS_NONZERO_DOUBLE( r, figtreeEvaluateIfgtAdaptive );
  FIGTREE_CHECK_POS_NONZERO_DOUBLE( epsilon, figtreeEvaluateIfgtAdaptive );
  FIGTREE_CHECK_NONNULL_PTR( g, figtreeEvaluateIfgtAdaptive );

  //Memory allocation
  int pMaxTotal = nchoosek(pMax - 1 + d, d);
  int * pMaxTotals = new int[pMax];
  for( int i = 0; i < pMax; i++ )
    pMaxTotals[i] = nchoosek( i + d, d );

  double hSquare=h*h;
  double * targetCenterMonomials = new double[pMaxTotal];
  double * dy = new double[d];
  double * C = new double[W*K*pMaxTotal];
  double * ry = new double[K];
  double * rySquare = new double[K];

  double rx = clusterRadii[0];
  for(int i = 0; i < K; i++)
  {
    ry[i] = r + clusterRadii[i];
    rySquare[i] = ry[i]*ry[i];
    rx = MAX( rx, clusterRadii[i] );
  } 

  //////////////////////////////////////////////////////////////////////////////
  // Evaluate  
  //////////////////////////////////////////////////////////////////////////////
  // for each cluster, compute max distances at which we can use a certain truncation number
  double * maxSourceDists2 = new double[pMax];
  figtreeSourceTruncationRanges( r, rx, h, epsilon, pMax, maxSourceDists2 );
  computeCAdaptivePoint( d, N, W, K, pMaxTotal, pMax, h, clusterIndex, x, q, clusterCenter, maxSourceDists2, pMaxTotals, C );
  delete [] maxSourceDists2;

  // for each cluster, compute distance ranges for each truncation number
  double * targetDists2 = new double[2*pMax];
  figtreeTargetTruncationRanges( r, rx, h, epsilon, pMax, targetDists2, targetDists2+pMax );

  //int * pHistogram = new int[pMax];
  //memset(pHistogram,0,sizeof(int)*pMax);

  memset( g, 0, sizeof(double)*M*W );
    //int targetBase = j*d;        
  for(int k = 0; k < K; k++)
  {
    for(int j = 0; j < M; j++)
    {
      //int centerBase = k*d;
      double targetCenterDistanceSquare = 0.0;
      for(int i = 0; i < d; i++)
      {
        dy[i] = y[j*d + i] - clusterCenter[k*d + i];
        targetCenterDistanceSquare += dy[i]*dy[i];
        if(targetCenterDistanceSquare > rySquare[k]) break;
      }

      if(targetCenterDistanceSquare <= rySquare[k])
      {
        int p = figtreeTargetTruncationNumber( targetCenterDistanceSquare, pMax, targetDists2, targetDists2+pMax );
        int pTotal = pMaxTotals[p-1];
        //pHistogram[p-1]++;
        computeTargetCenterMonomials( d, h, dy, p, targetCenterMonomials );
        double f=exp(-targetCenterDistanceSquare/hSquare);
        for(int w = 0; w < W; w++ )
        {        
          double * C_offset = C + (K*w + k)*pMaxTotal;
          for(int alpha = 0; alpha < pTotal; alpha++)
          {
            g[M*w + j] += *(C_offset++)*f*targetCenterMonomials[alpha];
          }
        }                      
      }
    }
  }

  //printf( "target p histogram: ");
  //for( int i = 0; i < pMax; i++ )
  //  printf( " %i", pHistogram[i] );
  //printf( "\n" );
  //delete [] pHistogram;

  //////////////////////////////////////////////////////////////////////////////
  // Release memory
  //////////////////////////////////////////////////////////////////////////////
  delete [] rySquare;
  delete [] ry;
  delete [] C;
  delete [] dy;
  delete [] targetCenterMonomials;

  delete [] targetDists2;
  delete [] pMaxTotals;

  return 0;
}

//------------------------------------------------------------------------------
// This function evaluates the IFGT by using a different truncation number for
//   each cluster to ensure that the total error contribution from each
//   cluster satisfies the error bound.  Thus, this is the cluster-wise adaptive
//   version of the IFGT.  
//   
// See 'Automatic online tuning for fast Gaussian summation,' by Morariu et al, 
//   NIPS 2008 for details.  
//
// NOTES: 1) Currently is not implemented to work for W>1.  For W>1, use the
//           point-wise adaptive version instead.
//
//        2) The method could be extended to use different truncation numbers 
//           for each target by splitting targets into concentric regions and
//           computing the cluster-wise truncation for each concentric region
//           separately (because each region will have a different max distance 
//           from the cluster center, the truncations will differ by region).  
//           However, the current implementation uses only one region for the 
//           targets, and varies the truncation by cluster for that region.
//
// Created by Vlad Morariu on 2008-06-04.
//------------------------------------------------------------------------------
int figtreeEvaluateIfgtAdaptiveCluster( int d, int N, int M, int W, double * x, 
                                        double h, double * q, double * y, 
                                        int pMax, int K, int * clusterIndex, 
                                        double * clusterCenter, double * clusterRadii,
                                        double r, double epsilon, int * clusterTruncations,
                                        double * g )
{
  // check input arguments
  FIGTREE_CHECK_POS_NONZERO_INT( d, figtreeEvaluateIfgtAdaptive );
  FIGTREE_CHECK_POS_NONZERO_INT( N, figtreeEvaluateIfgtAdaptive );
  FIGTREE_CHECK_POS_NONZERO_INT( M, figtreeEvaluateIfgtAdaptive );
  FIGTREE_CHECK_POS_NONZERO_INT( W, figtreeEvaluateIfgtAdaptive );
  FIGTREE_CHECK_NONNULL_PTR( x, figtreeEvaluateIfgtIfgtAdaptive );
  FIGTREE_CHECK_POS_NONZERO_DOUBLE( h, figtreeEvaluateIfgtAdaptive );
  FIGTREE_CHECK_NONNULL_PTR( q, figtreeEvaluateIfgtAdaptive );
  FIGTREE_CHECK_NONNULL_PTR( y, figtreeEvaluateIfgtAdaptive );
  FIGTREE_CHECK_POS_NONZERO_INT( pMax, figtreeEvaluateIfgtAdaptive );
  FIGTREE_CHECK_POS_NONZERO_INT( K, figtreeEvaluateIfgtAdaptive );
  FIGTREE_CHECK_NONNULL_PTR( clusterIndex, figtreeEvaluateIfgtAdaptive );
  FIGTREE_CHECK_NONNULL_PTR( clusterCenter, figtreeEvaluateIfgtAdaptive );
  FIGTREE_CHECK_NONNULL_PTR( clusterRadii, figtreeEvaluateIfgtAdaptive );
  FIGTREE_CHECK_POS_NONZERO_DOUBLE( r, figtreeEvaluateIfgtAdaptive );
  FIGTREE_CHECK_POS_NONZERO_DOUBLE( epsilon, figtreeEvaluateIfgtAdaptive );
  FIGTREE_CHECK_NONNULL_PTR( g, figtreeEvaluateIfgtAdaptive );

  //Memory allocation
  int pMaxTotal = nchoosek(pMax - 1 + d, d);
  int * pMaxTotals = new int[pMax];
  for( int i = 0; i < pMax; i++ )
    pMaxTotals[i] = nchoosek( i + d, d );

  double hSquare=h*h;
  double * targetCenterMonomials = new double[pMaxTotal];
  double * dy = new double[d];
  double * C = new double[W*K*pMaxTotal];
  double * ry = new double[K];
  double * rySquare = new double[K];

  double rx = clusterRadii[0];
  for(int i = 0; i < K; i++)
  {
    ry[i] = r + clusterRadii[i];
    rySquare[i] = ry[i]*ry[i];
    rx = MAX( rx, clusterRadii[i] );
  } 

  //////////////////////////////////////////////////////////////////////////////
  // Evaluate  
  //////////////////////////////////////////////////////////////////////////////
  // for each cluster, compute max distances at which we can use a certain truncation number
  computeCAdaptiveCluster( d, N, W, K, pMaxTotal, pMax, h, clusterIndex, x, q, clusterCenter, clusterTruncations, pMaxTotals, C );

  memset( g, 0, sizeof(double)*M*W );
    //int targetBase = j*d;        
  for(int k = 0; k < K; k++)
  {
    int p = clusterTruncations[k];
    int pTotal = pMaxTotals[p-1];
    for(int j = 0; j < M; j++)
    {
      //int centerBase = k*d;
      double targetCenterDistanceSquare = 0.0;
      for(int i = 0; i < d; i++)
      {
        dy[i] = y[j*d + i] - clusterCenter[k*d + i];
        targetCenterDistanceSquare += dy[i]*dy[i];
        if(targetCenterDistanceSquare > rySquare[k]) break;
      }

      if(targetCenterDistanceSquare <= rySquare[k])
      {
        computeTargetCenterMonomials( d, h, dy, p, targetCenterMonomials );
        double f=exp(-targetCenterDistanceSquare/hSquare);
        for(int w = 0; w < W; w++ )
        {        
          double * C_offset = C + (K*w + k)*pMaxTotal;
          for(int alpha = 0; alpha < pTotal; alpha++)
          {
            g[M*w + j] += *(C_offset++)*f*targetCenterMonomials[alpha];
          }
        }                      
      }
    }
  }

  //////////////////////////////////////////////////////////////////////////////
  // Release memory
  //////////////////////////////////////////////////////////////////////////////
  delete [] rySquare;
  delete [] ry;
  delete [] C;
  delete [] dy;
  delete [] targetCenterMonomials;

  delete [] pMaxTotals;

  return 0;
}

//------------------------------------------------------------------------------
// This function approximates Gauss Transform using Approximate Nearest 
// Neighbors.
// Originally constructor, Evaluate(), and destructor from 
// ImprovedFastGaussTransform.cpp of FIGTree code, by Vikas C. Raykar.
//
// Modified by Vlad Morariu on 2007-06-19.
//------------------------------------------------------------------------------
int figtreeEvaluateIfgtTree( int d, int N, int M, int W, double * x, 
                              double h, double * q, double * y, 
                              int pMax, int K, int * clusterIndex, 
                              double * clusterCenter, double * clusterRadii,
                              double r, double epsilon, double * g )
{
#ifdef FIGTREE_NO_ANN
  printf("This code was not compiled with support for ANN.  Please recompile\n");
  printf("with 'FIGTREE_NO_ANN' not defined to enable ANN support.\n");
  return -1;
#else
  // check input arguments
  FIGTREE_CHECK_POS_NONZERO_INT( d, figtreeEvaluateIfgtTree );
  FIGTREE_CHECK_POS_NONZERO_INT( N, figtreeEvaluateIfgtTree );
  FIGTREE_CHECK_POS_NONZERO_INT( M, figtreeEvaluateIfgtTree );
  FIGTREE_CHECK_POS_NONZERO_INT( W, figtreeEvaluateIfgtTree );
  FIGTREE_CHECK_NONNULL_PTR( x, figtreeEvaluateIfgtTree );
  FIGTREE_CHECK_POS_NONZERO_DOUBLE( h, figtreeEvaluateIfgtTree );
  FIGTREE_CHECK_NONNULL_PTR( q, figtreeEvaluateIfgtTree );
  FIGTREE_CHECK_NONNULL_PTR( y, figtreeEvaluateIfgtTree );
  FIGTREE_CHECK_POS_NONZERO_INT( pMax, figtreeEvaluateIfgtTree );
  FIGTREE_CHECK_POS_NONZERO_INT( K, figtreeEvaluateIfgtTree );
  FIGTREE_CHECK_NONNULL_PTR( clusterIndex, figtreeEvaluateIfgtTree );
  FIGTREE_CHECK_NONNULL_PTR( clusterCenter, figtreeEvaluateIfgtTree );
  FIGTREE_CHECK_NONNULL_PTR( clusterRadii, figtreeEvaluateIfgtTree );
  FIGTREE_CHECK_POS_NONZERO_DOUBLE( r, figtreeEvaluateIfgtTree );
  FIGTREE_CHECK_POS_NONZERO_DOUBLE( epsilon, figtreeEvaluateIfgtTree );
  FIGTREE_CHECK_NONNULL_PTR( g, figtreeEvaluateIfgtTree );

  //Memory allocation
  int pMaxTotal = nchoosek(pMax-1+d,d);
  double * targetCenterMonomials = new double[pMaxTotal];
  double * dy = new double[d];
  double * C = new double[W*K*pMaxTotal]; 
  double hSquare = h*h;

  //Find the maximum cluster radius
  double pcr_max = clusterRadii[0];
  for(int i = 0; i < K; i++)
  {
    if (clusterRadii[i] > pcr_max)
    {
      pcr_max = clusterRadii[i];
    }
  }
  double rSquare=(r+pcr_max)*(r+pcr_max);

  //Allocate storage using ANN procedures 
  ANNpointArray dataPts = annAllocPts(K,d);     // allocate data points
  ANNidxArray   nnIdx   = new ANNidx[K];        // allocate near neigh indices
  ANNdistArray  dists = new ANNdist[K];         // allocate near neighbor dists
  
  // Copy the cluster centers to the ANN data structure
  for (int k = 0; k < K; k++)
  {
    for ( int j = 0; j < d; j++ )
      dataPts[k][j]= clusterCenter[k*d + j];
  }

  // build search structure
  ANNkd_tree * kdTree = new ANNkd_tree(              
                                        dataPts,  // the data points
                                        K,        // number of points
                                        d,        // dimension of space
                                        1,
                                        ANN_KD_SUGGEST);

  ////////////////////////////////////////////////////////////////////
  // Evaluate
  ////////////////////////////////////////////////////////////////////
  computeC( d, N, W, K, pMaxTotal, pMax, h, clusterIndex, x, q, clusterCenter, C );  

  for(int j = 0; j < M; j++)
  {
    for( int w = 0; w < W; w++ )
    {
      g[M*w+j]=0.0;
    }

    int targetBase=j*d;        
    
    ANNpoint queryPt=&(y[targetBase]);

    int NN = kdTree->annkFRSearchUnordered(           // search
                                   queryPt,  // query point
                                   rSquare,  // squared radius
                                   N,        // number of near neighbors
                                   nnIdx,     // nearest neighbors (returned)
                                   dists,     // distance (returned)
                                   0.0 );

    if (NN>0)
    {
      for(int l = 0; l < NN; l++)
      {
        int k = nnIdx[l];
        int centerBase = k*d;
        double  targetCenterDistanceSquare = dists[l];
        for(int i = 0; i < d; i++)
        {
          dy[i] = y[targetBase + i] - clusterCenter[centerBase + i];
        }
        computeTargetCenterMonomials( d, h, dy, pMax, targetCenterMonomials );
        double e = exp(-targetCenterDistanceSquare/hSquare);
        for(int w = 0; w < W; w++ )
        {        
          for(int alpha = 0; alpha < pMaxTotal; alpha++)
          {
            g[M*w + j] += (C[(K*w+k)*pMaxTotal + alpha]*e*targetCenterMonomials[alpha]);
          }    
        }
      }
    }
  }

  ////////////////////////////////////////////////////////////////////
  // Release Memory
  ////////////////////////////////////////////////////////////////////
  delete [] targetCenterMonomials;
  delete [] dy;
  delete [] C;

  annDeallocPts(dataPts);
  delete [] nnIdx;              
  delete [] dists;
  delete kdTree;
  annClose();        

  return 0;
#endif
}

//------------------------------------------------------------------------------
// Gauss Transform computed using the ANN library.
// Given a specified epsilon, the code computes the Gauss transform by summing 
// the sources only within a certain radius--whose contribution is at least 
// epsilon. The neighbors are found using the ANN library.
// http://www.cs.umd.edu/~mount/ANN/.
// Originally constructor, Evaluate(), and destructor from GaussTransformTree.cpp 
// of FIGTree code, by Vikas C. Raykar.
//
// Modified by Vlad Morariu on 2007-06-20
// Modified by Vlad Morariu on 2008-06-10 to use ANN fixed radius search with
//   nearest neighbors unordered (saves a significant amt of time), and with
//   one call instead of two, as a result.
//------------------------------------------------------------------------------
int figtreeEvaluateDirectTree( int d, int N, int M, double * x, double h, 
                               double * q, double * y, double epsilon, double * g )
{
#ifdef FIGTREE_NO_ANN
  printf("This code was not compiled with support for ANN.  Please recompile\n");
  printf("with compiler flag 'FIGTREE_NO_ANN' not set to enable ANN support.\n");
  return -1;
#else
  // check input arguments 
  FIGTREE_CHECK_POS_NONZERO_INT( d, figtreeEvaluateDirectTreeUnordered );
  FIGTREE_CHECK_POS_NONZERO_INT( N, figtreeEvaluateDirectTreeUnordered );
  FIGTREE_CHECK_POS_NONZERO_INT( M, figtreeEvaluateDirectTreeUnordered );
  FIGTREE_CHECK_NONNULL_PTR( x, figtreeEvaluateDirectTreeUnordered );
  FIGTREE_CHECK_POS_NONZERO_DOUBLE( h, figtreeEvaluateDirectTreeUnordered );
  FIGTREE_CHECK_NONNULL_PTR( q, figtreeEvaluateDirectTreeUnordered );
  FIGTREE_CHECK_NONNULL_PTR( y, figtreeEvaluateDirectTreeUnordered );
  FIGTREE_CHECK_POS_NONZERO_DOUBLE( epsilon, figtreeEvaluateDirectTreeUnordered );
  FIGTREE_CHECK_NONNULL_PTR( g, figtreeEvaluateDirectTreeUnordered );

  double hSquare = h*h;
  double epsANN = 0.0;

  // Compute the cutoff radius
  double r = h*sqrt(log(1/epsilon));
  double rSquare=r*r;

  // Allocate storage using ANN procedures 
  ANNpointArray dataPts = annAllocPts(N,d);  // allocate data points
  ANNidxArray   nnIdx   = new ANNidx[N];     // allocate near neigh indices
  ANNdistArray  dists   = new ANNdist[N];    // allocate near neighbor dists
  
  // Copy the source points to the ANN data structure
  for (int i = 0; i < N; i++)
  {
    for ( int j = 0; j < d; j++ )
      dataPts[i][j]= x[i*d + j];
  }

  // build search structure
  ANNkd_tree * kdTree = new ANNkd_tree(              
                                        dataPts, // the data points
                                        N,       // number of points
                                        d,       // dimension of space
                                        1,
                                        ANN_KD_SUGGEST );

  ///////////////////////////////////////////////////////////////////////
  // Evaluate
  ///////////////////////////////////////////////////////////////////////
  for(int j = 0; j < M; j++)
  {
    g[j] = 0.0;  
    int targetBase = j*d;          
    ANNpoint queryPt = &(y[targetBase]);

    int NN = kdTree->annkFRSearchUnordered( // fixed radius search
                            queryPt,        // query point
                            rSquare,        // squared radius
                            N,              // number of near neighbors
                            nnIdx,          // nearest neighbors (returned)
                            dists,          // distance (returned)
                            epsANN );
    for(int l = 0; l < NN; l++)
    {
      int i = nnIdx[l];
      double sourceTargetDistanceSquare = dists[l];
      g[j] += (q[i]*exp(-sourceTargetDistanceSquare/hSquare));
    }
  }

  //////////////////////////////////////////////////////////////////////////////
  // Free memory
  //////////////////////////////////////////////////////////////////////////////
  annDeallocPts(dataPts);
  delete [] nnIdx;              
  delete [] dists;
  delete kdTree;
  annClose();   

  return 0;
#endif
}


//------------------------------------------------------------------------------
// This function evaluates the IFGT by using a different for each source assuming
//   the worst-case placement of any target and for each target assuming the 
//   worst-case placement of any source.  Because the error bound is guaranteed
//   to be satisfied for each source-target point pair, it is a point-wise adaptive
//   version of the IFGT.
//
// This function uses a tree for finding nearby cluster centers.
//
// See 'Automatic online tuning for fast Gaussian summation,' by Morariu et al, 
//   NIPS 2008 for details.  
//
// NOTES: Unlike the cluster-wise adaptive version, this does work for W>1.
//
// Created by Vlad Morariu on 2008-06-04.
//------------------------------------------------------------------------------
int figtreeEvaluateIfgtTreeAdaptivePoint( int d, int N, int M, int W, double * x, 
                                      double h, double * q, double * y, 
                                      int pMax, int K, int * clusterIndex, 
                                      double * clusterCenter, double * clusterRadii,
                                      double r, double epsilon,
                                      double * g )
{
#ifdef FIGTREE_NO_ANN
  printf("This code was not compiled with support for ANN.  Please recompile\n");
  printf("with compiler flag 'FIGTREE_NO_ANN' not set to enable ANN support.\n");
  return -1;
#else
  // check input arguments
  FIGTREE_CHECK_POS_NONZERO_INT( d, figtreeEvaluateIfgtTreeAdaptivePoint );
  FIGTREE_CHECK_POS_NONZERO_INT( N, figtreeEvaluateIfgtTreeAdaptivePoint );
  FIGTREE_CHECK_POS_NONZERO_INT( M, figtreeEvaluateIfgtTreeAdaptivePoint );
  FIGTREE_CHECK_POS_NONZERO_INT( W, figtreeEvaluateIfgtTreeAdaptivePoint );
  FIGTREE_CHECK_NONNULL_PTR( x, figtreeEvaluateIfgtIfgtTreeAdaptivePoint );
  FIGTREE_CHECK_POS_NONZERO_DOUBLE( h, figtreeEvaluateIfgtTreeAdaptivePoint );
  FIGTREE_CHECK_NONNULL_PTR( q, figtreeEvaluateIfgtTreeAdaptivePoint );
  FIGTREE_CHECK_NONNULL_PTR( y, figtreeEvaluateIfgtTreeAdaptivePoint );
  FIGTREE_CHECK_POS_NONZERO_INT( pMax, figtreeEvaluateIfgtTreeAdaptivePoint );
  FIGTREE_CHECK_POS_NONZERO_INT( K, figtreeEvaluateIfgtTreeAdaptivePoint );
  FIGTREE_CHECK_NONNULL_PTR( clusterIndex, figtreeEvaluateIfgtTreeAdaptivePoint );
  FIGTREE_CHECK_NONNULL_PTR( clusterCenter, figtreeEvaluateIfgtTreeAdaptivePoint );
  FIGTREE_CHECK_NONNULL_PTR( clusterRadii, figtreeEvaluateIfgtTreeAdaptivePoint );
  FIGTREE_CHECK_POS_NONZERO_DOUBLE( r, figtreeEvaluateIfgtTreeAdaptivePoint );
  FIGTREE_CHECK_POS_NONZERO_DOUBLE( epsilon, figtreeEvaluateIfgtTreeAdaptivePoint );
  FIGTREE_CHECK_NONNULL_PTR( g, figtreeEvaluateIfgtTreeAdaptivePoint );

  //Memory allocation
  int pMaxTotal = nchoosek(pMax - 1 + d, d);
  int * pMaxTotals = new int[pMax];
  for( int i = 0; i < pMax; i++ )
    pMaxTotals[i] = nchoosek( i + d, d );

  double hSquare=h*h;
  double * targetCenterMonomials = new double[pMaxTotal];
  double * dy = new double[d];
  double * C = new double[W*K*pMaxTotal];
  double * ry = new double[K];
  double * rySquare = new double[K];

  double rx = clusterRadii[0];
  for(int i = 0; i < K; i++)
  {
    ry[i] = r + clusterRadii[i];
    rySquare[i] = ry[i]*ry[i];
    rx = MAX( rx, clusterRadii[i] );
  } 

  //
  // Build tree on cluster centers
  //
  double rSquare = (r+rx)*(r+rx);
  //Allocate storage using ANN procedures 
  ANNpointArray dataPts = annAllocPts(K,d);     // allocate data points
  ANNidxArray   nnIdx   = new ANNidx[K];        // allocate near neigh indices
  ANNdistArray  dists = new ANNdist[K];         // allocate near neighbor dists
  
  // Copy the cluster centers to the ANN data structure
  for (int k = 0; k < K; k++)
  {
    for ( int j = 0; j < d; j++ )
      dataPts[k][j]= clusterCenter[k*d + j];
  }

  // build search structure
  ANNkd_tree * kdTree = new ANNkd_tree(              
                                        dataPts,  // the data points
                                        K,        // number of points
                                        d,        // dimension of space
                                        1,
                                        ANN_KD_SUGGEST);


  //////////////////////////////////////////////////////////////////////////////
  // Evaluate  
  //////////////////////////////////////////////////////////////////////////////
  // for each cluster, compute max distances at which we can use a certain truncation number
  double * maxSourceDists2 = new double[pMax];
  figtreeSourceTruncationRanges( r, rx, h, epsilon, pMax, maxSourceDists2 );
  computeCAdaptivePoint( d, N, W, K, pMaxTotal, pMax, h, clusterIndex, x, q, clusterCenter, maxSourceDists2, pMaxTotals, C );
  delete [] maxSourceDists2;

  // for each cluster, compute distance ranges for each truncation number
  double * targetDists2 = new double[2*pMax];
  figtreeTargetTruncationRanges( r, rx, h, epsilon, pMax, targetDists2, targetDists2+pMax );

  //int * pHistogram = new int[pMax];
  //memset(pHistogram,0,sizeof(int)*pMax);

  memset( g, 0, sizeof(double)*M*W );

  for(int j = 0; j < M; j++)
  {
    int targetBase = j*d;
    ANNpoint queryPt=&(y[targetBase]);

    int NN = kdTree->annkFRSearchUnordered(           // search
                                   queryPt,  // query point
                                   rSquare,  // squared radius
                                   K,        // number of near neighbors
                                   nnIdx,     // nearest neighbors (returned)
                                   dists,     // distance (returned)
                                   0.0 );
    for(int l = 0; l < NN; l++)
    {
      int k = nnIdx[l];

      int centerBase = k*d;
      double targetCenterDistanceSquare = dists[l];
      if(targetCenterDistanceSquare <= rySquare[k])
      {
        int p = figtreeTargetTruncationNumber( targetCenterDistanceSquare, pMax, targetDists2, targetDists2+pMax );
        int pTotal = pMaxTotals[p-1];
        //pHistogram[p-1]++;
        for(int i = 0; i < d; i++)
        {
          dy[i] = y[targetBase + i] - clusterCenter[centerBase + i];
        }
        computeTargetCenterMonomials( d, h, dy, p, targetCenterMonomials );
        double f=exp(-targetCenterDistanceSquare/hSquare);
        for(int w = 0; w < W; w++ )
        {        
          double * C_offset = C + (K*w + k)*pMaxTotal;
          for(int alpha = 0; alpha < pTotal; alpha++)
          {
            g[M*w + j] += *(C_offset++)*f*targetCenterMonomials[alpha];
          }
        }                      
      }
    }
  }

  //printf( "target p histogram: ");
  //for( int i = 0; i < pMax; i++ )
  //  printf( " %i", pHistogram[i] );
  //printf( "\n" );
  //delete [] pHistogram;

  //////////////////////////////////////////////////////////////////////////////
  // Release memory
  //////////////////////////////////////////////////////////////////////////////
  delete [] rySquare;
  delete [] ry;
  delete [] C;
  delete [] dy;
  delete [] targetCenterMonomials;

  delete [] targetDists2;
  delete [] pMaxTotals;

  annDeallocPts(dataPts);
  delete [] nnIdx;              
  delete [] dists;
  delete kdTree;
  annClose();   
  return 0;
#endif
}

//------------------------------------------------------------------------------
// This function evaluates the IFGT by using a different truncation number for
//   each cluster to ensure that the total error contribution from each
//   cluster satisfies the error bound.  Thus, this is the cluster-wise adaptive
//   version of the IFGT. 
//
// This function uses a tree for finding nearby cluster centers.
//   
// See 'Automatic online tuning for fast Gaussian summation,' by Morariu et al, 
//   NIPS 2008 for details.  
//
// NOTES: 1) Currently is not implemented to work for W>1.  For W>1, use the
//           point-wise adaptive version instead.
//
//        2) The method could be extended to use different truncation numbers 
//           for each target by splitting targets into concentric regions and
//           computing the cluster-wise truncation for each concentric region
//           separately (because each region will have a different max distance 
//           from the cluster center, the truncations will differ by region).  
//           However, the current implementation uses only one region for the 
//           targets, and varies the truncation by cluster for that region.
//
// Created by Vlad Morariu on 2008-06-04.
//------------------------------------------------------------------------------
int figtreeEvaluateIfgtTreeAdaptiveCluster( int d, int N, int M, int W, double * x, 
                                        double h, double * q, double * y, 
                                        int pMax, int K, int * clusterIndex, 
                                        double * clusterCenter, double * clusterRadii,
                                        double r, double epsilon, int * clusterTruncations,
                                        double * g )
{
#ifdef FIGTREE_NO_ANN
  printf("This code was not compiled with support for ANN.  Please recompile\n");
  printf("with compiler flag 'FIGTREE_NO_ANN' not set to enable ANN support.\n");
  return -1;
#else
  // check input arguments
  FIGTREE_CHECK_POS_NONZERO_INT( d, figtreeEvaluateIfgtTreeAdaptiveCluster );
  FIGTREE_CHECK_POS_NONZERO_INT( N, figtreeEvaluateIfgtTreeAdaptiveCluster );
  FIGTREE_CHECK_POS_NONZERO_INT( M, figtreeEvaluateIfgtTreeAdaptiveCluster );
  FIGTREE_CHECK_POS_NONZERO_INT( W, figtreeEvaluateIfgtTreeAdaptiveCluster );
  FIGTREE_CHECK_NONNULL_PTR( x, figtreeEvaluateIfgtIfgtTreeAdaptiveCluster );
  FIGTREE_CHECK_POS_NONZERO_DOUBLE( h, figtreeEvaluateIfgtTreeAdaptiveCluster );
  FIGTREE_CHECK_NONNULL_PTR( g, figtreeEvaluateIfgtTreeAdaptiveCluster );
  FIGTREE_CHECK_NONNULL_PTR( y, figtreeEvaluateIfgtTreeAdaptiveCluster );
  FIGTREE_CHECK_POS_NONZERO_INT( pMax, figtreeEvaluateIfgtTreeAdaptiveCluster );
  FIGTREE_CHECK_POS_NONZERO_INT( K, figtreeEvaluateIfgtTreeAdaptiveCluster );
  FIGTREE_CHECK_NONNULL_PTR( clusterIndex, figtreeEvaluateIfgtTreeAdaptiveCluster );
  FIGTREE_CHECK_NONNULL_PTR( clusterCenter, figtreeEvaluateIfgtTreeAdaptiveCluster );
  FIGTREE_CHECK_NONNULL_PTR( clusterRadii, figtreeEvaluateIfgtTreeAdaptiveCluster );
  FIGTREE_CHECK_POS_NONZERO_DOUBLE( r, figtreeEvaluateIfgtTreeAdaptiveCluster );
  FIGTREE_CHECK_POS_NONZERO_DOUBLE( epsilon, figtreeEvaluateIfgtTreeAdaptiveCluster );
  FIGTREE_CHECK_NONNULL_PTR( g, figtreeEvaluateIfgtTreeAdaptiveCluster );

  //Memory allocation
  int pMaxTotal = nchoosek(pMax - 1 + d, d);
  int * pMaxTotals = new int[pMax];
  for( int i = 0; i < pMax; i++ )
    pMaxTotals[i] = nchoosek( i + d, d );

  double hSquare=h*h;
  double * targetCenterMonomials = new double[pMaxTotal];
  double * dy = new double[d];
  double * C = new double[W*K*pMaxTotal];
  double * ry = new double[K];
  double * rySquare = new double[K];

  double rx = clusterRadii[0];
  for(int i = 0; i < K; i++)
  {
    ry[i] = r + clusterRadii[i];
    rySquare[i] = ry[i]*ry[i];
    rx = MAX( rx, clusterRadii[i] );
  } 

  //
  // Build tree on cluster centers
  //
  double rSquare = (r+rx)*(r+rx);
  //Allocate storage using ANN procedures 
  ANNpointArray dataPts = annAllocPts(K,d);     // allocate data points
  ANNidxArray   nnIdx   = new ANNidx[K];        // allocate near neigh indices
  ANNdistArray  dists = new ANNdist[K];         // allocate near neighbor dists
  
  // Copy the cluster centers to the ANN data structure
  for (int k = 0; k < K; k++)
  {
    for ( int j = 0; j < d; j++ )
      dataPts[k][j]= clusterCenter[k*d + j];
  }

  // build search structure
  ANNkd_tree * kdTree = new ANNkd_tree(              
                                        dataPts,  // the data points
                                        K,        // number of points
                                        d,        // dimension of space
                                        1,
                                        ANN_KD_SUGGEST);

  //////////////////////////////////////////////////////////////////////////////
  // Evaluate  
  //////////////////////////////////////////////////////////////////////////////
  // for each cluster, compute max distances at which we can use a certain truncation number
  computeCAdaptiveCluster( d, N, W, K, pMaxTotal, pMax, h, clusterIndex, x, q, clusterCenter, clusterTruncations, pMaxTotals, C );

  memset( g, 0, sizeof(double)*M*W );
  
  for(int j = 0; j < M; j++)
  {
    int targetBase = j*d;
    ANNpoint queryPt=&(y[targetBase]);

    int NN = kdTree->annkFRSearchUnordered(           // search
                                   queryPt,  // query point
                                   rSquare,  // squared radius
                                   K,        // number of near neighbors
                                   nnIdx,     // nearest neighbors (returned)
                                   dists,     // distance (returned)
                                   0.0 );
    for(int l = 0; l < NN; l++)
    {
      int k = nnIdx[l];
      int centerBase = k*d;
      int p = clusterTruncations[k];
      int pTotal = pMaxTotals[p-1];

      double targetCenterDistanceSquare = dists[l];
      if(targetCenterDistanceSquare <= rySquare[k])
      {
        for(int i = 0; i < d; i++)
        {
          dy[i] = y[targetBase + i] - clusterCenter[centerBase + i];
        }
        computeTargetCenterMonomials( d, h, dy, p, targetCenterMonomials );
        double f=exp(-targetCenterDistanceSquare/hSquare);
        for(int w = 0; w < W; w++ )
        {        
          double * C_offset = C + (K*w + k)*pMaxTotal;
          for(int alpha = 0; alpha < pTotal; alpha++)
          {
            g[M*w + j] += *(C_offset++)*f*targetCenterMonomials[alpha];
          }
        }                      
      }
    }
  }

  //////////////////////////////////////////////////////////////////////////////
  // Release memory
  //////////////////////////////////////////////////////////////////////////////
  delete [] rySquare;
  delete [] ry;
  delete [] C;
  delete [] dy;
  delete [] targetCenterMonomials;
  delete [] pMaxTotals;

  annDeallocPts(dataPts);
  delete [] nnIdx;              
  delete [] dists;
  delete kdTree;
  annClose();  
  return 0;
#endif
}

int figtreeCalcMinMax( int d, int n, double * x, double * mins, double * maxs, int update )
{
  // check input arguments
  FIGTREE_CHECK_POS_NONZERO_INT( d, figtreeCalcMinMax );
  FIGTREE_CHECK_POS_NONZERO_INT( n, figtreeCalcMinMax );
  FIGTREE_CHECK_NONNULL_PTR( x, figtreeCalcMinMax );
  FIGTREE_CHECK_NONNULL_PTR( mins, figtreeCalcMinMax );
  FIGTREE_CHECK_NONNULL_PTR( maxs, figtreeCalcMinMax );

  // use first sample values as current min and max if we're not updating 
  //   some previously computed min and max values.  
  if( update != 1 && n > 0 )
  {
    for( int i = 0; i < d; i++ )
    {
      mins[i] = x[i];
      maxs[i] = x[i];
    }
  }

  // go through each sample in x and update mins and maxs for each dimension
  for( int i = 0; i < n; i++ )
  {
    for( int j = 0; j < d; j++ )
    {
      mins[j] = MIN( mins[j], x[i*d+j] );
      maxs[j] = MAX( maxs[j], x[i*d+j] );
    }
  }

  return 0;
}

int figtreeCalcMaxRange( double d, double * mins, double * maxs, double * maxRange )
{
  FIGTREE_CHECK_POS_NONZERO_INT( d, figtreeCalcMaxRange );
  FIGTREE_CHECK_NONNULL_PTR( mins, figtreeCalcMaxRange );
  FIGTREE_CHECK_NONNULL_PTR( maxs, figtreeCalcMaxRange );
  FIGTREE_CHECK_NONNULL_PTR( maxRange, figtreeCalcMaxRange );

  double maxRangeTemp = maxs[0] - mins[0];
  for( int i = 0; i < d; i++ )
    maxRangeTemp = MAX( maxRangeTemp, maxs[i] - mins[i] );
  *maxRange = maxRangeTemp;
  return 0;
}


///////////////////////////////////////////////////////////////////////////////////
//
//
//  Methods for automatic selection of evaluation method.
//  Functions created by Vlad Morariu 05-02-2008.
//  Functions modified by Vlad Morariu 06-06-2008.
//  Functions modified by Vlad Morariu 11-02-2008.  
//    - Added code to perform sub-sampling w/o replacement, and cleaned up code a 
//      little for release.
//
///////////////////////////////////////////////////////////////////////////////////

//------------------------------------------------------------------------------
// The functions below are used to estimate avg number of neighbors in the source set
//   given a query point from the target set.  
//
// Created by Vlad Morariu 2008-06-04.
//------------------------------------------------------------------------------
#ifndef FIGTREE_NO_ANN
inline void figtreeGetAverageNumNeighbors( ANNkd_tree * kdTree, int d, int M, double * y, double r, int Msample, double * avgNbrs, double * avgAnnFlops )
{
  int numNbrs = 0;
  int flopsAccum = 0, flops =0;
  double rSquare = r*r;
  double epsANN = 0.0;

  for(int i = 0; i < Msample; i++)
  {
    int idx = rand() % M;
    ANNpoint queryPt = &(y[ idx*d ]);

    int NN = kdTree->annkFRSearchUnorderedFlops( // fixed radius search
                                   queryPt,  // query point
                                   rSquare,  // squared radius
                                   0,        // number of near neighbors
                                   NULL,     // nearest neighbors (returned)
                                   NULL,     // distance (returned)
                                   epsANN,
                                   &flops );
    numNbrs += NN;
    flopsAccum += flops;
  }

  *avgNbrs = numNbrs / (double)Msample;
  *avgAnnFlops = flopsAccum / (double)Msample;
}

inline void figtreeEstimatedNeighborSources( int d, int M, double * y, double h, double epsilon, ANNkd_tree * sourcesKdTree, int Msample, double * avgNbrSources, double * avgAnnFlopsSources )
{
  // Compute the cutoff radius
  double r = h*sqrt(log(1/epsilon));

  // estimate avg number of neighbors
  figtreeGetAverageNumNeighbors( sourcesKdTree, d, M, y, r, Msample, avgNbrSources, avgAnnFlopsSources ); 
}

inline void figtreeEstimatedNeighborClusters( int d, int M, double * y, int K, double * clusterRadii, double r, ANNkd_tree * clustersKdTree, int Msample, double * avgNbrClusters, double * avgAnnFlopsClusters )
{
  //Find the maximum cluster radius
  double pcrMax = clusterRadii[0];
  for(int i = 0; i < K; i++)
  {
    if (clusterRadii[i] > pcrMax)
    {
      pcrMax = clusterRadii[i];
    }
  }
  double rMax =(r+pcrMax);

  // estimate numbers of pts with more than 1 neighbor and avg number of neighbors
  figtreeGetAverageNumNeighbors( clustersKdTree, d, M, y, rMax, Msample, avgNbrClusters, avgAnnFlopsClusters ); 
}
#endif
inline void figtreeEstimatedNeighborClustersNoAnn( int d, int N, int M, double h, double * y, 
                                                   int K, double * clusterCenter, double * clusterRadii, double r, int Msample,
                                                   double * avgNbrClustersNoAnn, double * avgFindCentersFlops )
{
  double flops = 0;
  double avgNbrClusters = 0;
  double * dy = new double[d];
  double * ry = new double[K];
  double * rySquare = new double[K];

  for(int i = 0; i < K; i++)
  {
    ry[i] = r + clusterRadii[i];
    rySquare[i] = ry[i]*ry[i];
  }   

  for(int j = 0; j < Msample; j++)
  {
    int targetBase = (rand()%M)*d;        
    for(int k = 0; k < K; k++)
    {
      int centerBase = k*d;
      double targetCenterDistanceSquare = 0.0;
      for(int i = 0; i < d; i++)
      {
        dy[i] = y[targetBase + i] - clusterCenter[centerBase + i];
        targetCenterDistanceSquare += dy[i]*dy[i];
        flops += 3;
        if(targetCenterDistanceSquare > rySquare[k]) break;
      }

      if( targetCenterDistanceSquare <= rySquare[k] )
        avgNbrClusters++;
    }
  }

  *avgFindCentersFlops = flops/Msample;
  *avgNbrClustersNoAnn = avgNbrClusters/Msample;

  delete [] dy;
  delete [] ry;
  delete [] rySquare;
}


//------------------------------------------------------------------------------
//
// The functions below estimate the number of floating point operations (flops) 
//   for different parts of the figtree code.
// Where possible the estimates are made directly by counting the number of 
//   floating point operations in the code itself, but for things such as
//   building the kd-Tree and clustering, we use the theoretical complexity.
//
// Created by Vlad Morariu 2008-05-03 to 2008-06-04.
// Modified by Vlad Morariu 2008-12-05 
//   Changed floating op estimation functions to reflect revised versions of code.
//------------------------------------------------------------------------------
inline double figtreeEstimatedFlopsComputeSourceCenterMonomials( int d, int pMaxTotal )
{
  return (d + (double)pMaxTotal);
}

inline double figtreeEstimatedFlopsComputeTargetCenterMonomials( int d, int pMaxTotal )
{
  return (d + (double)pMaxTotal);
}

inline double figtreeEstimatedFlopsComputeConstantSeries( int pMaxTotal )
{
  return (2.0*pMaxTotal);
}

inline double figtreeEstimatedFlopsComputeC( int d, int N, int W, int K, int pMaxTotal, double flopsExp )
{
  double computeSourceCenterMonomials = figtreeEstimatedFlopsComputeSourceCenterMonomials( d, pMaxTotal );
  double computeConstantSeries = figtreeEstimatedFlopsComputeConstantSeries( pMaxTotal );
  return N*(3*d + computeSourceCenterMonomials + W*(2+flopsExp+2*pMaxTotal)) + computeConstantSeries + W*K*pMaxTotal;
}

inline double figtreeEstimatedFlopsBuildTree( int d, int N )
{
  return d*N*log((double)N);
}

inline double figtreeEstimatedFlopsKCenterClustering( int d, int N, int K )
{
  return 3*d*(N*log((double)K) + K*(double)K); // not using the NlogK version yet
}

inline double figtreeEstimatedFlopsDirect( int d, int N, int M, int W, double flopsExp )
{
  return W*(1 + ((double)M)*((double)N)*(d*3 + 3 + flopsExp));
}

inline double figtreeEstimatedFlopsDirectTree( int d, int N, int M, int W, double avgNbrSources, double avgAnnFlopsSources, double flopsExp )
{
  return W*(7 + ((double)M)*( avgAnnFlopsSources + avgNbrSources*(3+flopsExp) ));
}

inline double figtreeEstimatedFlopsIfgt( int d, int N, int M, int W, int K, int pMaxTotal, double avgNbrClustersNoAnn, double avgFindCentersFlops, double flopsExp )
{
  // cost of computing coeffs from sources
  double computeC = figtreeEstimatedFlopsComputeC( d, N, W, K, pMaxTotal, flopsExp );
  double source = 2*K + 1 + computeC;

  // cost of computing target monomials and evaluating
  double computeTargetCenterMonomials = figtreeEstimatedFlopsComputeTargetCenterMonomials( d, pMaxTotal );
  double target = ((double)M)*(avgFindCentersFlops + avgNbrClustersNoAnn*(computeTargetCenterMonomials + 1 + flopsExp + W*3*((double)pMaxTotal)));

  return source + target;
}

inline double figtreeEstimatedFlopsIfgtTree( int d, int N, int M, int W, int K, int pMaxTotal, double avgNbrClusters, double avgAnnFlopsClusters, double flopsExp )
{
  // cost of computing coeffs from sources
  double computeC = figtreeEstimatedFlopsComputeC( d, N, W, K, pMaxTotal, flopsExp );
  double source = 4 + computeC;

  // cost of computing target monomials and evaluating
  double computeTargetCenterMonomials = figtreeEstimatedFlopsComputeTargetCenterMonomials( d, pMaxTotal );
  double target =  ((double)M)*( avgAnnFlopsClusters + avgNbrClusters*(d + computeTargetCenterMonomials + 1 + flopsExp + W*3*((double)pMaxTotal) ) );

  return source + target;
}

//------------------------------------------------------------------------------
// This function chooses the evaluation method for figtree, given the input
// parameters and data.
//
// Created by Vlad Morariu 2008-05-03 to 2008-06-04.
//------------------------------------------------------------------------------
int figtreeChooseEvaluationMethod( int d, int N, int M, int W, double * x, double h,
                                   double * y, double epsilon, int ifgtParamMethod, int verbose,
                                   int * bestMethod, double * flops, void * data_struct )
{
  int ret = 0;           // return value (0 is no error, -1 is error)

  //
  // Parameters for estimating avg number of neighbors, flops needed to find neighbors using tree,
  // and setting k limit.
  //
  // In future releases of the code, users should be able to choose these parameters as they affect the
  // quality of method selection.  For now, users will have to recompileif they want to change these 
  // parameters.
  //
  int Msample = M_SAMPLE; // how many of the target points do we sample to estimate avg number of neighbors and flops?
  int Nss = MAX( MIN(N,N_SS_MIN), (int)pow(N,N_SS_POW) ); // number of subsampled sources
  double kLimitToAvgNbrRatio = K_LIMIT_TO_AVG_NBR_RATIO;  // set k limit to this times avg number of source neighbors
  double flopsExp = FLOPS_EXP;                            // floating point ops per exp() call

  // other values computed fromthe parameters above
  double ss = N/(double)Nss;       // how much subsampling to do for building kd-Tree
  int kLimit = N, kMax;            // kLimit is the highest K that we allow, kMax is the max K-value chosen by param selection

  // initialize flop estimates; -1 indicates flop estimation was not completed
  double flopsDirect = figtreeEstimatedFlopsDirect( d, N, M, W, flopsExp );
  double flopsDirectTree = -1;
  double flopsIfgt = -1;
  double flopsIfgtTree = -1;

  // avg number of neighboring clusters for IFGT (no tree), and avg flops spent finding
  //   the neighboring clusters
  double avgNbrClustersNoAnn, avgFindCentersFlops;

  // this datastructure will hold the kd-Trees, cluster centers, etc
  FigtreeData data = figtreeCreateData();  // obtain data structure filled with 0's

#ifndef FIGTREE_NO_ANN

  // these variables will represent statistics of source distribution (avg neighbors)
  double avgNbrSources, avgNbrClusters;

  // these variables contain avg cost of each kd-Tree query for finding neighbors
  double avgAnnFlopsSources, avgAnnFlopsClusters;

  // get a random sampling of Nss points, w/o replacement
  int * shuffled_indexes = new int[N];
  for( int i = 0; i < N; i++ )
    shuffled_indexes[i] = i;
  std::random_shuffle( shuffled_indexes, shuffled_indexes+N );  // could also use random_sample, if available

  // Allocate storage using ANN procedures 
  data.annSources = annAllocPts(Nss,d);  // allocate data points

  // Copy the source points to the ANN data structure
  for (int i = 0; i < Nss; i++)
  {
    for ( int j = 0; j < d; j++ )
      data.annSources[i][j]= x[shuffled_indexes[i]*d + j];
  }

  delete [] shuffled_indexes;

  // Build search structure
  data.annSourcesKdTree = new ANNkd_tree(              
                                        data.annSources, // the data points
                                        Nss,             // number of points
                                        d,               // dimension of space
                                        1,
                                        ANN_KD_SUGGEST );

  figtreeEstimatedNeighborSources( d, M, y, h, epsilon, data.annSourcesKdTree, Msample, &avgNbrSources, &avgAnnFlopsSources);
  flopsDirectTree  = figtreeEstimatedFlopsBuildTree( d, N ) + figtreeEstimatedFlopsDirectTree( d, N, M, W, ss*avgNbrSources, ss*avgAnnFlopsSources*(log((double)N)/log((double)Nss)), flopsExp );

  // set kLimit based on how many source points, on average, are within radius of each query pt
  kLimit = MIN((int)(kLimitToAvgNbrRatio*ss*avgNbrSources),N);

#endif // #ifndef FIGTREE_NO_ANN

  if( kLimit > 0 )
  {
    double maxRange=0;

    // calculate R, the diameter of the hypercube that encompasses sources and targets
    double * mins = new double[d];
    double * maxs = new double[d];
    figtreeCalcMinMax( d, N, x, mins, maxs, 0 );
    figtreeCalcMinMax( d, M, y, mins, maxs, 1 );
    figtreeCalcMaxRange( d, mins, maxs, &maxRange );
    delete [] mins;
    delete [] maxs;

    if( ifgtParamMethod == FIGTREE_PARAM_NON_UNIFORM )
      ret = figtreeChooseParametersNonUniform( d, N, x, h, epsilon, kLimit, maxRange, &kMax, &data.pMax, &data.r, NULL );
    else
      ret = figtreeChooseParametersUniform( d, h, epsilon, kLimit, maxRange, &kMax, &data.pMax, &data.r, NULL );
    if( ret < 0 )
    {
      printf("figtree: figtreeChooseParameters%sUniform() failed.\n", 
             ((ifgtParamMethod == FIGTREE_PARAM_NON_UNIFORM) ? "Non" : ""));
      return ret;
    }

    verbose && printf("figtreeChooseParameters%sUniform() chose p=%i, k=%i.\n", 
                      ((ifgtParamMethod == FIGTREE_PARAM_NON_UNIFORM) ? "Non" : ""), data.pMax, kMax );

    // do k-center clustering
    data.clusterIndex = new int[N];
    data.numPoints    = new int[kMax]; 
    data.clusterCenters = new double[d*kMax];
    data.clusterRadii = new double[kMax];

    ret = figtreeKCenterClustering( d, N, x, kMax, &data.K, &data.rx, data.clusterIndex, 
                                 data.clusterCenters, data.numPoints, data.clusterRadii );
    if( ret >= 0 )
    {
      double errorBound = epsilon + 1;
      // choose truncation number again now that clustering is done
      ret = figtreeChooseTruncationNumber( d, h, epsilon, data.rx, maxRange, &data.pMax, &errorBound );
      if( ret >= 0 )
      {
        // these are the params we would evaluate IFGT with
        verbose && printf( "Eval IFGT(h= %3.2e, pMax= %i, K= %i, r= %3.2e, rx= %3.2e, epsilon= %3.2e, bound = %3.2e)\n", 
                           h, data.pMax, data.K, data.r, data.rx, epsilon, errorBound);

        double pMaxTotalDouble = nchoosek_double(data.pMax-1+d,d);
        if( errorBound <= epsilon && pMaxTotalDouble < INT_MAX )
        {
          data.pMaxTotal = (int)pMaxTotalDouble;

          //
          // Estimate number of flops for performing IFGT
          //
          figtreeEstimatedNeighborClustersNoAnn( d, N, M, h, y, data.K, data.clusterCenters, data.clusterRadii, data.r, Msample, &avgNbrClustersNoAnn, &avgFindCentersFlops );
          flopsIfgt = figtreeEstimatedFlopsIfgt( d, N, M, W, data.K, data.pMaxTotal, avgNbrClustersNoAnn, avgFindCentersFlops, flopsExp );

#ifndef FIGTREE_NO_ANN
          // 
          // Estimate number of flops for performing IFGT with kd-Tree
          //

          // Allocate storage using ANN procedures 
          data.annClusters = annAllocPts(data.K,d);  // allocate data points

          // Copy the source points to the ANN data structure
          for (int i = 0; i < data.K; i++)
          {
            for ( int j = 0; j < d; j++ )
              data.annClusters[i][j]= data.clusterCenters[i*d + j];
          }

          // build search structure
          data.annClustersKdTree = new ANNkd_tree(              
                                                data.annClusters, // the data points
                                                data.K,           // number of points
                                                d,                // dimension of space
                                                1,
                                                ANN_KD_SUGGEST );
        
          figtreeEstimatedNeighborClusters( d, M, y, data.K, data.clusterRadii, data.r, data.annClustersKdTree, Msample, &avgNbrClusters, &avgAnnFlopsClusters );     
          flopsIfgtTree    = figtreeEstimatedFlopsIfgtTree( d, N, M, W, data.K, data.pMaxTotal, avgNbrClusters, avgAnnFlopsClusters, flopsExp );
#endif
        }
      }
      else // if figtreeChooseTruncationNumber fails...
      {
        printf("figtree: figtreeChooseTruncationNumber() failed.\n");
      }
    } // if figtreeKCenterClustering fails...
    else
    {
      printf("figtree: figtreeKCenterClustering() failed.\n");
    }
  }

  // save clustering and parameter selection data if desired
  if( data_struct != NULL )
  {
    FigtreeData * data_out = (FigtreeData*)data_struct;
    // copy all data except for ANN data structures
    *data_out = data;
    data_out->annClusters = NULL;
    data_out->annClustersKdTree = NULL;
    data_out->annSources = NULL;
    data_out->annSourcesKdTree = NULL;

    // set copied data to NULL, or else it will be deallocated
    data.clusterCenters = NULL;
    data.clusterIndex = NULL;
    data.clusterRadii = NULL;
    data.numPoints = NULL;
  }

  figtreeReleaseData( &data );

  // return flops per method if user supplies non-NULL pointer
  if( flops != NULL )
  {
    flops[ FIGTREE_EVAL_DIRECT      ] = flopsDirect;
    flops[ FIGTREE_EVAL_DIRECT_TREE ] = flopsDirectTree;
    flops[ FIGTREE_EVAL_IFGT        ] = flopsIfgt;
    flops[ FIGTREE_EVAL_IFGT_TREE   ] = flopsIfgtTree;
  }

  // choose best method if user supplies non-NULL pointer
  if( bestMethod!= NULL )
  {
    double bestFlops = flopsDirect;
    *bestMethod = FIGTREE_EVAL_DIRECT;
    
    if( flopsDirectTree != -1 && flopsDirectTree < bestFlops )
    {
      *bestMethod = FIGTREE_EVAL_DIRECT_TREE;
      bestFlops = flopsDirectTree;
    }

    if( flopsIfgt!= -1 && flopsIfgt < bestFlops )
    {
      *bestMethod = FIGTREE_EVAL_IFGT;
      bestFlops = flopsIfgt;
    }

    if( flopsIfgtTree != -1 && flopsIfgtTree < bestFlops )
    {
      *bestMethod = FIGTREE_EVAL_IFGT_TREE;
      bestFlops = flopsIfgtTree;
    }
  }

#ifndef FIGTREE_NO_ANN
  annClose();
#endif

  return ret;
}
