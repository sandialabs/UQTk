// File: figtree.h
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
//   Rename library to FIGTree (and some other function renamimg)
//
// Modified: 05-03-08 by Vlad Morariu 
//   Add initial version of method selection code.
//
// Modified: 05-27-08 by Vlad Morariu 
//   Change figtreeChooseParameters* and figtreeChooseTruncationNumber functions
//   to return the predicted errorBound.  This can then be used to check if
//   the parameters chosen will satisfy the desired error bound (they may not 
//   since we enforce a limit on pMax (the truncation number).
//
// Modified: 05-29-08 to 06-10-08 by Vlad Morariu
//   Make changes for 'Automatic online tuning for fast Gaussian summation,' by
//   Morariu et al, NIPS 2008.  Changes include automatic parameter selection 
//   improvements, point-wise and cluster-wise truncation number selection, and
//   automatic method selection.  These changes make this approach much easier
//   to use since the user does not need to choose parameters.
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

//------------------------------------------------------------------------------
// The code was written by Vlad Morariu, Vikas Raykar, and Changjiang Yang 
// and is copyrighted under the Lesser GPL: 
//
// Copyright (C) 2008 Vlad Morariu, Vikas Raykar, and Changjiang Yang 
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
#ifndef FIGTREE_H
#define FIGTREE_H

#ifdef WIN32
  //----------------------------------------------------------------------
  //  To compile the code into a windows DLL, you must define the 
  //  symbol FIGTREE_DLL_EXPORTS. 
  // 
  //  To compile the code statically into a windows executable 
  //    (i.e. not using a separate DLL) define FIGTREE_DLL_STATIC.
  //----------------------------------------------------------------------
  #ifdef FIGTREE_STATIC
    #define FIGTREE_DLL_API  // since FIGTREE_STATIC is defined, code is statically 
                          // linked and no exports or imports are needed
  #else
    #ifdef FIGTREE_DLL_EXPORTS
      #define FIGTREE_DLL_API __declspec(dllexport)
    #else
      #define FIGTREE_DLL_API __declspec(dllimport)
    #endif
  #endif
#else
  //----------------------------------------------------------------------
  // FIGTREE_DLL_API is ignored for all other systems
  //----------------------------------------------------------------------
  #define FIGTREE_DLL_API
#endif

// we want to export C functions if using DLL
#if defined(__cplusplus) && !defined(FIGTREE_STATIC)
extern "C" {
#endif




//------------------------------------------------------------------------------
//
// Some useful constants for the figtree() function call
//
//------------------------------------------------------------------------------

// Constants for choosing evaluation method
#define FIGTREE_EVAL_DIRECT        0  // direct evaluation of gauss transform
#define FIGTREE_EVAL_IFGT          1  // truncated taylor series evaluation
#define FIGTREE_EVAL_DIRECT_TREE   2  // direct evaluation
#define FIGTREE_EVAL_IFGT_TREE     3  // ifgt+tree
#define FIGTREE_EVAL_AUTO          4  // automatically chooses one of the four
#define FIGTREE_EVAL_SIZE          5  // total number of eval methods

// Constants for choosing parameter selection method for IFGT
#define FIGTREE_PARAM_UNIFORM      0  // estimate params assuming sources are 
                                      //   uniformly distributed (not recommended 
                                      //   for non-uniform data)
#define FIGTREE_PARAM_NON_UNIFORM  1  // estimate params by using actual source 
                                      //   distribution (runs k-center clustering 
                                      //   twice, but speedup during evaluation
                                      //   more than makes up for it in general)
#define FIGTREE_PARAM_SIZE         2  // total number of param estimation methods

// Constants for choosing truncation selection method for IFGT
#define FIGTREE_TRUNC_MAX          0  // use worst case truncation number for all pts
#define FIGTREE_TRUNC_POINT        1  // use point-wise error bounds for individual truncations
#define FIGTREE_TRUNC_CLUSTER      2  // use cluster-wise error bounds for individual truncations
#define FIGTREE_TRUNC_SIZE         3  // total number of truncation selection methods

//------------------------------------------------------------------------------
//
// Main functions exported by the library.  For more control over internal
//   workings of figtree, see figtree_internal.h and figtree.cpp.
//
//------------------------------------------------------------------------------


// Note: All matrix pointers are assumed to point to a contiguous one 
//   dimensional array containing the entries of the matrx in row major format.
//   Thus, for an M x N matrix and a pointer ptr to its data, 
//   ptr[0] ... ptr[N-1] contains the first row of the matrix, and so on. 

// Evaluates gauss transform in one shot (chooses parameters and does clustering
// if necessary, and evaluates).
//
// Input
//    * d --> data dimensionality.
//    * N --> number of source points.
//    * M --> number of target points.
//    * W --> number of weights that will be used for each source point. 
//            This is useful if one needs multiple gauss transforms that have
//            the same sources, targets, and bandwidth, but different
//            weights/strengths (q). By computing coefficients for all W weight
//            sets at once, we avoid duplicating much of the overhead.  However,
//            more memory is needed to store a set of coefficients for each set
//            of weights.
//    * x --> N x d matrix of N source points in d dimensions.
//    * h --> the source scale or bandwidth.
//    * q --> W x N matrix of the source strengths.
//    * y --> M x d matrix of M target points in d dimensions.
//    * epsilon --> desired error
//    * evalMethod --> the evaluation method to use in evaluating gauss 
//            transform. Can be FIGTREE_EVAL_[DIRECT,IFGT,DIRECT_TREE,
//            IFGT_TREE], defined above. epsilon is needed for all but 
//            DIRECT method.  Parameter selection is done only in the IFGT
//            or IFGT_TREE case.  
//    * ifgtParamMethod --> the method to use for determining parameters.
//            Can be FIGTREE_PARAM_UNIFORM or FIGTREE_PARAM_NON_UNIFORM. 
//    * ifgtTruncMethod --> the method to use for determining where to truncate
//            Taylor series for each source and target point.  Can be
//            FIGTREE_TRUNC_MAX, FIGTREE_TRUNC_POINT, and FIGTREE_TRUNC_CLUSTER
//    * verbose --> if nonzero, prints parameters chosen for evaluation
//
// Output
//    * g --> W x M vector of the Gauss Transform evaluated at each target point. 
//            The ith row is the result of the transform using the ith set of 
//            weights.
FIGTREE_DLL_API 
int figtree( int d, int N, int M, int W, double * x, double h, 
             double * q, double * y, double epsilon, double * g,
             int evalMethod = FIGTREE_EVAL_AUTO,
             int ifgtParamMethod = FIGTREE_PARAM_NON_UNIFORM, 
             int ifgtTruncMethod = FIGTREE_TRUNC_CLUSTER,
             int verbose = 0 );

// Chooses between evaluation methods
// Input
//    * d --> data dimensionality.
//    * N --> number of source points.
//    * M --> number of target points.
//    * W --> number of weights that will be used for each source point. 
//            This is useful if one needs multiple gauss transforms that have
//            the same sources, targets, and bandwidth, but different
//            weights/strengths (q). By computing coefficients for all W weight
//            sets at once, we avoid duplicating much of the overhead.  However,
//            more memory is needed to store a set of coefficients for each set
//            of weights.
//    * x --> N x d matrix of N source points in d dimensions.
//    * h --> the source scale or bandwidth.
//    * y --> M x d matrix of M target points in d dimensions.
//    * epsilon --> desired error
//    * paramMethod --> the method to use for determining parameters.
//            Can be FIGTREE_PARAM_UNIFORM or FIGTREE_PARAM_NON_UNIFORM.
//    * verbose --> if nonzero, prints parameters chosen for evaluation
//
// Output
//    * bestEvalMethod --> if non-NULL, the evaluation method to use in evaluating gauss 
//            transform. Can be FIGTREE_EVAL_[DIRECT,IFGT,DIRECT_TREE,
//            IFGT_TREE], defined above. 
//
//    * flops --> if non-NULL, a double array of size FIGTREE_EVAL_SIZE, indexed by eval method
//                type where flops[evalMethod] = estimated number of flops if we use this method.
//                For evalMethod=FIGTREE_EVAL_IFGT or evalMethod=FIGTREE_EVAL_IFGT_TREE,
//                it is possible that flops[.] = -1 in the case that it is decided early (before
//                finishing the estimation) that these two methods will be slower than
//                FIGTREE_EVAL_DIRECT or FIGTREE_EVAL_DIRECT_TREE
//    * data_struct --> a structure that keeps some structures (k-center clustering, other params)
//                that were computed during method selection that can be reused 
//
FIGTREE_DLL_API
int figtreeChooseEvaluationMethod( int d, int N, int M, int W, double * x, double h,
                                 double * y, double epsilon, int paramMethod=FIGTREE_PARAM_NON_UNIFORM, int verbose=0,
                                 int * bestEvalMethod=0, double * flops=0, void * data_struct=0 );


// Gonzalez's farthest-point clustering algorithm.
//
// Input
//    * d --> data dimensionality.
//    * N --> number of source points.
//    * x --> N x d matrix of N source points in d dimensions 
//       (in one contiguous array, row major format where each row is a point).
//    * kMax --> maximum number of clusters.
//
// Output
//    * K --> actual number of clusters (less than kMax if duplicate pts exist)
//    * rx --> maximum radius of the clusters (rx).
//    * clusterIndex --> vector of length N where the i th element is the 
//                cluster number to which the i th point belongs. 
//                ClusterIndex[i] varies between 0 to K-1.
//    * clusterCenters --> K x d matrix of K cluster centers 
//                (contiguous 1-d array, row major format).
//    * numPoints --> number of points in each cluster.
//    * clusterRadii --> radius of each cluster.
FIGTREE_DLL_API 
int figtreeKCenterClustering( int d, int N, double * x, int kMax, int * K,
                              double * rx, int * clusterIndex, double * clusterCenters, 
                              int * numPoints, double * clusterRadii );

#if defined(__cplusplus) && !defined(FIGTREE_STATIC)
} // extern "C"
#endif

#endif //FIGTREE_H