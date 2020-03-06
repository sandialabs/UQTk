#ifndef FIGTREE_INTERNAL_H
#define FIGTREE_INTERNAL_H

#include "figtree.h"

// we want to export C functions if using DLL
#if defined(__cplusplus) && !defined(FIGTREE_STATIC)
extern "C" {
#endif

// Given the maximum cluster radius, this function computes the maximum 
// truncation number that guarantees results within the desired error bound.
//Input
//    * d --> dimension of the points.
//    * h --> the source bandwidth.
//    * epsilon --> the desired error.
//    * rx --> maximum cluster radius
//    * maxRange --> max dimension range.  The range along a dimension is 
//          the difference between the max and min values that can ever
//          occur along that dimension.  The max dimension range is the 
//          maximum range among all dimensions.  For example, if all 
//          points lie in unit hypercube, then maxRange = 1.
//
//Output
//    * pMax --> maximum truncation number for the Taylor series.
//    * errorBound --> the error bound (if desired error bound is not reached, then it will be epsilon+1)
FIGTREE_DLL_API 
int figtreeChooseTruncationNumber( int d, double h, double epsilon, 
                                   double rx, double maxRange, int * pMax, double * errorBound );

// Chooses parameters for IFGT and FIGTree by assuming that sources are
//   uniformly distributed in a unit cube.
//Input
//    * d --> dimension of the points.
//    * h --> the source bandwidth.
//    * epsilon --> the desired error.
//    * kLimit --> upper limit on the number of clusters, kLimit.
//    * maxRange --> max dimension range.  The range along a dimension is 
//          the difference between the max and min values that can ever
//          occur along that dimension.  The max dimension range is the 
//          maximum range among all dimensions.  For example, if all 
//          points lie in unit hypercube, then maxRange = 1.
//
//Note : [ Use roughly kLimit=round(40*sqrt(d)/h) ]
//Output
//    * K --> number of clusters.
//    * pMax --> maximum truncation number for the Taylor series.
//    * r --> source cutoff radius.
//    * errorBound --> the expected error bound (if desired error bound is not reached, then it will be epsilon+1)
FIGTREE_DLL_API 
int figtreeChooseParametersUniform( int d, double h, double epsilon, 
                                    int kLimit, double maxRange,
                                    int * K, int * pMax, double * r, double * errorBound );

// Chooses parameters for IFGT and FIGTree without assumption that sources are
// uniformly distributed.  In this case, k-center clustering is done as part
// of the parameter selection so that the radius of each cluster is not 
// estimated but computed directly for the sources.  This results in an 
// additional k-center clustering operation, but the speedup when sources are
// not uniformly distributed can be large because
// optimal parameters are much more accurately estimated.  Even when sources are
// uniformly distributed, the slowdown is often small compared to 
// figtreeChooseParametersUniform.
// 
// Input
//    * d --> dimension of the points.
//    * N --> number of source points.
//    * x --> N x d matrix of N source points in d dimensions.
//    * h --> the source bandwidth.
//    * epsilon --> the desired error.
//    * kLimit --> upper limit on the number of clusters, K.
// Note : Use kLimit=N to allow for optimal param estimation.
//    * maxRange --> max dimension range.  The range along a dimension is 
//          the difference between the max and min values that can ever
//          occur along that dimension.  The max dimension range is the 
//          maximum range among all dimensions.  For example, if all 
//          points lie in unit hypercube, then maxRange = 1.
// Output
//    * K --> number of clusters.
//    * pMax --> maximum truncation number for the Taylor series.
//    * r --> source cutoff radius.
//    * errorBound --> the expected error bound (if desired error bound is not reached, then it will be epsilon+1)
FIGTREE_DLL_API 
int figtreeChooseParametersNonUniform( int d, int N, double * x, 
                                       double h, double epsilon, int kLimit, double maxRange,
                                       int * K, int * pMax, double * r, double * errorBound );

// Computes exact gauss transform (within machine precision) using direct 
// evaluation. Provided for time/error comparison.
//
// Input
//    * d --> data dimensionality.
//    * N --> number of source points.
//    * M --> number of target points.
//    * x --> N x d matrix of N source points in d dimensions.
//    * h --> the source scale or bandwidth.
//    * q --> 1 x N or N x 1 vector of the source strengths.
//    * y --> M x d matrix of M target points in d dimensions.
//
// Output
//    * g --> 1 x M vector of the Gauss Transform evaluated at each target
//            point. 
FIGTREE_DLL_API 
int figtreeEvaluateDirect( int d, int N, int M, double * x, double h, 
                           double * q, double * y, double * g );

// Computes an approximation to Gauss Transform.  Implementation based on:
// Fast computation of sums of Gaussians in high dimensions. Vikas C. Raykar, 
//   C. Yang, R. Duraiswami, and N. Gumerov, CS-TR-4767, Department of computer
//   science, University of Maryland, College Park.
//
// Input
//    * d --> data dimensionality.
//    * N --> number of source points.
//    * M --> number of target points.
//    * W --> number of weights that will be used for each source point. 
//            This really does multiple transforms, with different weights each
//            time but with same sources and targets.  This saves a lot of time
//            since most of the work is not duplicated.  However, it requires 
//            more memory to store the coefficients for each set of weights.
//    * x --> N x d matrix of N source points in d dimensions.
//    * h --> the source scale or bandwidth.
//    * q --> W x N vector of the source strengths.
//    * y --> M x d matrix of M target points in d dimensions.
//    * pMax --> maximum truncation number for the Taylor series.
//    * K --> the number of clusters.
//    * clusterIndex --> N x 1 vector the i th element is the cluster number 
//            to which the i th point belongs. [ ClusterIndex[i] varies between
//            0 to K-1. ]
//    * clusterCenter --> K x d matrix of K cluster centers.
//    * clusterRadii  --> K x 1 matrix of the radius of each cluster.
//    * r --> cutoff radius
//    * epsilon --> desired error
//
//Output
//    * g --> W x M vector of the Gauss Transform evaluated at each target
//            point. Each row q is the result of the transform using the qth set
//            of weights.
FIGTREE_DLL_API 
int figtreeEvaluateIfgt( int d, int N, int M, int W, double * x, 
                         double h, double * q, double * y, 
                         int pMax, int K, int * clusterIndex, 
                         double * clusterCenter, double * clusterRadii,
                         double r, double epsilon, double * g );

// Computes an approximation to Gauss Transform using approximate 
// nearest-neighbors. Same as figtreeEvaluateIfgt() but uses Approximate 
// Nearest-Neighbors library to find source clusters which influence each target
// (part of FIGTree). Parameters for this function can be computed using 
// figtreeChooseParameters[Non]Uniform and figtreeKCenterClustering. 
//
// Input
//    * d --> data dimensionality.
//    * N --> number of source points.
//    * M --> number of target points.
//    * W --> number of weights that will be used for each source point.
//            This really does multiple transforms, with different weights each
//            time but with same sources and targets.  This saves a lot of time
//            since most of the work is not duplicated.  However, it requires
//            more memory to store the coefficients for each set of weights.
//    * x --> N x d matrix of N source points in d dimensions.
//    * h --> the source scale or bandwidth.
//    * q --> W x N vector of the source strengths.
//    * y --> M x d matrix of M target points in d dimensions.
//    * pMax --> maximum truncation number for the Taylor series.
//    * K --> the number of clusters.
//    * clusterIndex --> N x 1 vector the i th element is the cluster number 
//            to which the i th point belongs. [ ClusterIndex[i] varies between
//            0 to K-1. ]
//    * clusterCenter --> K x d matrix of K cluster centers.
//    * clusterRadii  --> K x 1 matrix of the radius of each cluster.
//    * r --> cutoff radius
//    * epsilon --> desired error
//
// Output
//    * g --> W x M vector of the Gauss Transform evaluated at each target
//            point.  Each row of g is the result of the transform using one set
//            of weights.
FIGTREE_DLL_API 
int figtreeEvaluateIfgtTree( int d, int N, int M, int W, double * x, 
                             double h, double * q, double * y, 
                             int pMax, int K, int * clusterIndex, 
                             double * clusterCenter, double * clusterRadii,
                             double r, double epsilon, double * g );

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
FIGTREE_DLL_API
int figtreeEvaluateIfgtAdaptivePoint( int d, int N, int M, int W, double * x, 
                                      double h, double * q, double * y, 
                                      int pMax, int K, int * clusterIndex, 
                                      double * clusterCenter, double * clusterRadii,
                                      double r, double epsilon,
                                      double * g );

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
FIGTREE_DLL_API
int figtreeEvaluateIfgtAdaptiveCluster( int d, int N, int M, int W, double * x, 
                                        double h, double * q, double * y, 
                                        int pMax, int K, int * clusterIndex, 
                                        double * clusterCenter, double * clusterRadii,
                                        double r, double epsilon, int * clusterTruncations,
                                        double * g );

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
                                          double * g );

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
                                            double * g );

// Computes an approximation to Gauss Transform using approximate 
// nearest-neighbors.  Direct method (no taylor expansion is done), with tree 
// directly on samples (part of FIGTree). Requires Approximate 
// Nearest-Neighbor(ANN) library.
//
// Input
//    * d --> data dimensionality.
//    * N --> number of source points.
//    * M --> number of target points.
//    * x --> N x d matrix of N source points in d dimensions.
//    * h --> the source scale or bandwidth.
//    * q --> 1 x N vector of the source strengths.
//    * y --> M x d matrix of M target points in d dimensions.
//    * epsilon --> desired error
//
// Output
//    * g --> 1 x M vector of the Gauss Transform evaluated at each target
//            point.
FIGTREE_DLL_API 
int figtreeEvaluateDirectTree( int d, int N, int M, double * x, double h, 
                               double * q, double * y, double epsilon, double * g );

// Computes min and max values along each dimension.  Used to determine
//   the size of the hypercube (and the max distance R that any two pts 
//   can be from each other).
//
// Input
//    * d --> data dimensionality.
//    * n --> number of source points.
//    * x --> n x d matrix of n source points in d dimensions.
//    * mins --> d x 1 vector of minimum values; input values ignored if update == 0
//    * maxs --> d x 1 vector of maximum values; input values ignored if update == 0
//    * update --> if set to 1, then max[i] will contain 
//          max(max of values of all samples along dimension i, max[i] input value), and 
//          similarly for min[i].
//
// Output
//    * mins --> d x 1 vector of minimum values
//    * maxs --> d x 1 vector of maximum values
FIGTREE_DLL_API
int figtreeCalcMinMax( int d, int n, double * x, double * mins, double * maxs, int update=0 );

FIGTREE_DLL_API
int figtreeCalcMaxRange( double d, double * mins, double * maxs, double * R );

#if defined(__cplusplus) && !defined(FIGTREE_STATIC)
} // extern "C"
#endif

#endif // FIGTREE_INTERNAL_H
