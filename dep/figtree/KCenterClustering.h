//-------------------------------------------------------------------
// This code was modified by Vlad Morariu:
// 11/03/06: 
//   Removed references to Matlab to compile code into a library
// 01/24/07:
//   KCenterClustering now has the ability to increase the number of 
//   clusters incrementally, calculating the max cluster radius at each
//   iteration.  
// 02/07/07:
//   Clustering now stops when the max cluster radius
//   is zero (when number of  clusters has reached the number of 
//   unique points), and the number of ACTUAL clusters used is returned.
//-------------------------------------------------------------------

//-------------------------------------------------------------------
// The code was written by Changjiang Yang and Vikas Raykar
// and is copyrighted under the Lesser GPL: 
//
// Copyright (C) 2006  Changjiang Yang and Vikas Raykar
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
// The author may be contacted via email at:cyang(at)sarnoff(.)com
// vikas(at)umiacs(.)umd(.)edu
//-------------------------------------------------------------------

//----------------------------------------------------------------------------
// File    : KCenterClustering.h
// Purpose : Interface for the k-center clustering algorithm.
// Author  : Vikas C. Raykar (vikas@cs.umd.edu)
// Date    : April 25 2005, June 10 2005, August 23, 2005
//
//----------------------------------------------------------------------------
// Gonzalez's farthest-point clustering algorithm.
//
// June 10, 2005: 
// This version now returns the number points and the radius of each cluster.
//
// August 23, 2005:
// Speed up using the doubly circular list.
// The clusters far away are trimmed. The nodes inside the neighboring
// clusters which are within half sphere are trimmed.
//
//----------------------------------------------------------------------------
//
// INPUT 
// ----------------
//
// Dim               --> dimension of the points.
// NSources          --> number of sources.
// pSources          --> pointer to sources, (d*N).
// NumClusters       --> number of clusters.
//
// OUTPUT
// ----------------
//
// MaxClusterRadius  --> maximum radius of the clusters, (rx).
// pClusterIndex     --> vector of length N where the i th element is the
//                     cluster number to which the i th point belongs.
//                     pClusterIndex[i] varies between 0 to K-1. 
// pClusterCenters   --> pointer to the cluster centers, (d*K). 
// pNumPoints        --> pointer to the number of points in each cluster, (K).
// pClusterRadii     --> pointer to the radius of each cluster, (K). 
//----------------------------------------------------------------------------

#ifndef K_CENTER_CLUSTERING_H
#define K_CENTER_CLUSTERING_H

class KCenterClustering{
  public:    

    //Output parameters

    double MaxClusterRadius;  //maximum cluster radius
  
    //Functions

    //constructor 
    KCenterClustering(int Dim,
      int NSources,
      double *pSources,
      int *pClusterIndex,
      int NumClusters
    );

    //destructor
    ~KCenterClustering();
    
    //K-center clustering
    //Returns the number of actual clusters (it might have stopped early if all clusters have
    //    radius of 0 -- which means that the number of clusters has reached the number
    //    of unique pts)
    int Cluster();

    //Incremental k-center clustering
    //  nClusters - if non-NULL, value is set to the # of clusters at end of call
    //  maxRadius - if non-NULL, value is set to the max radius of all clusters
    void ClusterIncrement( int * nClusters, double * maxRadius );

    //Compute cluster centers and the number of points in each cluster
    //and the radius of each cluster.
    void ComputeClusterCenters( int NumClusters,
                                double *pClusterCenters,
                                int *pNumPoints,
                                double *pClusterRadii
                              );
    
  private:
    //Input Parameters

    int d;           // dimension of the points.
    int N;           // number of sources.
    double *px;      // pointer to sources, (d*N).
    int K;           // max number of clusters
    int *pci;        // pointer to a vector of length N where the i th element is the 
                     // cluster number to which the i th point belongs.
    double *dist_C;  // distances to the center.
    double *r;

    int *pCenters;   // indices of the centers.
    int *cprev;      // index to the previous node
    int *cnext;      // index to the next node
    int *far2c;      // farthest node to the center

    int numClusters; // added by Vlad to keep track of # of clusters

    //Functions
    double ddist(const int d, const double *x, const double *y);
    int idmax(int n, double *x);
    
};


#endif