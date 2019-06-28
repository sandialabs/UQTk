from sklearn.cluster import SpectralClustering, KMeans
from collections import deque
import numpy as np
import sys
import os
import math
import time

def genProposalCovariance(data, clusterLabels):
    lambdas = list()
    clusterMap = dict()
    ndims = data.shape[1]
    n_clusters = max(clusterLabels) + 1

    for i in range(n_clusters):
        samplesInCluster = len(list(filter(lambda x: x == i, clusterLabels)))
        lambdas.append(samplesInCluster)

        clusterMap[i] = np.empty((samplesInCluster, ndims))
        count = 0

        for idx in range(len(clusterLabels)):
            if (clusterLabels[idx] == i):
                clusterMap[i][count] = data[idx]
                count += 1

        # print("Size of cluster: ", i, " ", samplesInCluster)

    lambdas = list(map(lambda x: x / sum(lambdas), lambdas))

    proposal = np.empty((ndims, ndims))

    for i in range(n_clusters):
        cvmat = lambdas[i] * np.cov(clusterMap[i], rowvar = False)
        print("Cluster Covariance:\n", cvmat)
        proposal += cvmat

    np.savetxt("Cov.dat", proposal)
    # w, v = np.linalg.eig(proposal)
    # print("Covariance Eigenvalues: ", w)


def main(datafileName):
    # Sample dat file, rows are samples columns are dimensions
    datafile = datafileName

    data = np.genfromtxt(datafile, unpack=False)

    nsamples = data.shape[0]
    ndims = data.shape[1]

    n = 200; # Rule of thumb
    print("Initial KMeans Clusters: ", n)

    time0 = time.time()
    kmeans = KMeans(n_clusters = n, max_iter = 500).fit(data)
    print("KMeans complete")
    time1 = time.time()
    print("KMeans Time:", time1 - time0)

    cluster = SpectralClustering(eigen_solver='arpack',
                                 affinity='nearest_neighbors',
                                 n_neighbors = 7,
                                 n_init = 20).fit(kmeans.cluster_centers_)

    print("Spectral Clustering complete")
    time2 =time.time()
    print("Spectral Clustering Time:", time2 - time1)

    clusterLabels = deque()

    for idx in range(len(data)):
        clusterLabels.append(cluster.labels_[kmeans.labels_[idx]])

    n_clusters = max(cluster.labels_) + 1
    print("Number of Clusters: ", n_clusters)
    np.savetxt("ClusterLabels.dat", clusterLabels)

    genProposalCovariance(data, clusterLabels)


print("=="*40)
print("Cluster Analysis of: ", sys.argv[1])
print("=="*40)
main(sys.argv[1])
