#=====================================================================================
#
#                      The UQ Toolkit (UQTk) version 3.1.1
#                          Copyright (2021) NTESS
#                        https://www.sandia.gov/UQToolkit/
#                        https://github.com/sandialabs/UQTk
#
#     Copyright 2021 National Technology & Engineering Solutions of Sandia, LLC (NTESS).
#     Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government
#     retains certain rights in this software.
#
#     This file is part of The UQ Toolkit (UQTk)
#
#     UQTk is open source software: you can redistribute it and/or modify
#     it under the terms of BSD 3-Clause License
#
#     UQTk is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     BSD 3 Clause License for more details.
#
#     You should have received a copy of the BSD 3 Clause License
#     along with UQTk. If not, see https://choosealicense.com/licenses/bsd-3-clause/.
#
#     Questions? Contact the UQTk Developers at <uqtk-developers@software.sandia.gov>
#     Sandia National Laboratories, Livermore, CA, USA
#=====================================================================================

# Import Modules
try:
    import numpy as np
except ImportError:
    print("Numpy module could not be found")

################################################################################
def pcorr(spl):
    """
    Compute Pearson's correlation coefficient between random vectors

    Args:
       spl - 2D array of samples, first dimensions is the number of samples,
            second dimension is the number of random vectors

    Output:
       Returns a 2D array of correlation coefficients between pairs of random vectors;
       only entries 0<=j<i<no. of random vectors are populated

    References:
       http://en.wikipedia.org/wiki/Pearson_product-moment_correlation_coefficient

    Author:
       Cosmin Safta <csafta@sandia.gov>
     """
    nspl  = spl.shape[0]
    nvars = spl.shape[1]
    splmn  = np.mean(spl,axis=0)
    splstd = np.std(spl,axis=0)
    corrcoef = np.zeros((nvars,nvars))
    for i in range(1,nvars):
        for j in range(i):
            rxy=np.mean((spl[:,i]-splmn[i])*(spl[:,j]-splmn[j]))/(splstd[i]*splstd[j]);
            corrcoef[i,j] = rxy
    return corrcoef
################################################################################
def distcorr(spl):
    """
    Compute distance correlation between random vectors

    Args:
       spl - 2D array of samples, first dimensions is the number of samples,
             second dimension is the number of random vectors

    Output:
       Returns a 2D array of distance correlations between pairs of random vectors;
       only entries 0<=j<i<no. of random vectors are populated

    References:
       http://en.wikipedia.org/wiki/Distance_correlation

    Author:
       Cosmin Safta <csafta@sandia.gov>
    """
    nspl  = spl.shape[0]
    nvars = spl.shape[1]
    if (nspl>5000):
        print('Warning ! This might be a lengthy calculation: nspl='+nspl)
    As=[]
    for i in range(nvars):
        Amat = np.matrix(np.zeros((nspl,nspl)))
        for i1 in range(nspl):
            for j1 in range(nspl):
                Amat[i1,j1] = abs(spl[i1,i]-spl[j1,i])
        # compute means
        Alin = np.array([np.mean(Amat[i1,:]) for i1 in range(nspl)])
        Acol = np.array([np.mean(Amat[:,j1]) for j1 in range(nspl)])
        Amn = np.mean(Amat)
        # subtract/add means (linewise, columnwise, overall)
        Amat = Amat - Alin.reshape(nspl,1)
        Amat = (Amat.T - Acol.reshape(nspl,1)).T
        Amat = Amat+Amn
        As.append(Amat.copy())
    dCor = np.zeros((nvars,nvars))
    dVarX = [np.sqrt(np.sum(np.multiply(As[i],As[i]))/(nspl*nspl)) for i in range(nvars)]
    for i in range(1,nvars):
        for j in range(i):
            dCov   = np.sqrt(np.sum(np.multiply(As[i],As[j]))/(nspl*nspl))
            dCor[i,j] = dCov/np.sqrt(dVarX[i]*dVarX[j])
    As=[]
    return dCor
