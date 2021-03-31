#!/usr/bin/env python
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

"""
Regression-related tools
"""

import os
import sys

try:
    import numpy as np
except ImportError:
    print('Numpy was not found.')


from .multiindex import mi_addfront

#############################################################
#############################################################
#############################################################

def regression(xdata,ydata,mode,basisparams,regparams):
    """
    Polynomial regression given x- and y-data. A wrapper around the regression app

    Arguments:
        * xdata       : N x d array of x-data
        * ydata       : N x e array of y-data
        * mode        : m (only mean), ms (mean and stdev) or msc (mean,stdev and cov)
        * basisparams : tuple of two elements, (pctype,mindex)
                      : pctype - PC type
                      : mindex - multiindex array
        * regparams   : tuple of two elements, (method, methodpars)
                      : method - regression method, 'lsq' or 'wbcs'
                      : methodpars - parameters of the regression (regularization weights for wbcs)
                      : Note: regparams=np.array(methodpars) on output

    Returns:
        * cfs     : Coefficient vector
        * mindex  : Multiindex array
        * Sig     : Coefficient covariance matrix
        * used    : Indices of retained multiindices
    """

    # Read input settings
    pctype,mindex=basisparams
    method,methodpars=regparams

    # Turn regparams tuple into an array
    regparams=np.array(methodpars).reshape(-1,1)

    # Get the dimensionality
    dim=mindex.shape[1]

    # Save the appropriate files for the regression app
    np.savetxt('xdata.dat',xdata)
    np.savetxt('ydata.dat',ydata)
    np.savetxt('mindex.dat',mindex,fmt='%d')
    np.savetxt('regparams.dat',regparams,fmt='%24.16f')

    # Regularization
    lamstr=''
    regstr=''
    if method=='lsq':
        lamstr='-l 0.0'
    if method=='wbcs':
        regstr='-w regparams.dat'

    # Run the regression app
    cmd='regression -c 1.e-5 -x xdata.dat -y ydata.dat -b PC_MI -s '+pctype+' -p mindex.dat -m '+mode+' -r '+method + ' '+regstr+' '+lamstr +' > regr.log'
    print('Running %s'%(cmd))
    os.system(cmd)

    # Read the resulting files
    cfs=np.loadtxt('coeff.dat')
    used=np.loadtxt('selected.dat',dtype=int)
    mindex=np.loadtxt('mindex.dat',dtype=int).reshape(-1,dim)
    mindex=mindex[used]
    if (mode=='msc'):
        Sig=np.loadtxt('Sig.dat')
    else:
        Sig=[]

    # Return coefficient, multiindex, coef. covariance matrix, and indices of used basis terms
    return (cfs,mindex,Sig,used)

#############################################################
#############################################################
#############################################################

def regression_iter(xdata,ydata,mode,basisparams,regparams,iterparams):
    """
    Iterative regression involving multiindex growth.
    See inputs and outputs of regression(), with additional argument

    iterparams : a tuple (niter,eps,update_weights,update_mindex)
        * niter : Number of iterations
        * eps   : Nugget for iterative reweighting
        * update_weights : boolean flag whether to recompute the weights or not
        * update_mindex  : boolean flag whether to update multiindex or not

    """

    # Read the inputs
    pctype,mindex=basisparams
    method,methodpars=regparams
    niter,eps,update_weights,update_mindex=iterparams

    # Set the current parameters
    basisparams_cur=[pctype,mindex]
    regparams_cur=[method,methodpars]


    nrange=np.arange(mindex.shape[0])
    cur_used=nrange
    npc=mindex.shape[0]
    for i in range(niter):
        print('Iteration %d / %d'%(i+1,niter))
        print('Initial mindex size %d'%(basisparams_cur[1].shape[0]))
        cfs_cur,mindex_cur,Sig,used=regression(xdata,ydata,mode,basisparams_cur,regparams_cur)
        print('New mindex size     %d'%(mindex_cur.shape[0]))
        #tmp=cur_used[used]
        #cur_used=tmp.copy()
        npc_cur=mindex_cur.shape[0]
        # Update weights or not
        if (update_weights==True):
            regparams_cur[1]=1./(abs(cfs_cur)+eps)
        else:
            tmp=regparams_cur[1]
            regparams_cur[1]=tmp[list(used)] #read used.dat and replace it here

        # Update multiindex or not
        if (update_mindex==True and i<niter-1):
            mindex_new,mindex_add,mindex_f=mi_addfront(mindex_cur)

            mindex_cur=mindex_new.copy()
            basisparams_cur[1]=mindex_new
            regparams_new=np.ones(mindex_new.shape[0])/eps
            regparams_new[0:npc_cur]=regparams_cur[1]
            regparams_cur[1]=regparams_new

    # Return coefficient, multiindex, coef. covariance matrix, and indices of used basis terms
    return (cfs_cur,mindex_cur,Sig,used)
