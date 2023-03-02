#=====================================================================================
#
#                      The UQ Toolkit (UQTk) version 3.1.3
#                          Copyright (2023) NTESS
#                        https://www.sandia.gov/UQToolkit/
#                        https://github.com/sandialabs/UQTk
#
#     Copyright 2023 National Technology & Engineering Solutions of Sandia, LLC (NTESS).
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
try:
    import numpy as npy
except ImportError:
    print('gsalib requires numpy package -> Exit')
    quit()

import os.path

def genSpl_Si(nspl,ndim,abrng,**kwargs):
    # get default values for optional arguments
    splout  = kwargs.get('splout', "gsaSplSi.dat") # samples file
    matfile = kwargs.get('matfile',"mat12.npz")    # intermediary matrices
    verb    = kwargs.get('verb', 0)                # verbosity
    nd      = kwargs.get('nd',  18)                # no. of significant digits in samples output
    # Test nd values
    if (nd<6) or (nd>18):
        raise ValueError("Number of digits should be between 6 and 18")
    #------------------------------------------------------------------------------------
    # create nspl uniform samples in [a_i,b_i], i=1,ndim
    #------------------------------------------------------------------------------------
    if verb>0:
        print('Create ensemble of input parameters')
    mat1=npy.random.random_sample((nspl,ndim))
    mat2=npy.random.random_sample((nspl,ndim))
    for i in range(ndim):
        mat1[:,i] = mat1[:,i]*(abrng[i,1]-abrng[i,0])+abrng[i,0]
        mat2[:,i] = mat2[:,i]*(abrng[i,1]-abrng[i,0])+abrng[i,0]
    # save temporary matrices
    npy.savez(matfile, mat1=mat1, mat2=mat2)
    # assemble the big matrix for main sensitivities directly to a file
    if os.path.isfile(splout):
        os.remove(splout)
    f_handle = file(splout, 'a')
    npy.savetxt(f_handle, mat1, fmt="%."+str(nd)+"e", delimiter='  ', newline='\n')
    for idim in range(ndim):
        if verb>0:
            print(' - working on parameter %d'%(idim))
        matj=mat2.copy();
        matj[:,idim]=mat1[:,idim]
        npy.savetxt(f_handle, matj, fmt="%."+str(nd)+"e", delimiter='  ', newline='\n')
    npy.savetxt(f_handle, mat2, fmt="%."+str(nd)+"e", delimiter='  ', newline='\n')
    f_handle.close()

def genSens_Si(modeval,ndim,**kwargs):
    # get optional arguments
    verb    = kwargs.get('verb', 0)  # verbosity
    #------------------------------------------------------------------------------------
    # load model evaluations and compute main sensitivities
    #------------------------------------------------------------------------------------
    ymod = npy.genfromtxt(modeval)
    nspl = ymod.shape[0]/(ndim+2)
    if verb > 0:
        print('Compute sensitivities, no. of samples: %d'%(nspl))
    sobolSi = npy.zeros(ndim)
    yMat1   = ymod[:nspl]
    yMat2   = ymod[nspl*(ndim+1):]
    mean12  = npy.mean(yMat1*yMat2)
    vv1     = npy.var(yMat1,ddof=1)
    for idim in range(1,ndim+1):
        vari = npy.sum(yMat1*ymod[idim*nspl:(idim+1)*nspl])
        sobolSi[idim-1] = (vari/(nspl-1.0)-mean12)/vv1
        if verb > 1:
            print(' - parameter %d: %e'%(idim,sobolSi[idim-1]))
    if verb > 0:
        print(' - total first order sensitivity: %e'%(npy.sum(sobolSi)))
    return sobolSi

def genSpl_SiT(nspl,ndim,abrng,**kwargs):
    # get optional arguments
    splout  = kwargs.get('splout', "gsaSplSiT.dat") # samples file
    matfile = kwargs.get('matfile', "mat12.npz")    # intermediary matrices
    verb    = kwargs.get('verb', 0)                 # verbosity
    nd      = kwargs.get('nd', 18)                  # no. of significant digits in samples output
    # Test nd values
    if (nd<6) or (nd>18):
        raise ValueError("Number of digits should be between 6 and 18")
    #------------------------------------------------------------------------------------
    # create nspl uniform samples in [a_i,b_i], i=1,ndim
    #------------------------------------------------------------------------------------
    if verb>0:
        print('Create ensemble of input parameters')
    mat1=npy.random.random_sample((nspl,ndim))
    mat2=npy.random.random_sample((nspl,ndim))
    for i in range(ndim):
        mat1[:,i] = mat1[:,i]*(abrng[i,1]-abrng[i,0])+abrng[i,0]
        mat2[:,i] = mat2[:,i]*(abrng[i,1]-abrng[i,0])+abrng[i,0]
    # save temporary matrices
    npy.savez(matfile, mat1=mat1, mat2=mat2)
    # assemble the big matrix for main sensitivities
    if os.path.isfile(splout):
        os.remove(splout)
    f_handle = file(splout, 'a')
    npy.savetxt(f_handle, mat1, fmt="%."+str(nd)+"e", delimiter='  ', newline='\n')
    for idim in range(ndim):
        if verb>0:
            print(' - working on parameter %d'%(idim))
        matj=mat1.copy();
        matj[:,idim]=mat2[:,idim]
        npy.savetxt(f_handle, matj, fmt="%."+str(nd)+"e", delimiter='  ', newline='\n')
    npy.savetxt(f_handle, mat2, fmt="%."+str(nd)+"e", delimiter='  ', newline='\n')
    f_handle.close()
    return

def genSens_SiT(modeval,ndim,**kwargs):
    # get optional arguments
    verb      = kwargs.get('verb', 0)       # verbosity
    siTmethod = kwargs.get('type', 'type1') # sampling method
    #------------------------------------------------------------------------------------
    # load model evaluations and compute main sensitivities
    #------------------------------------------------------------------------------------
    ymod = npy.genfromtxt(modeval)
    nspl = ymod.shape[0]/(ndim+2)
    if verb > 0:
        print('Compute sensitivities, no. of samples: %d'%(nspl))
    sobolSiT = npy.zeros(ndim)
    yMat1  = ymod[:nspl]
    vv1    = npy.var(yMat1,ddof=1)
    Ey     = npy.average(yMat1)
    for idim in range(1,ndim+1):
        if (siTmethod == "type1"):
            ssqrs=0.0
            for i in range(nspl):
                ssqrs = ssqrs+yMat1[i]*ymod[idim*nspl+i]
            vari = ssqrs/(nspl-1.0)
            sobolSiT[idim-1] = 1-(vari-Ey**2)/(vv1);
        else:
            vari      = npy.sum(npy.power(yMat1-ymod[idim*nspl:(idim+1)*nspl],2))/nspl
            sobolSiT[idim-1] = vari/(2.0*vv1);
        if verb > 1:
            print(' - parameter %d: %e'%(idim,sobolSiT[idim-1]))
    if verb > 0:
        print(' - total main sensitivity: %e'%(npy.sum(sobolSiT)))
    return sobolSiT

def genSpl_SiTcust(nspl,ndim,abrng,collst,**kwargs):
    # get optional arguments
    splout  = kwargs.get('splout', "gsaSplSiT.dat") # samples file
    verb    = kwargs.get('verb', 0)                 # verbosity
    nd      = kwargs.get('nd', 18)                  # no. of significant digits in samples output
    if (nd<6) or (nd>18):
        raise ValueError('Number of digits should be between 6 and 18')
    #------------------------------------------------------------------------------------
    # create nspl uniform samples in [a_i,b_i], i=1,ndim
    #------------------------------------------------------------------------------------
    if verb>0:
        print('Create ensemble of input parameters')
    mat1=npy.random.random_sample((nspl,ndim))
    mat2=npy.random.random_sample((nspl,ndim))
    for i in range(ndim):
        mat1[:,i] = mat1[:,i]*(abrng[i,1]-abrng[i,0])+abrng[i,0]
        mat2[:,i] = mat2[:,i]*(abrng[i,1]-abrng[i,0])+abrng[i,0]
    # assemble the big matrix for main sensitivities
    if os.path.isfile(splout):
        os.remove(splout)
    f_handle = file(splout, 'a')
    npy.savetxt(f_handle, mat1, fmt="%."+str(nd)+"e", delimiter='  ', newline='\n')
    for j in collst:
        print(j)
        matj = mat1.copy()
        matj[:,j] = mat2[:,j].copy()
        npy.savetxt(f_handle, matj, fmt="%."+str(nd)+"e", delimiter='  ', newline='\n')
    f_handle.close()
    return

def genSens_SiTcust(modeval,ndim,collst,**kwargs):
    # get optional arguments
    verb      = kwargs.get('verb', 0)       # verbosity
    siTmethod = kwargs.get('type', 'type1') # sampling method
    #------------------------------------------------------------------------------------
    # load model evaluations and compute main sensitivities
    #------------------------------------------------------------------------------------
    ymod = npy.genfromtxt(modeval)
    nspl = ymod.shape[0]/(1+len(collst))
    print('No. of samples %d'%(nspl))
    if verb > 0:
        print('Compute sensitivities, no. of samples: %d'%(nspl))
    sobolSiT = npy.zeros(len(collst))
    yMat1  = ymod[:nspl]
    Ey     = npy.average(yMat1)
    vv1    = npy.var(yMat1,ddof=1)
    for j in range(len(collst)):
        print(j,collst[j])
        if (siTmethod == "type1"):
            ssqrs=0.0
            for i in range(nspl):
                ssqrs = ssqrs+yMat1[i]*ymod[(j+1)*nspl+i]
            vari = ssqrs/(nspl-1.0)
            sobolSiT[j] = 1.0-(vari-Ey**2)/(vv1);
        else:
            vari = npy.sum(npy.power(yMat1-ymod[(j+1)*nspl:(j+2)*nspl],2))/nspl
            sobolSiT[j] = vari/(2.0*vv1);
    if verb > 0:
        npy.set_printoptions(precision=4)
        print(' - total sensitivities: ')
        print(sobolSiT)
    return sobolSiT

def genSpl_Sij(ndim,**kwargs):
    # get optional arguments
    splout  = kwargs.get('splout', "gsaSplSij.dat") # samples file
    matfile = kwargs.get('matfile', "mat12.npz")    # intermediary matrices
    verb    = kwargs.get('verb', 0)                 # verbosity
    nd      = kwargs.get('nd', 18)                  # no. of significant digits in samples output
    if verb > 0:
        print('Load intermediary matrices of input parameters')
    if os.path.isfile(matfile):
        m12=npy.load(matfile)
    else:
        raise IOError('Could not load samples')
        quit()
    mat1=m12["mat1"]
    mat2=m12["mat2"]
    # assemble the big matrix for main sensitivities
    if os.path.isfile(splout):
        os.remove(splout)
    f_handle = file(splout, 'a')
    npy.savetxt(f_handle, mat1, fmt="%."+str(nd)+"e", delimiter='  ', newline='\n')
    for idim in range(ndim-1):
        for jdim in range(idim+1,ndim):
            if verb>1:
                print(' - working on pair %d,%d'%(idim,jdim))
            matj=mat2.copy();
            matj[:,idim]=mat1[:,idim]
            matj[:,jdim]=mat1[:,jdim]
            npy.savetxt(f_handle, matj, fmt="%."+str(nd)+"e", delimiter='  ', newline='\n')
    npy.savetxt(f_handle, mat2, fmt="%."+str(nd)+"e", delimiter='  ', newline='\n')
    f_handle.close()
    return

def genSens_Sij(sobolSi,modeval,**kwargs):
    # get optional arguments
    verb    = kwargs.get('verb', 0)     # verbosity
    #------------------------------------------------------------------------------------
    # joint sensitivities
    #------------------------------------------------------------------------------------
    ndim = len(sobolSi)
    ymod = npy.genfromtxt(modeval)
    nspl = ymod.shape[0]/(ndim*(ndim-1)/2+2)
    if verb > 0:
        print('No. of samples, no. of dimensions: %d,%d'%(nspl,ndim))
    sobolSij = npy.array(npy.zeros((ndim,ndim)))
    yMat1    = ymod[:nspl]
    yMat2    = ymod[nspl*(ndim*(ndim-1)/2+1):]
    mean12   = npy.mean(yMat1*yMat2)
    vv1      = npy.var(yMat1,ddof=1)
    ijd = 0
    for idim in range(1,ndim):
        for jdim in range(idim+1,ndim+1):
            ijd += 1;
            vari = npy.sum(yMat1*ymod[ijd*nspl:(ijd+1)*nspl])
            sobolSij[idim-1,jdim-1] = (vari/(nspl-1.0)-mean12)/vv1-sobolSi[idim-1]-sobolSi[jdim-1]
            if verb > 1:
                print(' - pair %d,%d: %e'%(idim,jdim,sobolSij[idim-1,jdim-1]))
    if verb > 0:
        print(' - total Sij: %e'%(npy.sum(sobolSij[:,:])))
    return sobolSij
