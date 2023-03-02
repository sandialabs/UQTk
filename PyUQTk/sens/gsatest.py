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

try:
    import matplotlib
    foundMPL = True
except ImportError:
    foundMPL = False

import os
from gsalib import genSpl_Si, genSens_Si, genSpl_Sij, genSens_Sij, genSpl_SiT, genSens_SiT
if foundMPL:
    import matplotlib.pyplot as plt

def func1(x):
    return npy.sum(x)

def func2(x):
    f = npy.sum(x)
    for i in range(len(x)-1):
        f = f+(i+1)*(i+1)*x[i]*x[i+1]
    return f

#--------------------------------------------------------------------------------
#       Setup
#--------------------------------------------------------------------------------
nruns= 10
ndim = 4
nspl = 10000
abr  = npy.zeros((ndim,2))
abr[:,1]=1.0; abr[0,1]=3.0

#--------------------------------------------------------------------------------
#      Theoretical values for Sobol indices
#--------------------------------------------------------------------------------
Simath=npy.array([0.14908, 0.14908, 0.41411, 0.222699])
Sijmath=npy.array([0.00552147, 0.00981595, 0.04969339])
SiTmath=Simath.copy()
for i in range(ndim-1):
    SiTmath[i]   = SiTmath[i]  +Sijmath[i]
    SiTmath[i+1] = SiTmath[i+1]+Sijmath[i]


#--------------------------------------------------------------------------------
#      First-order and joint Sobol indices
#--------------------------------------------------------------------------------
runSi=[]
runSij=[]
print('Computing first-order and joint Sobol indices')
for irun in range(nruns):
    print('  run %d out of %d'%(irun+1,nruns))
    # Generate samples for Si
    genSpl_Si(nspl,ndim,abr,nd=12)
    # Load samples, run model, save model evaluations
    gsaens=npy.genfromtxt('gsaSplSi.dat')
    modelRuns=npy.array([func2(gsaens[i]) for i in range(gsaens.shape[0])])
    npy.savetxt('modelSi.dat', modelRuns, fmt="%.12e", delimiter=' ', newline='\n')
    # Compute first order sensitivity indices
    Si=genSens_Si('modelSi.dat',ndim,verb=0)
    # Generate samples for Sij
    genSpl_Sij(ndim,matfile='mat12.npz',nd=12)
    # Load samples, run model, save model evaluations
    gsaens=npy.genfromtxt('gsaSplSij.dat')
    modelRuns=npy.array([func2(gsaens[i]) for i in range(gsaens.shape[0])])
    npy.savetxt('modelSij.dat', modelRuns, fmt="%.12e", delimiter=' ', newline='\n')
    # Compute joint sensitivity indices
    Sij=genSens_Sij(Si,'modelSij.dat',verb=0)
    runSi.append(Si)
    runSij.append(Sij)


runSi=npy.array(runSi)
runSij=npy.array(runSij)

if foundMPL:
    width = 0.4
    ind = npy.arange(ndim)+0.5-width/2.0
    fs1 = 24
    # Si
    fig=plt.figure(figsize=(8,6))
    ax=fig.add_axes([0.15,0.10,0.8,0.85])
    Simean = npy.array([npy.average(runSi[:,i]) for i in range(ndim)])
    Sistd  = npy.array([npy.std(runSi[:,i]) for i in range(ndim)])
    rects1 = ax.bar(ind, Simean, width, color='r', yerr=Sistd,error_kw=dict(linewidth=3, color='b',capsize=5) )
    plt.plot(ind+width/2.0,Simath,'o',ms=8,mfc='k')
    ax.set_xlim([0,ndim])
    ax.set_ylim([0,0.5])
    ax.set_ylabel(r'$S_i$',fontsize=fs1)
    ax.set_xticks(ind+width/2)
    ax.set_yticks([0,0.1,0.2,0.3,0.4,0.5])
    ax.set_yticklabels( ('$0$', '$0.1$', '$0.2$', '$0.3$', '$0.4$', '$0.5$') ,fontsize=fs1-6)
    ax.set_xticklabels( ('$x_1$', '$x_2$', '$x_3$', '$x_4$') ,fontsize=fs1)
    plt.savefig('gsaspl_Si.pdf')
    #plt.show()
    # Sij
    width = 0.4
    ind = npy.arange(ndim-1)+0.5-width/2.0
    fs1 = 24
    fig=plt.figure(figsize=(8,6))
    ax=fig.add_axes([0.15,0.10,0.8,0.85])
    Simean = npy.array([npy.average(runSij[:,i,i+1]) for i in range(ndim-1)])
    Sistd  = npy.array([npy.std(runSij[:,i,i+1]) for i in range(ndim-1)])
    rects1 = ax.bar(ind, Simean, width, color='r', yerr=Sistd,error_kw=dict(linewidth=3, color='b',capsize=5) )
    plt.plot(ind+width/2.0,Sijmath,'o',ms=8,mfc='k')
    ax.set_xlim([0,ndim-1])
    ax.set_ylim([0,0.07])
    ax.set_ylabel(r'$S_{ij}$',fontsize=fs1)
    ax.set_xticks(ind+width/2)
    ax.set_yticks([0,0.02,0.04,0.06])
    ax.set_yticklabels( ('$0$', '$0.02$', '$0.04$', '$0.06$') ,fontsize=fs1-6)
    ax.set_xticklabels( ('$(x_1,x_2)$', '$(x_2,x_3)$', '$(x_3,x_4)$') ,fontsize=fs1)
    plt.savefig('gsaspl_Sij.pdf')
    #plt.show()
else:
    # could not find matplotlib, saving Sobol indices to file
    npy.savez("gsaSi_Sij.npz",Si=runSi,Sij=runSij)

#--------------------------------------------------------------------------------
#      Total-order Sobol indices
#--------------------------------------------------------------------------------
runSiT_1=[]
runSiT_2=[]
print('Computing total-order Sobol indices')
for irun in range(nruns):
    print('  run %d out of %d',(irun+1,nruns))
    # Generate samples for SiT
    genSpl_SiT(nspl,ndim,abr,nd=12)
    # Load samples, run model, save model evaluations
    gsaens=npy.genfromtxt('gsaSplSiT.dat')
    modelRuns=npy.array([func2(gsaens[i]) for i in range(gsaens.shape[0])])
    npy.savetxt('modelSiT.dat', modelRuns, fmt="%.12e", delimiter=' ', newline='\n')
    # Compute total sensitivity indices
    runSiT_1.append(genSens_SiT('modelSiT.dat',ndim,type='type1',verb=0))
    runSiT_2.append(genSens_SiT('modelSiT.dat',ndim,type='type2',verb=0))

runSiT_1=npy.array(runSiT_1)
runSiT_2=npy.array(runSiT_2)

if foundMPL:
    width = 0.4
    ind = npy.arange(ndim)+0.5-width/2.0
    fs1 = 24
    # SiT
    fig=plt.figure(figsize=(8,6))
    ax=fig.add_axes([0.15,0.10,0.8,0.85])
    SiT1mn  = npy.array([npy.average(runSiT_1[:,i]) for i in range(ndim)])
    SiT1std = npy.array([npy.std(runSiT_1[:,i]) for i in range(ndim)])
    rects1 = ax.bar(ind, SiT1mn, width/2.0, color='r', yerr=SiT1std,error_kw=dict(linewidth=3, color='b',capsize=5),label="Est.1")
    SiT2mn  = npy.array([npy.average(runSiT_2[:,i]) for i in range(ndim)])
    SiT2std = npy.array([npy.std(runSiT_2[:,i]) for i in range(ndim)])
    rects2 = ax.bar(ind+width/2.0, SiT2mn, width/2.0, color='y', yerr=SiT2std,error_kw=dict(linewidth=3, color='b',capsize=5),label="Est.2")
    plt.plot(ind+width/2.0,SiTmath,'o',ms=8,mfc='k',label="Exact")
    ax.set_xlim([0,ndim])
    ax.set_ylim([0,0.55])
    ax.set_ylabel(r'$S_i^T$',fontsize=fs1)
    ax.set_xticks(ind+width/2)
    ax.set_yticks([0,0.1,0.2,0.3,0.4,0.5])
    ax.set_yticklabels( ('$0$', '$0.1$', '$0.2$', '$0.3$', '$0.4$', '$0.5$') ,fontsize=fs1-6)
    ax.set_xticklabels( ('$x_1$', '$x_2$', '$x_3$', '$x_4$'),fontsize=fs1)
    plt.legend(loc=2,prop={'size':fs1})
    plt.savefig('gsaspl_SiT.pdf')
    #plt.show()
else:
     # could not find matplotlib, saving Sobol indices to file
    npy.savez("gsaSiT.npz",SiT1=runSiT_1,SiT2=runSiT_2)
