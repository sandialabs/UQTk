#!/usr/bin/env python
#=====================================================================================
#
#                      The UQ Toolkit (UQTk) version @UQTKVERSION@
#                          Copyright (@UQTKYEAR@) NTESS
#                        https://www.sandia.gov/UQToolkit/
#                        https://github.com/sandialabs/UQTk
#
#     Copyright @UQTKYEAR@ National Technology & Engineering Solutions of Sandia, LLC (NTESS).
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

import os
import shutil
import sys
import numpy as np
import math
import random as rnd
from scipy import stats, mgrid, c_, reshape, random, rot90
import matplotlib.pyplot as plt
from utils import get_npc

import fileinput

from pylab import *

rc('legend',loc='upper left', fontsize=22)
rc('lines', linewidth=4, color='r')
rc('axes',linewidth=3,grid=True,labelsize=22)
rc('xtick',labelsize=20)
rc('ytick',labelsize=20)

# define uqtkbin
if os.environ.get("UQTK_INS") is None:
    print("Error: Need to set path to UQTk install direactory as environment variable UQTK_INS -> Abort")
    quit()

else:
    if ( not os.path.isdir(os.environ["UQTK_INS"]) ):
        print("\"",os.environ["UQTK_INS"],"\" is not a valid path -> Abort")
        quit()

uqtkbin=os.environ["UQTK_INS"]+"/bin"
pcerv=uqtkbin+"/pce_rv"


species  = sys.argv[1] # u v w
qoi      = sys.argv[2] # ave tf
pctype   = sys.argv[3]
ord      = int(sys.argv[4])
methlist = sys.argv[5:]

nsam = 100000 # number of random samples
npts = 200    # number of points for KDE estimates

if (species=='u'):
    spid=1
    lims=[0,0.5]
elif (species=='v'):
    spid=2
    lims=[0,0.15]
elif (species=='w'):
    spid=3
    lims=[0.4,1.0]

fig = plt.figure(figsize=(8,6))
ax=fig.add_axes([0.13,0.12,0.83,0.85])

for method in methlist:
    sol=np.loadtxt("solution_"+method+"_modes.dat")


    # Get the second half of the time series 
    tail=sol.shape[0]/2

    # Compute the average of the given species
    if (qoi=='ave'):
        ave=np.average(sol[-tail:,1:],axis=0)
        ncol=ave.shape[0]
        npc=ncol/3
        np.savetxt("pccf.dat",ave[npc*(spid-1):npc*spid])
    elif (qoi=='tf'):
        sp_tf=np.array(sol[-1,1:])
        ncol=sp_tf.shape[0]
        npc=ncol/3
        np.savetxt("pccf.dat",sp_tf[npc*(spid-1):npc*spid])

    # Find (in a very unpleasant way) the stochastic dimension
    for dim in range(1,7):
        npcc=get_npc(dim,ord)
        if (npc==npcc):
            break
    
    pcerv=uqtkbin+"/pce_rv"
    os.system(pcerv+" -w'PC' -f'pccf.dat' -x" + pctype + " -d1 -n" + str(nsam) +" -p"+str(dim)+" -o"+str(ord))
    spls=np.genfromtxt("rvar.dat")
    #xlin=np.linspace(spls.min(),spls.max(),npts) ;
    xlin=np.linspace(lims[0],lims[1],npts) ;
    kernlin=stats.kde.gaussian_kde(spls);
    pdflin1=kernlin(xlin);

    plt.plot(xlin,pdflin1,linewidth=2,label=method)

ax.set_xlim(lims)
if qoi == 'ave':
    ax.set_xlabel('QoI: average '+species,fontsize=24)
else:
    ax.set_xlabel('QoI: '+species+' @ final time',fontsize=24)

ax.set_ylabel("PDF("+species+")",fontsize=24)
#plt.title("Species: "+species+", QoI: "+qoi,fontsize=24)
if species == 'u':
    plt.legend(loc="upper left")
else:
    plt.legend(loc="upper right")
plt.savefig("forUQ_"+species+"_"+qoi+"_PCEdens.pdf")





