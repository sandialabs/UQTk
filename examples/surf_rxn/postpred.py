#!/usr/bin/env python
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

import os
import sys
import shutil
import fileinput
import numpy as npy
import matplotlib.pyplot as plt
from scipy import stats, mgrid, c_, reshape, random, rot90
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

from pylab import *

rc('legend',loc='upper left', fontsize=12)
rc('lines', linewidth=4, color='r')
rc('axes',linewidth=3,grid=True,labelsize=22)
rc('xtick',labelsize=20)
rc('ytick',labelsize=20)


plot_w=True
zoom=False
chdimlist=[map(int, line.split()) for line in open("chdimlist",'r').readline() if line.strip()]
datafile="inputdata.dat"
chainfile="chain.dat"; # chain file
nskip=5000;            # skip first 'nskip' states
nthin=100;            # pick every 'nthin' state
ndim=6
paramxml_file='surf_rxn.par.xml'
sigma=0.1;             # sigma ratio

# load chain and thin
sol = npy.genfromtxt(chainfile)
chain = sol[nskip+1::nthin,1:-2]
nsam=chain.shape[0]
chdim=chain.shape[1]


# run model
tst=1002
uu=npy.empty((tst,))
vv=npy.empty((tst,))
ww=npy.empty((tst,))
for ip in range(nsam):
    
    shutil.copyfile(paramxml_file,'surf_rxn.in.xml')
    ii=0

    for idim in chdimlist:
        for line in fileinput.input("surf_rxn.in.xml", inplace = 1):
            print line.replace('PAR_'+str(idim[0]), str(chain[ip,ii])),
        ii=ii+1
    
 
    
    os.system('./SurfRxnDet.x > det.log')
    cursol=npy.loadtxt("solution.dat")
    #uvw_ss=cursol[-1,1:]
    
    uu=npy.vstack((uu,cursol[:,1]))
    vv=npy.vstack((vv,cursol[:,2]))
    ww=npy.vstack((ww,cursol[:,3]))

uu=uu[1:,:]
vv=vv[1:,:]
ww=ww[1:,:]

times=cursol[:,0]
u_ave=npy.average(uu,axis=0)
u_std=npy.std(uu,axis=0)
v_ave=npy.average(vv,axis=0)
v_std=npy.std(vv,axis=0)
w_ave=npy.average(ww,axis=0)
w_std=npy.std(ww,axis=0)

u_pps=(u_std**2+sigma**2*u_ave**2 +sigma**2*u_std**2)**0.5
v_pps=(v_std**2+sigma**2*v_ave**2 +sigma**2*v_std**2)**0.5
w_pps=(w_std**2+sigma**2*w_ave**2 +sigma**2*w_std**2)**0.5

fig = plt.figure(figsize=(10,7))
ax=fig.add_axes([0.10,0.15,0.85,0.75])

plu=plt.fill_between(times,u_ave-u_pps,u_ave+u_pps,color='grey')
plv=plt.fill_between(times,v_ave-v_pps,v_ave+v_pps,color='grey')
if (plot_w):
    plw=plt.fill_between(times,w_ave-w_pps,w_ave+w_pps,color='grey')
plt.plot(times, u_ave, label='u')
plt.plot(times, v_ave, label='v')
if (plot_w):
    plt.plot(times, w_ave, label='w')


# Plot data
data=npy.genfromtxt(datafile)
ndata=data.shape[0]

plt.plot([times[-1]]*ndata, data[:,0], 'o', markersize=8, color='blue',label='u data')
plt.plot([times[-1]]*ndata, data[:,1], 'o', markersize=8, color='green',label='v data')

ax.set_title("Posterior predictive distributions",fontsize=22)
ax.set_xlim([0,tst*1.1])
plt.legend(loc="upper left")


ax.set_xlabel("Time")
ax.set_ylabel("Species concentration")


#if(zoom):
ax.set_ylim([0,0.4])
plt.savefig('postpred_zoom.eps')

ax.set_ylim([0,1])
plt.savefig('postpred.eps')



