#!/usr/bin/env python
#=====================================================================================
#                     The UQ Toolkit (UQTk) version @UQTKVERSION@
#                     Copyright (@UQTKYEAR@) Sandia Corporation
#                     http://www.sandia.gov/UQToolkit/
#
#     Copyright (2013) Sandia Corporation. Under the terms of Contract DE-AC04-94AL85000
#     with Sandia Corporation, the U.S. Government retains certain rights in this software.
#
#     This file is part of The UQ Toolkit (UQTk)
#
#     UQTk is free software: you can redistribute it and/or modify
#     it under the terms of the GNU Lesser General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
#
#     UQTk is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU Lesser General Public License for more details.
#
#     You should have received a copy of the GNU Lesser General Public License
#     along with UQTk.  If not, see <http://www.gnu.org/licenses/>.
#
#     Questions? Contact Bert Debusschere <bjdebus@sandia.gov>
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

rc('legend',loc='upper left', fontsize=22)
rc('lines', linewidth=4, color='r')
rc('axes',linewidth=3,grid=True,labelsize=22)
rc('xtick',labelsize=20)
rc('ytick',labelsize=20)

plotdata=True
datafile="inputdata.dat"
# Read data file into numpy array
odeData = npy.loadtxt("solution_ref.dat")

# Define stride to reduce data
stride = 1
# extract x coordinates
t = odeData[::stride,0]
u = odeData[::stride,1]
v = odeData[::stride,2]
w = odeData[::stride,3]

# create figure and axis
fig = plt.figure(figsize=(8,8))
ax = fig.add_axes([0.15, 0.15, 0.75, 0.75])

# set axis limits
#ax.set_xlim([0,t[-1]])
#ax.set_ylim([0,0.2])


plt.plot(u,v, color='purple',label='reference run')
if (plotdata):
    data=npy.genfromtxt(datafile)
    plt.plot(data[:,0], data[:,1], 'o', color='black',markersize=8, label='data')
    plt.title("Data and reference run", fontsize=22)
    plt.legend(loc="upper right")

plt.xlabel("Species u")
plt.ylabel("Species v")
plt.savefig("uv_data.eps")
plt.clf()                  

################################################
fig = plt.figure(figsize=(8,8))
ax = fig.add_axes([0.15, 0.15, 0.75, 0.75])

ax.set_xlim([0,t[-1]*1.1])
plt.plot(t,u,label='u')
plt.plot(t,v,label='v')
plt.plot(t,w,label='w')
if (plotdata):
    data=npy.genfromtxt(datafile)
    ndata=data.shape[0]
    plt.plot([t[-1]]*ndata, data[:,0], 'o', markersize=8, color='blue',label='u data')
    plt.plot([t[-1]]*ndata, data[:,1], 'o', markersize=8, color='green',label='v data')
    plt.title("Data and reference run", fontsize=22)
plt.legend(loc="upper center")

plt.xlabel("Time t")
plt.ylabel("Species concentrations")
plt.savefig("uvwt_data.eps")

