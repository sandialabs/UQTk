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
import sys 
import numpy as np
from   scipy import stats
import matplotlib
import matplotlib.pylab as plt

import plot_utils as ut

method=sys.argv[1]

# Read data file into numpy array
odeData = ut.ReadDataFile("solution_"+method+".dat")

# Define stride to reduce data
stride = 10
# extract x coordinates
xCoords = odeData[::stride,0]

# font and linewidth parameters

# lw is line width
# fs is font size
lw,fs = ut.SetPlotParams()
elw = lw/2 # error bar line width

# create figure and axis
fig = plt.figure(figsize=(6,4))
ax = fig.add_axes([0.15, 0.15, 0.75, 0.75])

# set axis limits
ax.set_xlim([0,1000])
ax.set_ylim([0,1])

# plot mean u, v, w vs. x and add error bars: +/- one std dev.
pleg = []
pleg.append(plt.errorbar(xCoords,odeData[::stride,1],odeData[::stride,2],linewidth=lw,elinewidth=elw)) # u, std dev u
pleg.append(plt.errorbar(xCoords,odeData[::stride,3],odeData[::stride,4],linewidth=lw,elinewidth=elw)) # v, std dev v
pleg.append(plt.errorbar(xCoords,odeData[::stride,5],odeData[::stride,6],linewidth=lw,elinewidth=elw)) # w, std dev w

plt.xlabel("Time [-]",fontsize=fs)
plt.ylabel("Species Mass Fractions [-]",fontsize=fs)

leg = plt.legend( (pleg[0][0], pleg[1][0], pleg[2][0]),
                  (r"$u$"    , r"$v$"    , r"$w$"    ), loc='upper right', ncol=3)

plt.title("Method " + method,fontsize=fs)

plt.savefig("forUQ_surf_rxn_"+method+"_mstd.pdf")
                  


