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
#     Questions? Contact the UQTk Developers at https://github.com/sandialabs/UQTk/discussions
#     Sandia National Laboratories, Livermore, CA, USA
#=====================================================================================

import os
import sys 
import numpy as np
from   scipy import stats
import matplotlib
import matplotlib.pylab as plt

import plot_utils as ut

# Read data file into numpy array
odeData = ut.ReadDataFile("solution_NISP_modes.dat")

# Define stride to reduce data
stride = 1

# extract x coordinates
xCoords = odeData[::stride,0]

# font and linewidth parameters

# lw is line width
# fs is font size
lw,fs = ut.SetPlotParams()

# create figure and axis
fig = plt.figure(figsize=(6,4))
ax = fig.add_axes([0.15, 0.15, 0.75, 0.75])

# set axis limits
ax.set_xlim([0,1000])
ax.set_ylim([-0.1,0.55])

# plot 6 modes of u
ord = 5
pleg = []
for i in range(ord+1):
    pleg.append(plt.plot(xCoords,odeData[::stride,i+1],linewidth=lw)) # u_i

plt.xlabel("Time [-]",fontsize=fs)
plt.ylabel("$u_i$ [-]",fontsize=fs)

leg = plt.legend( (pleg[0][0], pleg[1][0], pleg[2][0], pleg[3][0], pleg[4][0], pleg[5][0]),
                  (r"$u_0$"  , r"$u_1$"  , r"$u_2$"  , r"$u_3$"  , r"$u_4$", r"$u_5$"    ), loc='upper center', ncol=3)

plt.savefig("surf_rxn_NISP_umodes.pdf")

