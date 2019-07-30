#!/usr/bin/env python
#=====================================================================================
#                     The UQ Toolkit (UQTk) version @UQTKVERSION@
#                     Copyright (@UQTKYEAR@) Sandia Corporation
#                     http://www.sandia.gov/UQToolkit/
#
#     Copyright (@UQTKYEAR@) Sandia Corporation. Under the terms of Contract DE-AC04-94AL85000
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
"""
Plots a 2D structured grid using in data saved in
xgrid.dat/ygrid.dat files.
usage ./pl2Dsgrid directory
"""

import os
import sys
import numpy as np
from   scipy import stats
import matplotlib
import matplotlib.pylab as plt
from   pylab import *

from pyutils import readfile, column

if (len(sys.argv) > 1):
    dir  = sys.argv[1]
else:
    print("Need 1 arguments: directory")

x,nx = readfile(dir+"/xgrid.dat");
y,ny = readfile(dir+"/ygrid.dat");
x=column(x,0);
y=column(y,0);

fig=plt.figure(figsize=(4,4))
ax = fig.add_axes([0.05, 0.05, 0.9, 0.9])
ax.set_xticks(x)
ax.set_yticks(y)
ax.set_xticklabels([])
ax.set_yticklabels([])
grid(True)
ticklines = ax.get_xticklines()
ticklines.extend( ax.get_yticklines() )
gridlines = ax.get_xgridlines()
gridlines.extend( ax.get_ygridlines() )
ticklabels = ax.get_xticklabels()
ticklabels.extend( ax.get_yticklabels() )
for line in ticklines:
    line.set_linewidth(3)

for line in gridlines:
    line.set_linewidth(1)
    line.set_color('k')
    line.set_linestyle('-')

plt.savefig("grid2d.pdf")
