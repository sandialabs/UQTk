#!/usr/bin/env python
#=====================================================================================
#
#                      The UQ Toolkit (UQTk) version 3.1.4
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
#     Questions? Contact the UQTk Developers at https://github.com/sandialabs/UQTk/discussions
#     Sandia National Laboratories, Livermore, CA, USA
#=====================================================================================

import os
import sys
import numpy as np
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from   matplotlib.lines import Line2D
from   pylab import *
import scipy as spy
from   scipy import stats
import re

################################################################################
def version_cmp(version1, version2):
    """Function to compare version numbers, courtesy of
    https://stackoverflow.com/questions/1714027/version-number-comparison-in-python
    """
    def normalize(v):
        return [int(x) for x in re.sub(r'(\.0+)*$','', v).split(".")]

    def my_cmp(a,b): # Since Python 3 no longer has the cmp() function
        return (a > b) - (a < b)

    return my_cmp(normalize(version1), normalize(version2))
################################################################################

if ( len(sys.argv) > 1 ):
    sample_file_name=sys.argv[1]
else:
    print("Please specify the file name samples.a.dat or samples.loga.dat as argument")
    quit()

# Get scipy version to see if the method set_bandwidth is available
spver=spy.__version__
print("scipy version: ",spver)

# set_bandwidth is available as of scipy version 0.11
min_version = "0.11"
bandwidth_present = True
if version_cmp(min_version,spver) == 1:
    bandwidth_present = False

#
# Load data
#
samples=np.genfromtxt(sample_file_name)
xS=np.linspace(samples.min(),samples.max(),200)
kernsS=stats.gaussian_kde(samples)
pdf1=kernsS(xS)

if bandwidth_present:
    print("Performing KDE with a range of bandwidths:")
    kernsS.set_bandwidth(bw_method=kernsS.factor/2.0)
    pdf2=kernsS(xS);
    kernsS.set_bandwidth(bw_method=kernsS.factor*4.0)
    pdf3=kernsS(xS);

#
#  Make the figure
#
lw1=2
fs1=18
fig = plt.figure(figsize=(8,6))
ax=fig.add_axes([0.10, 0.10, 0.85, 0.85]) ;
l1=plt.plot(xS,pdf1,linewidth=lw1,label="optimal")
if bandwidth_present:
    l2=plt.plot(xS,pdf2,linewidth=lw1,label="optimal/2")
    l3=plt.plot(xS,pdf3,linewidth=lw1,label="optimal*2")

xlb=sample_file_name.replace("samples.","").replace(".dat","")
if xlb=="a":
  xlb="\hat{a}"
else:
  xlb="\log(\hat{a})"
ax.set_xlabel("$"+xlb+"$",fontsize=fs1)
ax.set_ylabel("$PDF("+xlb+")$",fontsize=fs1)
plt.legend(loc="upper left")
#ax.set_xlim([0.0,0.55])
#ax.set_ylim([0.0,5.0])
#ax.set_xticks([0,0.1,0.2,0.3,0.4,0.5])
#ax.set_yticks([0,1,2,3,4])
#leg=plt.legend( (l1[0],l2[0]),
#                ("Monte-Carlo samples","Polynomial Chaos"),'upper center' )
fig_file_name=sample_file_name+".pdf"
plt.savefig(fig_file_name)
