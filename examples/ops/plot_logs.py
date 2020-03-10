#!/usr/bin/env python
#=====================================================================================
#
#                      The UQ Toolkit (UQTk) version 3.1.0
#                          Copyright (2020) NTESS
#                        https://www.sandia.gov/UQToolkit/
#                        https://github.com/sandialabs/UQTk
#
#     Copyright 2020 National Technology & Engineering Solutions of Sandia, LLC (NTESS).
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
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from   matplotlib.lines import Line2D
from   pylab import *
import scipy as spy
from   scipy import stats

#if ( len(sys.argv) > 1 ):
#    sample_file_name=sys.argv[1]
#else:
#    print "Need file name as argument"
#    quit()


#
# Load data for a, loga_int and loga_tay; set up KDE
#
a_samples=np.genfromtxt("samples.a.dat")
xa=np.linspace(a_samples.min(),a_samples.max(),200)
kernsa=stats.kde.gaussian_kde(a_samples)
pdfa=kernsa(xa)

lai_samples=np.genfromtxt("samples.loga_int.dat")
xlai=np.linspace(lai_samples.min(),lai_samples.max(),200)
kernslai=stats.kde.gaussian_kde(lai_samples)
pdflai=kernslai(xlai)

lat_samples=np.genfromtxt("samples.loga_tay.dat")
xlat=np.linspace(lat_samples.min(),lat_samples.max(),200)
kernslat=stats.kde.gaussian_kde(lat_samples)
pdflat=kernslat(xlat)

#
#  Make the figure
#
lw1=2
fs1=18

# Plot pdf of a
fig = plt.figure(figsize=(8,6))
ax=fig.add_axes([0.10, 0.10, 0.85, 0.85]) ;
l1=plt.plot(xa,pdfa,linewidth=lw1,label="a")

ax.set_xlabel("$a$",fontsize=fs1)
ax.set_ylabel("$PDF(a)$",fontsize=fs1)
plt.legend(loc="upper left")
fig_file_name="pdf_a.pdf"
plt.savefig(fig_file_name)

# Plot pdf of ln(a)
fig2 = plt.figure(figsize=(8,6))
ax=fig2.add_axes([0.10, 0.10, 0.85, 0.85]) ;
l1=plt.plot(xlai,pdflai,linewidth=lw1,label="ln(a), Int")
l2=plt.plot(xlat,pdflat,linewidth=lw1,label="ln(a), Tay")

ax.set_xlabel("$ln(a)$",fontsize=fs1)
ax.set_ylabel("$PDF(ln(a))$",fontsize=fs1)
plt.legend(loc="upper left")
fig_file_name="pdf_lna.pdf"
plt.savefig(fig_file_name)
