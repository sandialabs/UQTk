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

from __future__ import print_function #so the print statements in python2 looks right
import os
import shutil
import sys
import math
import getopt

try:
    import numpy as np
except ImportError:
    print("\ntmcmc_bimodal.py requires numpy package -> Exit\n")
    quit()

import random as rnd
from scipy.stats.mstats import mquantiles
import matplotlib as mpl
import matplotlib.pyplot as plt

import fileinput

try:
    import PyUQTk.inference.postproc as postp
except ImportError:
    print("\ntmcmc_bimodal.py requires the pyUQTk package.",
        "Make sure the UQTk root directory is included in the PYTHONPATH environment variable.\n")
    quit()

#from pylab import *
mpl.rcParams['legend.loc'     ] = 'upper right'
mpl.rcParams['legend.fontsize'] = 20
mpl.rcParams['lines.linewidth'] = 4
mpl.rcParams['lines.color'    ] = 'r'
mpl.rcParams['axes.linewidth' ] = 3
mpl.rcParams['axes.linewidth' ] = 3
mpl.rcParams['axes.grid'      ] = True
mpl.rcParams['axes.labelsize' ] = 22
mpl.rcParams['xtick.labelsize'] = 16
mpl.rcParams['ytick.labelsize'] = 16

###########################################################################
uqtkbin="../../bin"
pcequad=uqtkbin+"/pce_quad"
###########################################################################

help_string = """
Usage:
  tmcmc_bimodal.py [-h] [--noplots] [--stats]
what
  Use TMCMC to sample from a 3-dimensional posterior that's a product of a Gaussian prior and a bimodal likelihood;
  generate plots and statistics of the output
where
  -h = print help info
  --noplots = do not generate chain and posterior plots [defaults to False]
  --stats = generate MCMC chain statistics [defaults to False]
"""

#
# Process inputs
#
try:
    # opts,v_names = getopt.getopt(sys.argv[1:],"hi:s:",["nb="])
    opts,extra_arguments = getopt.getopt(sys.argv[1:],"h",["noplots","stats"])
except getopt.GetoptError as err:
    print(str(err))
    print(help_string)
    sys.exit(1)

# Default values
generate_plots = True # Generate chain and posterior pdf plots
generate_stats = False # By default, do not generate MCMC statistics

print(opts)

for o,a in opts:
    if o == "-h":
        print(help_string)
        sys.exit(0)
    elif o == "--noplots":
        generate_plots = False
    elif o == "--stats":
        generate_stats = True
    else:
        assert False, "unhandled command line parsing option"


# Other input settings
dense_plots = True      # Set to True to get dense format for posterior plots in triangular format
chainfile       = 'tmcmc_chain.dat'      # Name of the chain file


## Run the inference code #####################################################

print("Running TMCMC on bimodal posterior")
os.system('./tmcmc_bimodal.x >& logTMCMC.dat')

## Import data from MCMC file ###################################################
print("Loading in chain file ",chainfile )
all_samples, vnames = postp.extract_all_vars(chainfile,0,0,1,False)
n_all_vars = len(vnames)
print(all_samples.shape)

# Extract all MCMC chain variables in separate array
chn  = all_samples[:,0:1+n_all_vars]
nchn = chn.shape[0]

if generate_plots:
    ## Scatter plots from posterior ##################################################################
    for i in range(n_all_vars):
        for j in range(i):
            fig = plt.figure(figsize=(10,10))
            ax=fig.add_axes([0.12,0.12,0.8,0.8])
            plt.plot(chn[:,j+1],chn[:,i+1],'o',markeredgecolor='blue',markerfacecolor='b',markersize=5)
            ax.set_xlabel(vnames[j],fontsize=22)
            ax.set_ylabel(vnames[i],fontsize=22)
            plt.savefig('tmcmc_bimodal.chn_'+vnames[j]+'_'+vnames[i]+'.pdf')
            plt.clf()

    ## Plot posterior PDF 'triangle' ################################################
    np_kde = 100
    postp.plot_all_posteriors(chn[:,1:],vnames,np_kde,"tmcmc_bimodal.posteriors",0,dense_plots)

if generate_stats:
    postp.get_mcmc_stats(all_samples,vnames,"tmcmc_bimodal",0)

print("END: TMCMC example done: check out tmcmc_bimodal.*.pdf files for results.")
