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
  tmcmc_bimodal.py [-h] [--noplots]
what
  Use TMCMC to sample from a 3-dimensional posterior that's a product of a Gaussian prior and a bimodal likelihood;
  generate plots and statistics of the output
where
  -h = print help info
  --noplots = do not generate chain and posterior plots [defaults to False]
"""

#
# Process inputs
#
try:
    # opts,v_names = getopt.getopt(sys.argv[1:],"hi:s:",["nb="])
    opts,extra_arguments = getopt.getopt(sys.argv[1:],"h",["noplots"])
except getopt.GetoptError as err:
    print(str(err))
    print(help_string)
    sys.exit(1)

# Default values
generate_plots = True # Generate chain and posterior pdf plots

print(opts)

for o,a in opts:
    if o == "-h":
        print(help_string)
        sys.exit(0)
    elif o == "--noplots":
        generate_plots = False
    else:
        assert False, "unhandled command line parsing option"


# Other input settings
dense_plots = True      # Set to True to get dense format for posterior plots in triangular format
chainfile       = 'tmcmc_chain.dat'      # Name of the chain file


## Run the inference code #####################################################

print("Running TMCMC on bimodal posterior")
os.system('./tmcmc_umbridge.x')# > logTMCMC.dat')

## Import data from MCMC file ###################################################
print("Loading in chain file ",chainfile )
all_samples, vnames = postp.extract_all_vars(chainfile,0,0,1,False)
n_all_vars = len(vnames)
print(all_samples.shape)

# Extract all MCMC chain variables in separate array
chn  = all_samples[:,0:1+n_all_vars]
nchn = chn.shape[0]

if generate_plots:

    ## Scatter plots from intermediate pdfs ##################################################################

    indx = 0
    while 1:
        fname = './TMCMCIntermediates/samples.dat.'+str(indx)
        if os.path.exists(fname):
            samples_file = open(fname, 'r')

            # Extract first line to see how many columns we have
            first_line = samples_file.readline().rstrip('\n')
            first_line_items = [item for item in first_line.split()]
            n_cols = len(first_line_items)
            samples_file.seek(0)

            samps = []
            line_no = 0
            done = 0
            while not done:
                line = samples_file.readline()
                if (line == ""):
                    done = 1
                else:
                    line_no += 1
                    records = line.split()
                    num_records = [float(s) for s in records]
                    samps.append(num_records)
            samps = np.array(samps)

            for i in range(n_all_vars):
                for j in range(i):
                    fig = plt.figure(figsize=(10,10))
                    ax=fig.add_axes([0.12,0.12,0.8,0.8])
                    plt.plot(samps[:,j],samps[:,i],'o',markeredgecolor='blue',markerfacecolor='b',markersize=5)
                    plt.xlim(-4,4)
                    plt.ylim(-4,4)
                    ax.set_xlabel(vnames[j],fontsize=22)
                    ax.set_ylabel(vnames[i],fontsize=22)
                    plt.savefig('tmcmc_bimodal.chn_'+vnames[j]+'_'+vnames[i]+'.intermediate'+str(indx)+'.pdf')
                    plt.clf()
            indx = indx+10
        else:
            break

    ## Scatter plots from posterior ##################################################################
    for i in range(n_all_vars):
        for j in range(i):
            fig = plt.figure(figsize=(10,10))
            ax=fig.add_axes([0.12,0.12,0.8,0.8])
            plt.plot(chn[:,j+1],chn[:,i+1],'o',markeredgecolor='blue',markerfacecolor='b',markersize=5)
            ax.set_xlabel(vnames[j],fontsize=22)
            ax.set_ylabel(vnames[i],fontsize=22)
            plt.savefig('tmcmc_bimodal.chn_'+vnames[j]+'_'+vnames[i]+'.posterior.pdf')
            plt.clf()

    ## Plot posterior PDF 'triangle' ################################################
    np_kde = 100
    postp.plot_all_posteriors(chn[:,1:],vnames,np_kde,"tmcmc_bimodal.posteriors",0,dense_plots)

print("END: TMCMC example done: check out tmcmc_bimodal.*.pdf files for results.")
