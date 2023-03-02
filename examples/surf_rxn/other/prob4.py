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
import shutil
import sys
import numpy as np
import math
import random as rnd
from scipy.stats.mstats import mquantiles

import fileinput
import file_utils

# define uqtkbin
if os.environ.get("UQTK_SRC") is None:
    print "Error: Need to set path to uqtk src as environment variable UQTK_SRC -> Abort"
    quit()
else:
    if ( not os.path.isdir(os.environ["UQTK_SRC"]) ):
        print "\"",os.environ["UQTK_SRC"],"\" is not a valid path -> Abort"
        quit()

uqtkbin=os.environ["UQTK_SRC"]+"/src_cpp/bin"
pcequad=uqtkbin+"/pce_quad"


# Input settings 
npt=100                    # Number of data points
chdimlist=[1,2]
chstart=np.array([1.6, 20.75])
pcord=3                    # PC order of the chain samples 
pctype='HG'                # PC type
chainfile='chain.dat'      # Name of the chain file
bw=0.001                      # Bandwidth for Rosenblatt (if negative, the program picks the optimal value) 
n_burnin=5000              # Burnin samples

## Prepare the xml file
paramxml_file='surf_rxn.par.xml'
shutil.copyfile('prob4_surf_rxn.in.xml.templ',paramxml_file)
for line in fileinput.input(paramxml_file, inplace = 1):
    print line.replace('PCTYPE', pctype),
for line in fileinput.input(paramxml_file, inplace = 1):
    print line.replace('PCORDER', str(pcord)),
for line in fileinput.input(paramxml_file, inplace = 1):
    print line.replace('CHAINFILE', chainfile),
    

## Synthetic data generation

# Prepare a single run with default parameters below
params=np.array([1.6, 20.75, 0.04, 1.0, 0.36, 0.016])
ndim=params.shape[0]
for idim in range(ndim,0,-1):
    if (idim not in chdimlist):
        for line in fileinput.input(paramxml_file, inplace = 1):
            print line.replace('PAR_'+str(idim), str(params[idim-1])),

shutil.copyfile(paramxml_file,"surf_rxn.in.xml")
for idim in chdimlist:
    for line in fileinput.input("surf_rxn.in.xml", inplace = 1):
        print line.replace('PAR_'+str(idim), str(params[idim-1])),



        
print "Generating synthetic data"
os.system('./SurfRxnDet.x > det.log')
os.system('cp solution.dat solution_ref.dat')
os.system('tail -n1 solution.dat > ss.dat')
ss_sol=np.loadtxt("ss.dat")
# Corrupt with multiplicative noise
stdfrac=0.1
noisy_data = np.array([np.random.normal(ss_sol[1],ss_sol[1]*stdfrac, npt), np.random.normal(ss_sol[2],ss_sol[2]*stdfrac, npt)])
#noisy_data = np.array([np.random.normal(ss_sol[1],0.02, npt), np.random.normal(ss_sol[2],0.02, npt)])
np.savetxt("inputdata.dat",np.transpose(noisy_data))

# Run the inference code
ii=0
shutil.copyfile(paramxml_file,"surf_rxn.in.xml")
for idim in chdimlist:
    for line in fileinput.input("surf_rxn.in.xml", inplace = 1):
        print line.replace('PAR_'+str(idim), str(chstart[ii])),
    ii=ii+1

print "Running the parameter inference"
os.system('./SurfRxnInfer.x')

#
# Import data from MCMC file
#
# load chain file
# it expects first column is line id, then the last
# two columns are alpha and current log posterior 
print "Loading in chain file",chainfile
all_samples, vnames = file_utils.extract_all_vars(chainfile,n_burnin,0,1)
n_all_vars = len(vnames)
n_cols = len(all_samples[0,:])
# Extract all MCMC chain variables in separate array
d0 = all_samples[:,1:1+n_all_vars]
samfile="insamples.dat"
np.savetxt(samfile,d0)

# Find PC coefficients corresponding to the chain samples
print "Running KDE-Rosenblatt transformation to build input PCE"
os.system(pcequad+' -o' + str(pcord)+' -f ' + samfile + ' -x ' + pctype + ' -w' + str(bw) + ' > pcequad.log')

# Prepare the xml file for forward propagation
for line in fileinput.input("surf_rxn.in.xml", inplace = 1):
    print line.replace('uncertain', 'det'),
for line in fileinput.input("surf_rxn.in.xml", inplace = 1):
    print line.replace('infer','uncertain'),

# Forward propagation
print "Running NISP forward propagation"
os.system('./SurfRxnNISP.x > nisp.log')
os.system('cp solution.dat solution_NISP.dat')


# Postprocessing
dirname='prob4'
# This needs to be fixed as it crashes the script if the directory prob4 does not exist
shutil.rmtree(dirname)
os.mkdir(dirname,0755)

print "Generating plots in directory "+dirname

# Plot data and the reference run for illustrations
os.system('./plot_uvdata.py')
shutil.copy('uvwt_data.eps',dirname)
shutil.copy('uv_data.eps',dirname)

# Plot the forward propagated time series, utility borrowed from prob 2
os.system('./plSurfRxnMstd.py NISP')
shutil.copyfile('prob2_surf_rxn_NISP_mstd.pdf',dirname+os.sep+'PCprop.pdf')

# Plot posterior predictive
# A horrible way of doing this, writing the dimensions to be inferred in a file
# for posterior plotting routines
with open("chdimlist",'w') as f:
    for i in range(len(chdimlist)):
        f.write(str(chdimlist[i])+' ')
os.system('./postpred.py '+' '.join(str(e) for e in chdimlist))
shutil.copy('postpred.eps',dirname)
shutil.copy('postpred_zoom.eps',dirname)

# Plot chains
os.system('rm chain_param_*.eps')
os.system('./plot_chain_samples.py')
os.system('cp chain_param_*.eps '+dirname)
