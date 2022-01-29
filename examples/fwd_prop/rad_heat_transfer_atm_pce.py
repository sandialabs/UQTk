#!/usr/bin/env python
#=====================================================================================
#
#                      The UQ Toolkit (UQTk) version 3.1.2
#                          Copyright (2022) NTESS
#                        https://www.sandia.gov/UQToolkit/
#                        https://github.com/sandialabs/UQTk
#
#     Copyright 2022 National Technology & Engineering Solutions of Sandia, LLC (NTESS).
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
from __future__ import print_function #so the print statements in python2 looks right

import sys
import argparse

try:
	import numpy as np
except ImportError:
	print("Numpy module could not be found")

try:
	import matplotlib.pyplot as plt
	from matplotlib import rc
	rc('mathtext', default='regular')
except ImportError:
	print("Matplotlib not found")

try:
	from scipy import stats
except ImportError:
	print("Scipy stats module could not be found")

import sys

try:
    import PyUQTk.uqtkarray as uqtkarray
except ImportError:
    print("PyUQTk array module not found")
try:
    import PyUQTk.quad as uqtkquad
except ImportError:
    print("PyUQTk quad module not found")
try:
    import PyUQTk.pce as uqtkpce
except ImportError:
    print("PyUQTk PCE module not found")
try:
    import PyUQTk.tools as uqtktools
except ImportError:
    print("PyUQTk tools module not found")

try:
    import PyUQTk.PyPCE.pce_tools as pce_tools
except ImportError:
    print("PyUQTk pce_tools module not found")

try:
    from heat_transfer_pce_tools import *
except ImportError:
    print("File with PCE utilities not found.")


#####################################################################

#####################################################################
# Parse input arguments
usage_str = """Script produce a graph comparing PDFs of heat flux generated using
	NISP full and sparse quadrature methods and a Monte Carlo sampling method
	Look at manual example Forward Propagation of Uncertainty with PyUQTk for more explanation"""
parser = argparse.ArgumentParser(description=usage_str)
parser.add_argument("--no_verbose", help="To trun off intermediate print statements.", action='store_false')
parser.add_argument("-r", "--no_compute_rad", help="To not compute radiation. [Default to include radiation]", action='store_false')
parser.add_argument("--nord", dest="nord", type=int,
                    default=3, help="the order of the PCE [default = 3]")
parser.add_argument("--ndim", dest="ndim", type=int,
                    default=None, help="the number of dimensions of the PCE [default: 7 with radiation, 6 without radiation]")
parser.add_argument("--pc_type", dest="pc_type", type=str,
                    default="HG", help="indicates the polynomial type and weighting function [default = HG]")
parser.add_argument("-a", dest="pc_alpha", type=float,
                    default=0.0, help="Free parameter greater than -1. Used with Gamma-Laguerre and Beta-Jacobi PC types. [default = 0.0]")
parser.add_argument("-b", dest="pc_beta", type=float,
                    default=1.0, help="Free parameter greater than -1. Used with Gamma-Laguerre and Beta-Jacobi PC types. [default = 1.0]")
parser.add_argument("--param", dest="param", type=int, default=None,
					help="The parameter used for quadrature point generation. Equal to the number of quadrature points per dimension for full quadrature or the level for sparse quadrature methods. This parameter is generally set to nord + 1 in order to have the right polynomial exactness. [Default is nord + 1]")
parser.add_argument("--n_MC", dest="n_MC", type=int, default=100000,
					help = "Number of random samples to use in MC sampling (of the full problem or of the PCE of the solution) [default = 100000]")
args = parser.parse_args()

#save parsed arguments in variables
main_verbose = args.no_verbose
compute_rad = args.no_compute_rad
nord = args.nord
ndim = args.ndim
pc_type = args.pc_type
pc_alpha = args.pc_alpha
pc_beta = args.pc_beta
param = args.param
n_MC = args.n_MC

#set param if not set in command line
if(not param):
	param = nord + 1

if main_verbose > 0:
    print("\nConsidering heat transfer with ", end= ' ')
    if not compute_rad:
        print("no ",end= ' ')
    print("radiation.")

# Number of dimensions of the PCE (number of uncertain variables)
if(not ndim):
	if compute_rad:
	    ndim = 7
	else:
	    ndim = 6


# Nominal values of the parameters used to calculate the heat flux
#These values are the same as shown in tha table in example forward propgation of uncertainty with PyUQTk in the manual
Ti = 293.0 # Room temperature in K
To = 273.0 # Outside temperature in K
dw = 0.01  # Window thickness in m
kw = 1.0   # Window conductivity in W/mK
hi = 2.0   # Inner convective heat transfer coefficient in W/m^2K
ho = 6.0   # Outer wall convective heat transfer coefficient in W/m^2K
TA = 150.0 # Atmospheric temperature in K. Only used to calculate heat flux
		   # when it is assumed that radiative heat transfer occurs

# Set the standard deviations of the Gaussian uncertainty of these uncertain parameters
std_Ti=Ti*0.005
std_To=To*0.005
std_dw=dw*0.01
std_kw=kw*0.05
std_hi=hi*0.15
std_ho=ho*0.15
std_TA=TA*0.1

#
######### Forward Propagation using Monte Carlo sampling #########
#

# Generate random samples of the uncertain parameters
samp_Ti=np.random.normal(Ti, std_Ti, n_MC)
samp_To=np.random.normal(To, std_To, n_MC)
samp_dw=np.random.normal(dw, std_dw, n_MC)
samp_kw=np.random.normal(kw, std_kw, n_MC)
samp_hi=np.random.normal(hi, std_hi, n_MC)
samp_ho=np.random.normal(ho, std_ho, n_MC)
if compute_rad:
    samp_TA=np.random.normal(TA, std_TA, n_MC)

# List to store values of Q calculated from the random samples of the parameters
Q_samples=[]

# Evaluate heat transfer model for all parameter samples
if not compute_rad:
    Q_evals=fwd_model(samp_Ti,samp_To,samp_dw,samp_kw,samp_hi,samp_ho,main_verbose)
else:
    Q_evals=fwd_model_rad(samp_Ti,samp_To,samp_dw,samp_kw,samp_hi,samp_ho,samp_TA,main_verbose)

# Perform KDE on Q_samples
xpts_MC, PDF_data_MC= KDE(Q_evals)

#
######### Forward Propagation using PCEs and full quadrature ##########
#

# Instantiate PC Object with full quadrature methods
print("\nInstantiating PC Object\n")

pc_model = uqtkpce.PCSet("NISP", nord,ndim,pc_type, pc_alpha,pc_beta)
pc_model.SetQuadRule(pc_type, 'full', param)
npce = pc_model.GetNumberPCTerms() # Number of terms in the PCE
if main_verbose > 0:
	print("The number of terms in each PCE is",npce)
if main_verbose > 10:
	pc_model.PrintMultiIndexNormSquared()

print("\nInstantiation complete")

#Get numpy array of quadrature points
qdpts, totquat = pce_tools.UQTkGetQuadPoints(pc_model)

# Convert Quadrature points in \xi_i to equivalent samples of input parameters
# (taking advantage of the fact that inputs are assumed to be Gaussian)
# This is equivalent to evaluating 1st order PC expansions for the input parameters.
Ti_samples = Ti + std_Ti * qdpts[:,0]
To_samples = To + std_To * qdpts[:,1]
dw_samples= dw + std_dw * qdpts[:,2]
kw_samples= kw + std_kw * qdpts[:,3]
hi_samples= hi + std_hi * qdpts[:,4]
ho_samples= ho + std_ho * qdpts[:,5]
if compute_rad:
    TA_samples= TA + std_TA * qdpts[:,6]

# Evaluate Forward model for sampled parameters
if not compute_rad:
    Q_evals=fwd_model(Ti_samples,To_samples, dw_samples, kw_samples,hi_samples,ho_samples,main_verbose)
else:
    Q_evals=fwd_model_rad(Ti_samples,To_samples, dw_samples, kw_samples,hi_samples,ho_samples,TA_samples,main_verbose)

# Do the actual Galerkin Projection
c_k = pce_tools.UQTkGalerkinProjection(pc_model,Q_evals)

pc_model.SeedBasisRandNumGen(123)
#Draw samples of PCE evaulations
pce_evals = pce_tools.UQTkDrawSamplesPCE(pc_model, c_k, n_MC)

#Peform kernel density estimation
xpts_pce, PDF_data_pce= KDE(pce_evals)

#
##### Forward Propagation using PCEs and sparse quadrature ######
#

# Instantiate PC Object with sparse quadrature methods
print("\nInstantiating PC Object\n")

pc_model2 = uqtkpce.PCSet("NISP", nord,ndim,pc_type, pc_alpha,pc_beta)
pc_model2.SetQuadRule(pc_type, 'sparse', param)
npce2 = pc_model2.GetNumberPCTerms() # Number of terms in the PCE
if main_verbose > 0:
	print("The number of terms in each PCE is",npce2)
if main_verbose > 10:
	pc_model.PrintMultiIndexNormSquared()

print("\nInstantiation complete")

#Get numpy array of quadrature points
qdpts2, totquat2= pce_tools.UQTkGetQuadPoints(pc_model2)

# Convert Quadrature points in \xi_i to equivalent samples of input parameters
# (taking advantage of the fact that inputs are assumed to be Gaussian)
# This is equivalent to evaluating 1st order PC expansions for the input parameters.
Ti_samples2 = Ti + std_Ti * qdpts2[:,0]
To_samples2 = To + std_To * qdpts2[:,1]
dw_samples2= dw + std_dw * qdpts2[:,2]
kw_samples2= kw + std_kw * qdpts2[:,3]
hi_samples2= hi + std_hi * qdpts2[:,4]
ho_samples2= ho + std_ho * qdpts2[:,5]
if compute_rad:
    TA_samples2= TA + std_TA * qdpts2[:,6]

# Evaluate Forward model for sampled parameters
if not compute_rad:
    Q_evals2=fwd_model(Ti_samples2,To_samples2, dw_samples2, kw_samples2,hi_samples2,ho_samples2,main_verbose)
else:
    Q_evals2=fwd_model_rad(Ti_samples2,To_samples2, dw_samples2, kw_samples2,hi_samples2,ho_samples2,TA_samples2,main_verbose)

# Do the actual Galerkin Projection
c_k2 = pce_tools.UQTkGalerkinProjection(pc_model2,Q_evals2)


#Draw samples of PCE evaulations
pce_evals2 = pce_tools.UQTkDrawSamplesPCE(pc_model2, c_k2, n_MC)

#Peform kernel density estimation
xpts_pce2, PDF_data_pce2= KDE(pce_evals2)

#
# Display and plot summary information
#

# print(statements to indicate number of samples used
print("\nMonte Carlo sampling used %s points" %(n_MC))
print("Full quadrature method used %s points"%(totquat))
print("Sparse quadrature method used %s points"%(totquat2))

# Plot the three PDF curves on the same figure
plt.figure(figsize=(10,10))
plt.plot(xpts_pce, PDF_data_pce, linewidth=2, color='r', label='NISP full quadrature method')
plt.plot(xpts_MC, PDF_data_MC, linewidth=2, color='b', label='Monte Carlo Sampling')
plt.plot(xpts_pce2, PDF_data_pce2, linewidth=2, color='g', label= 'NISP sparse quadrature method')

# Label Axes
plt.xlabel("Total Heat Flux ($W/m^2$)", size=16)
plt.ylabel("PDF", size=16)
# Add title
plt.suptitle("Heat Transfer Through a Window", size=20)
# Change tick size
plt.tick_params(axis='both', labelsize=14)
# Pad tick labels
plt.gca().tick_params(pad=6)
# Create legend
plt.legend(loc='upper left', prop={'size':12})
# Save figure
fig_name="heat_flux_pce.pdf"
plt.savefig(fig_name)
print("\nheat_flux_pce.pdf has been saved.")
# Show figure
plt.show()
