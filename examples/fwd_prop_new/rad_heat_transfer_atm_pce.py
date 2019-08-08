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

import sys

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
    from heat_transfer_pce_tools import *
except ImportError:
    print("File with PCE utilities not found.")

try:
    import PyUQTk.PyPCE.pce_tools as pce_tools
except ImportError:
    print("could not import pce_tools")
    sys.exit(1)

np.random.seed(123)
#####################################################################

#####################################################################
# Some general settings
main_verbose = 1
compute_rad = True # Whether or not radiation is considered
if main_verbose > 0:
    print("\nConsidering heat transfer with ", end= ' ')
    if not compute_rad:
        print("no ",end= ' ')
    print("radiation.")

nord = 3 # Order of the PCE
# Number of dimensions of the PCE (number of uncertain variables)
if compute_rad:
    ndim = 7
else:
    ndim = 6
pc_type = "HG" # Use Wiener-Hermite PCE since inputs have a normal distribution
pc_alpha = 0.0
pc_beta = 1.0
param= nord+1   # Parameter for quadrature point generation
                # Equal to number of quad points per dimension for full quadrature or level for sparse quadrature
n_MC = 100000   # Number of random samples to use in MC sampling (of the full problem or of the PCE of the solution)


# Nominal values of the parameters used to calculate the heat flux
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
# Generate germ samples
germ_samples=np.random.normal(0,1, (n_MC,ndim))
# Evaluate the PCE at the germ samples
pce_evals=pce_tools.UQTkEvalPC(pc_model,c_k,germ_samples[0,:])
#Peform kernel density estimation
print("shape pce_evals = ", np.shape(pce_evals))
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

# Generate germ samples for KDE
germ_samples2=np.random.normal(0,1, (n_MC,ndim))
# Evaluate the PCE at the germ samples
pce_evals2=pce_tools.UQTkEvalPC(pc_model2,c_k2,germ_samples2[0,:])
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
