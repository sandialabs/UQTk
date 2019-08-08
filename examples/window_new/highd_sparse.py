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
from __future__ import print_function #so the print statements in python2 look right
try:
    from datetime import datetime
except ImportError:
    print("datetime could not be found")

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
    from tools_conductance_dp_pce_wrad import *
except ImportError:
    print("File with PCE utilities not found.")

try:
    import PyUQTk.PyPCE.pce_tool as pce_tools
except ImportError:
    print("PyUQTk pce_tools module not found")


start_time = datetime.now()
####################################################################
# The uncertain, Gaussian parameters used to calculate the heat flux

Ti = 293.0 # Room temperature (K)
To = 273.0 # Outside temperature (K)
hi = 2.0   # Coefficient of convection inside (W/m^2*s)
ho = 6.0   # Coefficient of convection outside (W/m^2*s)
dw = 0.005 # Width of the glass pane (m)
da = 0.01  # Width of the gap between the panes (m)
kw = 1.0   # Conduction constant for glass (W/m*K)
ka = 0.024 # Conduction constant for air (W/m*K)
Ts = 200.0   # Sky temperature (K)
mu = 1.86e-5 # Viscosity of air (kg/m*s)
rho = 1.29 # Density of air (kg/m^3)
beta = 3.67e-3 # Coefficient of thermal exspansion (1/K)
# Set the standard deviations of these uncertain parameters
std_Ti=Ti*0.005
std_To=To*0.01
std_dw=dw*0.01
std_da=da*0.01
std_kw=kw*0.05
std_ka=ka*0.05
std_hi=hi*0.15
std_ho=ho*0.15
std_Ts=Ts*0.1
std_mu=mu*0.05
std_rho=rho*0.05
std_beta=beta*0.05

# Number of random samples for KDE of distribution
n_kde = 100000
######### Forward Propagation using PCEs ##########
# Set verbose to 1 if you want intermediate print statements, otherwise set to 0
verbose = 0
nord = 3 # Order of the PCE
ndim = 12 # Number of dimensions of the PCE
pc_type = "HG"
pc_alpha = 0.0
pc_beta = 1.0
param= nord+1 # Parameter for quadrature point generation
# Equal to number of quad points per dimension for full quadrature or level for sparse quadrature

print("Performing forward UQ using sparse quadrature for a ", ndim, " dimensional PCE of order ",nord)

# Instantiate PC Object with sparse quadrature methods
pc_model2 = uqtkpce.PCSet("NISPnoq", nord,ndim,pc_type,pc_alpha,pc_beta)
pc_model2.SetQuadRule(pc_type,'sparse',param)
npce2 = pc_model2.GetNumberPCTerms() # Number of terms in the PCE
if verbose > 0:
    print("The number of terms in each PCE is ",npce2)
    pc_model2.PrintMultiIndexNormSquared()

print("\nInstantiation complete")

##### Sparse quadrature methods ######

#Get numpy array of quadrature points
qdpts2, totquat2= pce_tools.UQTkGetQuadPoints(pc_model2)

print("Sparse quadrature method requires model evaluations at ", totquat2, " points")
print("This may take a while if the chosen order is > 3")

# Convert Quadrature points in \xi_i to equivalent samples of input parameters
# (taking advantage of the fact that inputs are assumed to be Gaussian)
Ti_samples2 = Ti + std_Ti * qdpts2[:,0]
To_samples2 = To + std_To * qdpts2[:,1]
dw_samples2 = dw + std_dw * qdpts2[:,2]
da_samples2 = da + std_da * qdpts2[:,3]
kw_samples2 = kw + std_kw * qdpts2[:,4]
ka_samples2 = ka + std_ka * qdpts2[:,5]
hi_samples2 = hi + std_hi * qdpts2[:,6]
ho_samples2 = ho + std_ho * qdpts2[:,7]
Ts_samples2 = Ts + std_Ts * qdpts2[:,8]
mu_samples2 = mu + std_mu * qdpts2[:,9]
rho_samples2 = rho + std_rho * qdpts2[:,10]
beta_samples2 = beta + std_beta * qdpts2[:,11]

# Evaluate Forward model for sampled parameters
Q_evals2=fwd_model(Ti_samples2,To_samples2,Ts_samples2,dw_samples2,da_samples2,kw_samples2,\
ka_samples2,hi_samples2,ho_samples2,mu_samples2,rho_samples2,beta_samples2)
# Do the actual Galerkin Projection
c_k2 = pce_tools.UQTkGalerkinProjection(pc_model2,Q_evals2)
# Generate germ samples
germ_samples2=np.random.normal(0,1, (n_kde,ndim))
# Evaluate the PCE at the germ samples
pce_evals2=evaluate_pce(pc_model2,c_k2,germ_samples2)
# Peform kernel density estimation
xpts_pce2, PDF_data_pce2= KDE(pce_evals2)


# Plot the PDF curve
plt.figure(figsize=(10,10))
plt.plot(xpts_pce2, PDF_data_pce2, linewidth=2, color='g', label= 'NISP sparse quadrature method')
# Label Axes
plt.xlabel(r'Total Heat Flux ($W/m^2$)', size=16)
plt.ylabel("PDF", size=16)
# Add title
plt.suptitle("Heat Transfer Through a Dual Pane Window - PCE order %s"%(nord), size=20)
# Change tick size
plt.tick_params(axis='both', labelsize=14)
# Pad tick labels
plt.gca().tick_params(pad=6)
# Create legend
plt.legend(loc='upper left', prop={'size':12})
# Save figure
fig_name="td_heat_flux_pce.pdf"
plt.savefig(fig_name)
print("\ntd_heat_flux_pce.pdf has been saved.")
end_time = datetime.now()
print('Duration: {}'.format(end_time - start_time))
# Show figure
plt.show()
