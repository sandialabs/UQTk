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
#     Questions? Contact the UQTk Developers at <uqtk-developers@software.sandia.gov>
#     Sandia National Laboratories, Livermore, CA, USA
#=====================================================================================
from __future__ import print_function #so print statements in python2 look right
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
    import PyUQTk.utils as uqtkutils
except ImportError:
    print("PyUQTk tools module not found")
try:
    from fwd_prop_tools import *
except ImportError:
    print("File with PCE utilities not found.")

try:
    import PyUQTk.PyPCE.pce_tools as pce_tools
except ImportError:
    print("PyUQTk pce_tools module not found.")

start_time = datetime.now()
####################################################################
# The uncertain, Gaussian parameters used to calculate the heat flux


hi = 2.0   # Coefficient of convection inside (W/m^2*s)
ho = 6.0   # Coefficient of convection outside (W/m^2*s)
kw = 1.0   # Conduction constant for glass (W/m*K)
ka = 0.024 # Conduction constant for air (W/m*K)
Ts = 200.0   # Sky temperature (K)

# Set the standard deviations of these uncertain parameters
std_Ts=Ts*0.1
std_hi=hi*0.15
std_ho=ho*0.15
std_kw=kw*0.05
std_ka=ka*0.05

######### Forward Propagation using Monte Carlo sampling #########

# Number of random samples
n = 100000
# Generate random samples of the uncertain parameters
samp_Ts=np.random.normal(Ts, std_Ts, n)
samp_hi=np.random.normal(hi, std_hi, n)
samp_ho=np.random.normal(ho, std_ho, n)
samp_kw=np.random.normal(kw, std_kw, n)
samp_ka=np.random.normal(ka, std_ka, n)
# Create list of sample arrays
samples=[samp_Ts,samp_hi,samp_ho,samp_kw,samp_ka]

Q_samples=[]
# List to store values of Q assuming radiative heat transfer occurs
Q_r_samples=[]
# Calculate values of heat flux Q (assuming no radiative heat transfer)
# for the different sample values and append to the list
for i in range(n):
    (Q,T1,T2,T3,T4)=compute_heat_flux(samp_hi[i],samp_ho[i],samp_kw[i],samp_ka[i])

    Q_samples.append(Q)
# Calculate values of heat flux Q assuming radiative heat transfer to atmosphere and append to list
# For the required estimates of Q,T1,T2,T3 and T4 needed to solve the nonlinear system,
# we use the values obtained by solving the system assuming no radiative heat transfer
    Q2=r_heat_flux(samp_Ts[i],samp_hi[i],samp_ho[i],samp_kw[i],samp_ka[i],(Q,T1,T2,T3,T4))
    Q_r_samples.append(Q2)
# Convert Q_r_samples into numpy array
Q_r_samples=np.array(Q_r_samples)


# Perform KDE on Q_r_samples
xpts_Q2, PDF_data_Q2= KDE(Q_r_samples)
#############################################

######### Forward Propagation using PCEs ##########
# Set verbose to 1 if you want intermediate print statements, otherwise set to 0
verbose = 0
nord = 4 # Order of the PCE
ndim = 5 # Number of dimensions of the PCE
pc_type = "HG"
pc_alpha = 0.0
pc_beta = 1.0
param= nord+1 # Parameter for quadrature point generation
# Equal to number of quad points per dimension for full quadrature or level for sparse quadrature

# Instantiate both PC objects

# Instantiate PC Object with full quadrature methods
print("\nInstantiating PC Objects\n")

pc_model = uqtkpce.PCSet("NISPnoq", nord,ndim,pc_type, pc_alpha,pc_beta)
pc_model.SetQuadRule(pc_type, 'full', param)
npce = pc_model.GetNumberPCTerms() # Number of terms in the PCE
if verbose > 0:
    print("The number of terms in each PCE is",npce)
    pc_model.PrintMultiIndexNormSquared()

# Instantiate PC Object with sparse quadrature methods
pc_model2 = uqtkpce.PCSet("NISPnoq", nord,ndim,pc_type, pc_alpha,pc_beta)
pc_model2.SetQuadRule(pc_type, 'sparse', param)
npce2 = pc_model2.GetNumberPCTerms() # Number of terms in the PCE
if verbose > 0:
    print("The number of terms in each PCE is",npce2)
    pc_model2.PrintMultiIndexNormSquared()
print("\nInstantiation complete")

########## Full quadrature methods ###########

#Get numpy array of quadrature points
qdpts, totquat= pce_tools.UQTkGetQuadPoints(pc_model)
# Convert Quadrature points in \xi_i to equivalent samples of input parameters
# (taking advantage of the fact that inputs are assumed to be Gaussian)
Ts_samples = Ts + std_Ts * qdpts[:,0]
hi_samples = hi + std_hi * qdpts[:,1]
ho_samples = ho + std_ho * qdpts[:,2]
kw_samples = kw + std_kw * qdpts[:,3]
ka_samples = ka + std_ka * qdpts[:,4]

# Evaluate Forward model for sampled parameters
Q_evals=fwd_model(Ts_samples,hi_samples,ho_samples,kw_samples,ka_samples,)
# Do the actual Galerkin Projection
c_k = pce_tools.UQTkGalerkinProjection(pc_model,Q_evals)
# Generate germ samples
germ_samples=np.random.normal(0,1, (n,ndim))
# Evaluate the PCE at the germ samples
pce_evals=evaluate_pce(pc_model,c_k,germ_samples)

#Peform kernel density estimation
xpts_pce, PDF_data_pce = KDE(pce_evals)

##### Sparse quadrature methods ######

#Get numpy array of quadrature points
qdpts2, totquat2= pce_tools.UQTkGetQuadPoints(pc_model2)

# Convert Quadrature points in \xi_i to equivalent samples of input parameters
# (taking advantage of the fact that inputs are assumed to be Gaussian)
Ts_samples2 = Ts + std_Ts * qdpts2[:,0]
hi_samples2 = hi + std_hi * qdpts2[:,1]
ho_samples2 = ho + std_ho * qdpts2[:,2]
kw_samples2 = kw + std_kw * qdpts2[:,3]
ka_samples2 = ka + std_ka * qdpts2[:,4]

# Evaluate Forward model for sampled parameters
Q_evals2=fwd_model(Ts_samples2,hi_samples2,ho_samples2,kw_samples2,ka_samples2,)
# Do the actual Galerkin Projection
c_k2 = pce_tools.UQTkGalerkinProjection(pc_model2,Q_evals2)
# Generate germ samples
germ_samples2=np.random.normal(0,1, (n,ndim))
# Evaluate the PCE at the germ samples
pce_evals2=evaluate_pce(pc_model2,c_k2,germ_samples2)
#Peform kernel density estimation
xpts_pce2, PDF_data_pce2= KDE(pce_evals2)
# Print statements to indicate number of samples used
print("\nMonte Carlo sampling used ", n, " points")
print("Full quadrature method used ", totquat, " points")
print("Sparse quadrature method used ", totquat2, " points")
# Produce the inverse relationship graph
plot_mi_dims(pc_model,c_k,ndim)
# Plot the three PDF curves on the same figure
plt.figure(figsize=(10,10))
plt.plot(xpts_pce, PDF_data_pce, linewidth=2, color='r', label='NISP full quadrature method')
plt.plot(xpts_Q2, PDF_data_Q2, linewidth=2, color='b', label='Monte Carlo Sampling')
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
fig_name="fd_heat_flux_pce.pdf"
plt.savefig(fig_name)
print("\nfd_heat_flux_pce.pdf has been saved.")

end_time = datetime.now()
print('Duration: {}'.format(end_time - start_time))

# Show figure
plt.show()
