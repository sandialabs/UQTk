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
import sys
import math

try:
    import numpy as np
except ImportError:
    print("Numpy module not found")

try:
    import matplotlib.pyplot as plt
except ImportError:
    print("Matplotlib not found")

try:
    import d_springs_tools
except ImportError:
    print("File with d_springs_tools not found")

sys.path.append("../../PyUQTk/adaptation_tools")
sys.path.append("../../PyUQTk/pce_tools")
sys.path.append("../../PyUQTk/pyuqtkarray")
sys.path.append("../../PyUQTk/pyuqtkarray_tools")
sys.path.append("../../PyUQTk/quad")
sys.path.append("../../PyUQTk/pce")
sys.path.append("../../PyUQTk/tools")
sys.path.append("../../PyUQTk/")

#try:
import PyUQTk.PyPCE.adaptation_tools as adp
from PyUQTk.PyPCE import pce_tools
#except ImportError:
    #print("PyUQTk pce_tools or adaptation_tools module not found")



try:
    import _uqtkarray as uqtkarray
    import _quad as uqtkquad
    import _pce as uqtkpce
    import _tools as uqtktools
except ImportError:
    print("PyUQTk array, quad, PCE, or tools module not found")
####################################################

def forward_propagation(x1, x2, x3, x4, x5, x6, x7, \
    std_x1, std_x2, std_x3, std_x4, std_x5, std_x6, std_x7, a, b, \
    nord, ndim, pc_type,param, R, main_verbose=0, sf="sparse", pc_alpha=0.0, pc_beta=1.0):
    # Obtain PC model
    pc_model = uqtkpce.PCSet("NISP", nord,ndim,pc_type, pc_alpha,pc_beta)
    # Set quadrature rule and obtain quadrature points
    pc_model.SetQuadRule(pc_type, sf, param)
    npce = pc_model.GetNumberPCTerms() # Number of terms in the PCE
    qdpts, totquat= pce_tools.UQTkGetQuadPoints(pc_model)
    # map quadrature from eta space to xi space
    qdpts_xi = adp.eta_to_xi_mapping(qdpts, R)
    # map quadrature from xi space to input parameters
    xx = np.zeros((np.shape(qdpts_xi)[0],np.shape(qdpts_xi)[1]))
    xx[:,0] = x1 + qdpts_xi[:,0]*std_x1
    xx[:,1] = x2 + qdpts_xi[:,1]*std_x2
    xx[:,2] = x3 + qdpts_xi[:,2]*std_x3
    xx[:,3] = x4 + qdpts_xi[:,3]*std_x4
    xx[:,4] = x5 + qdpts_xi[:,4]*std_x5
    xx[:,5] = x6 + qdpts_xi[:,5]*std_x6
    xx[:,6] = x7 + qdpts_xi[:,6]*std_x7
    # evaluate QoI of the input parameters
    Q_evals = d_springs_tools.fwd_model(xx, a, b, main_verbose)
    c_k = pce_tools.UQTkGalerkinProjection(pc_model,Q_evals)
    return pc_model, c_k, totquat

##########################################################
main_verbose = 1
nord = 3   # order of PCE
nord0 = 1  # order to obtain Gaussian coefficients
ndim = 7   # dimension
pc_type = "HG"  # Polynomial type
pc_alpha = 0.0
pc_beta = 1.0
param = nord+1  # quadrature level
param0 = 1      # quadrature level for first order expansion
method = 3      # sort by importance method to obtain rotation matrix
sf = "sparse"   # full or sparse quadrature
n_MC = 100000   # # of MC samples

a = 1.0         # parameters in the forward model
b = 0.5         # parameters in the forward model

mu1 = 5.0    # mean value for x1 to x4
mu2 = 4.0    # mean value for x5 to x7
x1 = mu1     # Gaussian
x2 = mu1     # Gaussian
x3 = mu1     # Gaussian
x4 = mu1     # Gaussian
x5 = mu2     # Gaussian
x6 = mu2     # Gaussian
x7 = mu2     # Gaussian

sigma1 = 0.6  # standard deviation for x1 to x4
sigma2 = 0.5  # standard deviation for x1 to x4
std_x1 = sigma1
std_x2 = sigma1
std_x3 = sigma1
std_x4 = sigma1
std_x5 = sigma2
std_x6 = sigma2
std_x7 = sigma2

##################################################################
######### Forward Propagation using Monte Carlo sampling #########
##################################################################
if main_verbose > 0:
    print("\nStaring Monte Carlo Sampling method")

xx = np.zeros((n_MC,ndim))
xx[:,0] = np.random.normal(x1, std_x1, n_MC)
xx[:,1] = np.random.normal(x2, std_x2, n_MC)
xx[:,2] = np.random.normal(x3, std_x3, n_MC)
xx[:,3] = np.random.normal(x4, std_x4, n_MC)
xx[:,4] = np.random.normal(x5, std_x5, n_MC)
xx[:,5] = np.random.normal(x6, std_x6, n_MC)
xx[:,6] = np.random.normal(x7, std_x7, n_MC)
Q_evals = d_springs_tools.fwd_model(xx, a, b, main_verbose)
xpts_MC, PDF_data_MC= d_springs_tools.KDE(Q_evals)


##################################################################
##### Forward Propagation using PCEs and sparse quadrature #######
##################################################################
if main_verbose > 0:
    print("\nStaring NISP sparse quadrature method")

R0 = np.eye(ndim)
pc_model2, c_k2, totquat2 = forward_propagation(x1, x2, x3, x4, x5, x6, x7, \
    std_x1, std_x2, std_x3, std_x4, std_x5, std_x6, std_x7, a, b, \
    nord, ndim, pc_type, param, R0, main_verbose, sf="sparse", pc_alpha=0.0, pc_beta=1.0)

# Generate germ samples
germ_samples2=np.random.normal(0,1, (n_MC,ndim))
# Evaluate the PCE at the germ samples
pce_evals2=d_springs_tools.EvaluatePCE(pc_model2,c_k2,germ_samples2)
#Peform kernel density estimation
xpts_pce2, PDF_data_pce2= d_springs_tools.KDE(pce_evals2)

###################################################################
##### Forward Propagation using Adaptative PCEs ###################
###################################################################
if main_verbose > 0:
    print("\n###################################")
    print("# Staring basis adaptation method #")
    print("###################################")
    print("\nObtain first order coefficients ")
# Get first order coefficients
totquat_adapt = 0
R0 = np.eye(ndim)
pc_model0, c_k0, totquat0 = forward_propagation(x1, x2, x3, x4, x5, x6, x7, \
    std_x1, std_x2, std_x3, std_x4, std_x5, std_x6, std_x7, a, b, \
    nord0, ndim, pc_type, param0, R0, main_verbose, sf="sparse", pc_alpha=0.0, pc_beta=1.0)
# mean value correction using the known 0 order coefficient, which is obtained from first order PCE
coeff0 = c_k0[0]
# Obtain the rotation matrix
R = adp.gauss_adaptation(c_k0[1:ndim+1], ndim, method)
totquat_adapt += totquat0

# perform adaptation starting from 1 dimension, until convergent
converge = False
l2_error_thred = 0.01
s = 1
norm_eta = []
l2_error_eta_succ = []
while not converge:
    if main_verbose > 0:
        print('\n\nAdaptation dimension := ', s)
    # 1 dimension PCE only make sense to using full quadratures
    if s==1:
        pc_model3, c_k3, totquat3 = forward_propagation(x1, x2, x3, x4, x5, x6, x7, \
            std_x1, std_x2, std_x3, std_x4, std_x5, std_x6, std_x7, a, b, \
            nord, s, pc_type,param, R, main_verbose, sf="full", pc_alpha=0.0, pc_beta=1.0)
        # mean value correction
        c_k3[0] = coeff0
        # obtain norm of coefficients in eta space
        norm_eta.append(np.linalg.norm(c_k3))
        totquat_adapt += totquat3

    else:
        pc_model3, c_k3, totquat3 = forward_propagation(x1, x2, x3, x4, x5, x6, x7, \
            std_x1, std_x2, std_x3, std_x4, std_x5, std_x6, std_x7, a, b, \
            nord, s, pc_type,param, R, main_verbose, sf="sparse", pc_alpha=0.0, pc_beta=1.0)
        totquat_adapt += totquat3
        # mean value correction
        c_k3[0] = coeff0
        # obtain l2 error of coefficients of current dimensional (s) expansion and previous dimensional (s-1)
        l2_error_eta_tmp, C1 = adp.l2_error_eta(coeffs, c_k3, s-1, s, nord, pc_type, param, sf, pc_alpha=0.0, pc_beta=1.0)
        # obtain norm of coefficients in eta space
        norm_eta.append(np.linalg.norm(c_k3))
        l2_error_eta_succ.append(l2_error_eta_tmp)
        # if l2 error is less than threashold then exit
        if l2_error_eta_tmp < l2_error_thred:
            converge = True
            if main_verbose > 0:
                print('Number of dimensions in the adaptative PCEs is :', s)

            break

    coeffs = c_k3
    s+=1
    if s > ndim: # if s > ndim then no reduction is achieved
        s = ndim
        if main_verbose > 0:
            print('Full dimension used for adaptation :(')

        break

#print('norm_eta = ', norm_eta)
print('l2_error_eta_succ = ', l2_error_eta_succ)

##################################################################################
##### Projecting Coefficients of Adaptative PCEs from Eta Space to Xi space ######
#####                        Verification use only                          ######
##################################################################################
# using s dimension eta
if main_verbose > 0:
    print("\nStaring verification with full dimension PCE")

s0 = s
if s0<ndim:
    # Obtain coefficients of s dimesnional adaptations
    if s0==1:
        pc_model4, c_k4, totquat4 = forward_propagation(x1, x2, x3, x4, x5, x6, x7, \
            std_x1, std_x2, std_x3, std_x4, std_x5, std_x6, std_x7, a, b, \
            nord, s0, pc_type,param, R, main_verbose=0, sf="full", pc_alpha=0.0, pc_beta=1.0)
    else:
        pc_model4, c_k4, totquat4 = forward_propagation(x1, x2, x3, x4, x5, x6, x7, \
            std_x1, std_x2, std_x3, std_x4, std_x5, std_x6, std_x7, a, b, \
            nord, s0, pc_type,param, R, main_verbose=0, sf="sparse", pc_alpha=0.0, pc_beta=1.0)
    c_k4[0] = coeff0
    coeffs = c_k4
    # Obtain coefficients of ndim dimensional adaptations
    pc_model4, c_k4, totquat4 = forward_propagation(x1, x2, x3, x4, x5, x6, x7, \
    std_x1, std_x2, std_x3, std_x4, std_x5, std_x6, std_x7, a, b, \
    nord, ndim, pc_type,param, R, main_verbose=0, sf="sparse", pc_alpha=0.0, pc_beta=1.0)
    # C1 is the coefficients of s dimensional adaptation projected to full dimension eta space
    l2_error_eta_tmp, C1 = adp.l2_error_eta(coeffs, c_k4, s0, ndim, nord, pc_type, param, sf, pc_alpha=0.0, pc_beta=1.0)
else:
    C1 = c_k3
    C1[0] = coeff0
# transfer C1 from eta (ndim dimension) space to xi space
c_xi = adp.transf_coeffs_xi(C1, nord, ndim, pc_type, param, R, sf, pc_alpha=0.0, pc_beta=1.0)
# plot the coefficients of eta on xi space and compare with full dimesional expansion on xi space
plt.figure(figsize=(10,10))
plt.semilogy(np.abs(c_k2),'b',label='full dimension PCEs', linewidth=2)
plt.semilogy(np.abs(c_xi),'r',label=str(s0)+' dimension adapt PCEs', linewidth=2)
plt.tick_params(axis='both', labelsize=20)
plt.xlabel("PCE terms", size=20)
plt.ylabel("Coefficients of PCE terms", size=20)
plt.legend(loc='upper right', prop={'size':20})
#plt.title('Coefficients in '+r'$\xi$'+' space')
plt.grid(True)
figname = 'Coefficients_in_xi_space_'+str(s0)+'d'+str(ndim)+'d.pdf'
#plt.savefig(figname)
#plt.show()


# Generate germ samples for KDE
germ_samples3=np.random.normal(0,1, (n_MC,s))
# Evaluate the PCE at the germ samples
pce_evals3=d_springs_tools.EvaluatePCE(pc_model3,c_k3,germ_samples3)
#Peform kernel density estimation
xpts_pce3, PDF_data_pce3= d_springs_tools.KDE(pce_evals3)

# Print statements to indicate number of samples used
print("\nMonte Carlo sampling used %s points" %(n_MC))
print("Sparse quadrature method used %s points"%(totquat2))
print("Adaptation method used %s points"%(totquat_adapt))

plt.figure(figsize=(10,10))
plt.plot(xpts_MC, PDF_data_MC, linewidth=2, color='b', label= 'Monte Carlo Sampling')
plt.plot(xpts_pce2, PDF_data_pce2, linewidth=2, color='r', label= 'NISP sparse quadrature method')
plt.plot(xpts_pce3, PDF_data_pce3, linewidth=2, color='g', label= 'NISP '+str(s)+'d linear adaptive method')

# Label Axes
plt.xlabel("Effective modulus", size=16)
plt.ylabel("PDF of effective modulus", size=16)
# Add title
#plt.suptitle("Effective module of d spring series", size=20)
# Change tick size
plt.tick_params(axis='both', labelsize=14)
# Pad tick labels
plt.gca().tick_params(pad=6)
# Create legend
plt.legend(loc='upper right', prop={'size':12})
plt.grid(True)
# Save figure
fig_name="eff_module_d_spring_series_adapt_mt3.pdf"
#plt.savefig(fig_name)
# Show figure
plt.show()
