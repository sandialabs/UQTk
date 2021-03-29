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
from __future__ import print_function # To make print() in Python 2 behave like in Python 3

# include path for PyUQTk.

import os
src = os.getenv('UQTK_SRC')

import sys
sys.path.append('../pyuqtkarray/')
sys.path.append('../pyuqtkarray_tools/')
sys.path.append('../quad/')
sys.path.append('../pce/')
sys.path.append('../tools')
sys.path.append('../adaptation_tools/')
sys.path.append('../pce_tools/')



try:
    import pyuqtkarray as uqtkarray
    import pyuqtkarray_tools
except ImportError:
    print("PyUQTk array modules not found")

try:
    import _quad as uqtkquad
except ImportError:
    print("PyUQTk quad modules not found")

try:
    import _pce as uqtkpce
    import pce_tools
except ImportError:
    print("PyUQTk pce or pce_tools modules not found")

try:
    import _tools as uqtktools
except ImportError:
    print("PyUQTk tools modules not found")

try:
    import adaptation_tools
except ImportError:
    print("PyUQTk adaptation_tools module not found")

try:
	import numpy as np
except ImportError:
	print("Need numpy to test PyUQTk")

'''
This test use basis adaptation method to do expansion of

f(x) = sum(xi) + 0.25 sum(xi)^2 + 0.025 sum(xi)^3

And 1 dimension adaptation has analytical coefficients

C = [0.25 d, 0.075 d^(3/2)+d^(1/2), 0.25d, 0.025 d^(3/2)] for eta = sum(xi)
or
C = [0.25d, -0.075d^(3/2)-d^(1/2), 0.25d, -0.025d^(3/2)] for eta = -sum(xi)

for 3rd order expansion. Where d is the dimenison of x
'''

####  Numerical result of projected coefficients from eta space to xi space, dimension=5 ####
c_xiRef = np.array([[1.25, 1.375, 1.375, 1.375, 1.375, 1.375, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, \
 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.15, 0.15, 0.15,\
 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15,\
 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15,\
 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15]])

####################################################
def forward_propagation(mu, sigma, \
    nord, ndim, pc_type,param, R, main_verbose, sf="sparse"):
    # obtain pc model
    pc_model = uqtkpce.PCSet("NISP", nord,ndim,pc_type, 0.0, 1.0)
    # set quadrature rule and obtain quadratures
    pc_model.SetQuadRule(pc_type, sf, param)
    npce = pc_model.GetNumberPCTerms() # Number of terms in the PCE
    qdpts, totquat= pce_tools.UQTkGetQuadPoints(pc_model)
    # map the quadrature from eta space to xi space
    qdpts_xi = eta_to_xi_mapping(qdpts, R)
    # map germs to input parameters
    xx = mu + sigma * qdpts_xi
    # evaluate at input parameters
    Q_evals = fwd_model(xx, main_verbose)
    # Using Galerkin projection to compute coefficients
    c_k = pce_tools.UQTkGalerkinProjection(pc_model,Q_evals)
    return pc_model, c_k, totquat

def fwd_model(xx, main_verbose=0):
    '''
    Third order polynomial
    '''
    n_samples = np.shape(xx)[0]
    d = np.shape(xx)[1]
    f_samples = []
    for i_s in range(n_samples):
        if main_verbose>0:
            if n_samples < 10:
                print("  Sample",i_s,"out of",n_samples)
            elif i_s % (n_samples/10) == 0:
                print("  Sample",i_s,"out of",n_samples)

        f = np.sum(xx[i_s,:]) + 0.25 * np.sum(xx[i_s,:])**2 + 0.025 * np.sum(xx[i_s,:])**3
        f_samples.append(f)
    return np.array(f_samples)

########################################################################
####                         initial setup                         #####
########################################################################
main_verbose = 0
nord = 3   # order of PCE
nord0 = 1  # order to obtain Gaussian coefficients
ndim = 5   # dimension
pc_type = "HG"
param = nord+1 # level of quadratures
param0 = 1     # level of quadratures for first order expansion (obtain first order coefficients)
sf = "sparse"  # sparse or full quadrature
mu = 0.0       # mean value of input parameters
sigma = 1.0    # standard deviation of input parameters


########################################################################
####      Check the function that creates rotation matrix          #####
########################################################################
# reference 1 dimension PCE coefficients
CRef = np.array([[0.25*ndim, 0.075*ndim**(3.0/2.0)+np.sqrt(ndim),\
    0.25*ndim, 0.025*ndim**(3.0/2.0)]])
# first order expansion to obtain Gaussian coefficients
R0 = np.eye(ndim)
pc_model0, c_k0, totquat0 = forward_propagation(mu, sigma, \
	nord0, ndim, pc_type, param0, R0,main_verbose, sf)

tol = 1e-6                 # tolerance to check l2 errors
for method in range(4):    # loop over 4 different methods
    # Using different method to obtain rotation matrix and perform 1 dimensional adaptation
    R = gauss_adaptation(c_k0[1:ndim+1], ndim, method)
    pc_model1, c_k1, totquat1 = forward_propagation(mu, sigma, \
        nord, 1, pc_type,param, R, main_verbose, sf="full")

    # compare the first order coefficients with reference ones
    c_ktmp = np.copy(c_k1)
    if R[0,0] < 0:
    	c_ktmp[1] = -c_ktmp[1]
    	c_ktmp[3] = -c_ktmp[3]
    l2 = np.linalg.norm(c_ktmp-CRef)/np.linalg.norm(CRef)
    assert l2<tol, "Basis adaptation fails to obtain correct 1 dimension coefficients :-("

########################################################################
####                    Check other functions                      #####
########################################################################
# Using method =3 to perform full dimesnion adaptation
R = gauss_adaptation(c_k0[1:ndim+1], ndim, method)
pc_model2, c_k2, totquat2 = forward_propagation(mu, sigma, \
    nord, ndim, pc_type,param, R, main_verbose, sf)
# test l2_error_eta by 1 dimensional and full dimension expansions in eta space
# and obtain coefficients C1, which is the projection of 1 dimensional expansion
# coeffients in full dimensional expansion
l2_error_eta, C1 = l2_error_eta(c_k1, c_k2, 1, ndim, nord, pc_type, param, sf, 0.0, 1.0)
# test transf_coeffs_xi by transfer C1 to xi space, which is compared with pre-computed results
c_xi = transf_coeffs_xi(C1, nord, ndim, pc_type, param, R, sf, 0.0, 1.0)
l22 = np.linalg.norm(c_xi- c_xiRef)/np.linalg.norm(c_xiRef)

assert l2_error_eta<tol, "l2_error_eta function fails to obtain correct values :-("
assert l22<tol, "transf_coeffs_xi function fails to obtain correct values :-("
