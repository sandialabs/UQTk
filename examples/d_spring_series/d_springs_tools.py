#!/usr/bin/env python
#=====================================================================================
#
#                      The UQ Toolkit (UQTk) version 3.1.1
#                          Copyright (2021) NTESS
#                        https://www.sandia.gov/UQToolkit/
#                        https://github.com/sandialabs/UQTk
#
#     Copyright 2021 National Technology & Engineering Solutions of Sandia, LLC (NTESS).
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
import math

try:
	import numpy as np
except ImportError:
	print("Numpy module not found")

try:
	from scipy import stats
except ImportError:
	print("Scipy stats module not found")

try:
	import PyUQTk.uqtkarray as uqtkarray
	import PyUQTk.quad as uqtkquad
	import PyUQTk.pce as uqtkpce
	import PyUQTk.tools as uqtktools
except ImportError:
	print("PyUQTk array, quad, PCE, or tools module not found")

#################################################################

def KDE(fcn_evals):
	"""
	Performs kernel density estimation
	Input:
		fcn_evals: numpy array of evaluations of the forward model (values of heat flux Q)
	Output:
		xpts_pce: numpy array of points at which the PDF is estimated.
		PDF_data_pce: numpy array of estimated PDF values.
	"""
	# Perform KDE on fcn_evals
	kern_pce=stats.kde.gaussian_kde(fcn_evals)
	# Generate points at which to evaluate the PDF
	xpts=np.linspace(fcn_evals.min(),fcn_evals.max(),200)
	# Evaluate the estimated PDF at these points
	PDF_data=kern_pce(xpts)
	return xpts, PDF_data

def fwd_model(xx, a, b, main_verbose=0):
    '''
    d-spring series
    '''
    n_samples = np.shape(xx)[0]
    d = np.shape(xx)[1]
    f_samples = []
    for i_s in range(n_samples):
        if main_verbose >0:
            if n_samples < 10:
                print("  Sample",i_s,"out of",n_samples)
            elif i_s % (n_samples/10) == 0:
                print("  Sample",i_s,"out of",n_samples)
        num = 1
        sum_den = 0
        for i in range(d):
            num *= (1 + a * xx[i_s, i] + b * xx[i_s, i] ** 2)
            den = 1
            for j in range(d):
                if i != j:
                    den *= (1 + a * xx[i_s, j] + b * xx[i_s, j] ** 2)
            sum_den += den
        f = 1.0*d/(1+b) * num/sum_den
        f_samples.append(f)
    return np.array(f_samples)

################################################################################
def EvaluatePCE(pc_model,pc_coeffs,germ_samples):
    """
    Evaluate PCE at a set of samples of the germ of this PCE
    Input:
        pc_model: PC object with into about PCE
        pc_coeffs: numpy array with PC coefficients of the RVs to be evaluated.
                   Each column corresponds to one RV.
        germ_samples: numpy array with samples of the PCE grem at which the RVs
                      are to be evaluated. Each line is one sample. The number
                      of colums is the number of RVs.

    Output:
        Numpy array with PCE evaluations
    """

    # Get data set dimensions etc.
    n_test_samples = germ_samples.shape[0]
    ndim = germ_samples.shape[1]
    npce = pc_model.GetNumberPCTerms()

    # Put PC germ samples in a UQTk array
    std_samples_uqtk = uqtkarray.dblArray2D(n_test_samples, ndim)
    std_samples_uqtk.setnpdblArray(np.asfortranarray(germ_samples))

    # Numpy array to store all RVs evaluated from sampled PCEs
    rvs_sampled = np.zeros(n_test_samples)

    # Evaluate PCE for RVs in each dimension
    # Create and fill UQTk array for PC coefficients
    c_k_1d_uqtk = uqtkarray.dblArray1D(npce,0.0)
    for ip in range(npce):
        c_k_1d_uqtk[ip] = pc_coeffs[ip]

    # Create UQTk array to store outputs in
    rv_from_pce_uqtk = uqtkarray.dblArray1D(n_test_samples,0.0)

    # Evaluate the PCEs for reach input RV at those random samples
    pc_model.EvalPCAtCustPoints(rv_from_pce_uqtk,std_samples_uqtk,c_k_1d_uqtk)

    # Put evaluated samples in full 2D numpy array
    for isamp in range(n_test_samples):
        rvs_sampled[isamp] = rv_from_pce_uqtk[isamp]

    # return numpy array of PCE evaluations
    return rvs_sampled
