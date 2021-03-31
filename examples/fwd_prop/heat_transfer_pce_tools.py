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
from __future__ import print_function #so the print statements in python2 looks right

import math

try:
    import numpy as np
except ImportError:
    print("Numpy module could not be found")

try:
    from scipy import stats
except ImportError:
    print("Scipy stats module could not be found")

try:
    from scipy import optimize
except ImportError:
    print("Scipy optimize module could not be found")

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
#################################################################

def compute_heat_flux(Ti,To,dw,kw,hi,ho):
    """
    Computes heat flux, outside wall temp, and inside wall temp
    assuming no radiative heat transfer. Solves a linear system of 3 equations (the forward model).

    Input: Samples (scalars) of the 6 uncertain, Gaussian parameters.
    Output: Heat flux Q, inside window temperature T1, and outer window temperature T2

    """
    # Commonly used dimensionless terms
    dhko = dw*ho/kw
    dhki = dw*hi/kw

    # Compute T2, the outside window temp
    T2 = (dhki*Ti/(1.0+dhki) + dhko*To)/(dhko + dhki/(1+dhki))

    # Compute T1, the inside window temp
    T1 = (dhki*Ti + T2)/(dhki + 1.0)

    # Compute Q, the conductive heat flux
    Q=kw*(T1-T2)/dw


    return (Q,T1,T2)

###################################################################################
def r_heat_flux(Ti,To,dw,kw,hi,ho,TA, estimates):
    """
    Function to compute Q,T1,and T2 assuming radiative heat transfer occurs.
    Assumes radiative heat transfer to atmosphere and requires solving a nonlinear system of equations

    Input: Samples of the 7 uncertain, Gaussian parameters
           estimates: For the required estimates of Q,T1, and T2 needed to solve the nonlinear system,
           we use the values obtained by solving the system assuming no radiative heat transfer

    Output: Heat Flux Q

    """
    def equations(variable):
        Q,T1,T2 = variable
        # Returns the three forward model equations given values of Q,T1, and T2
        # Used by optimize.fsolve function to solve the nonlinear system
        e=0.95 # Emissivity of uncoated glass
        SBC=5.67037321e-8 # Stefan-Boltzmann constant
        f1=hi*(Ti-T1)-Q
        f2=(kw/dw)*(T1-T2)-Q
        f3=ho*(T2-To)+e*SBC*(T2**4-TA**4)-Q
        return (f1,f2,f3)
    # Solve the nonlinear system of 3 equations using the estimates
    Q,T1,T2= optimize.fsolve(equations,estimates)
    return(Q)

###################################################################################
def fwd_model(Ti_samples,To_samples, dw_samples, kw_samples,hi_samples,ho_samples,sub_verbose=0):
    """
    Evaluates the forward model
    Input:
        Samples of the 6 uncertain, Gaussian inputs
        sub_verbose: verbosity level [default = 0]
    Output:
        Q_evals: numpy array of evaluations of the forward model

    """
    # Determine number of samples
    n_samples=len(Ti_samples)

    # List to store values of Q calculated from
    # the random samples of the parameters
    Q_samples=[]

    # Calculate values of heat flux Q (assuming no radiative heat transfer)
    # for the different sample values and append to the list
    if sub_verbose > 0:
        print("\nComputing heat transfer for",n_samples,"parameter samples:")

    for i_s in range(n_samples):
        if sub_verbose > 0 and i_s % (n_samples/10) == 0:
            print("  Sample",i_s,"out of",n_samples)

        (Q,T1,T2)=compute_heat_flux(Ti_samples[i_s], To_samples[i_s], dw_samples[i_s],\
            kw_samples[i_s], hi_samples[i_s], ho_samples[i_s])
        Q_samples.append(Q)

    #Convert Q_samples to numpy array
    Q_evals = np.array(Q_samples)
    return Q_evals

###################################################################################
def fwd_model_rad(Ti_samples,To_samples, dw_samples, kw_samples,hi_samples,ho_samples,TA_samples,sub_verbose=0):
    """
    Evaluates the forward model with radiation
    Input:
        Samples of the 7 uncertain, Gaussian inputs
        sub_verbose: verbosity level [default = 0]
    Output:
        Q_evals: numpy array of evaluations of the forward model

    """
    # Determine number of samples
    n_samples=len(Ti_samples)

    # List to store values of Q (assuming no radiative heat transfer) calculated from
    # the random samples of the parameters
    Q_samples=[]

    # Calculate values of heat flux Q (assuming no radiative heat transfer)
    # for the different sample values and append to the list
    if sub_verbose > 0:
        print("\nComputing heat transfer with radiation for",n_samples,"parameter samples:")

    for i_s in range(n_samples):
        if sub_verbose > 0 and i_s % (n_samples/10) == 0:
            print("  Sample",i_s,"out of",n_samples)

        # First compute without radiation (to be used as initial guesses)
        (Q,T1,T2)=compute_heat_flux(Ti_samples[i_s], To_samples[i_s], dw_samples[i_s],\
            kw_samples[i_s], hi_samples[i_s], ho_samples[i_s])

        # Calculate values of heat flux Q assuming radiative heat transfer to atmosphere and append to list
        # For the required estimates of Q,T1, and T2 needed to solve the nonlinear system,
        # we use the values obtained by solving the system assuming no radiative heat transfer
        Q2=r_heat_flux(Ti_samples[i_s], To_samples[i_s], dw_samples[i_s], kw_samples[i_s],\
            hi_samples[i_s], ho_samples[i_s], TA_samples[i_s], (Q,T1,T2))
        Q_samples.append(Q2)

    #Convert Q_samples to numpy array
    Q_evals = np.array(Q_samples)
    return Q_evals

###################################################################################
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
