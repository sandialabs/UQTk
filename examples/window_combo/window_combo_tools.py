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

import sys
sys.path.append('../../PyUQTk/pyuqtkarray/')
sys.path.append('../../PyUQTk/quad/')
sys.path.append('../../PyUQTk/pce/')
sys.path.append('../../PyUQTk/tools/')
sys.path.append('../../PyUQTk')


try:
    import _uqtkarray as uqtkarray
except ImportError:
    print("PyUQTk array module not found")
try:
    import _quad as uqtkquad
except ImportError:
    print("PyUQTk quad module not found")
try:
    import _pce as uqtkpce
except ImportError:
    print("PyUQTk PCE module not found")
try:
    import _tools as uqtktools
except ImportError:
    print("PyUQTk tools module not found")
#################################################################

def compute_heat_flux(samples):
    """
    Computes heat flux, outside wall temp, and inside wall temp
    assuming no radiative heat transfer. Solves a linear system of 3 equations (the forward model).

    Input: An array of uncertain, Gaussian parameters.
    Output: Heat flux Q, inside window temperature T1, and outer window temperature T2

    """
    Ti=samples[0]
    To=samples[1]
    dw=samples[2]
    kw=samples[3]
    hi=samples[4]
    ho=samples[5]

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
################################################################################
def compute_heat_flux2(samples):
    """
    Computes heat flux, and temperature of four window pane surfaces, neglecting conduction
    in air gap between panes, and radiative heat transfer. Uses a solved system of 5 equations.

    Input: An array of samples of the uncertain, Gaussian parameters.
    Output: Heat flux Q, inside window temperature T1, temprature of the interior of first window pane temperature T2,
    temperature of the interior of second window pane T3, and outer window temperature T4.
    """

    hi=samples[1]
    ho=samples[2]
    kw=samples[3]
    ka=samples[4]


    dw = 0.005 # Width of the glass pane (m)
    da = 0.01  # Width of the gap between the panes (m)
    Ti = 293.0 # Room temperature (K)
    To = 273.0 # Outside temperature (K)

    #Consolidation of variables
    A = ho*ka*kw
    B = ho*da*kw
    C = hi*ka*kw
    D = ho*ka*dw

    # Compute T1
    T1 = (A*To+((C+hi*B+2*D*hi)*Ti))/(C+hi*B+2*D*hi+A)
    # Compute T2
    T2 = (A*To+D*hi*To+((C+B*hi+D*hi)*Ti))/(C+hi*B+hi*2*D+A)
    # Compute T3
    T3 = (To*(B*hi+D*hi+A)+(C+hi*D)*Ti)/(C+hi*(B+2*D)+A)
    # Compute T4
    T4 = (((hi*B+2*D*hi+A)*To)+C*Ti)/(C+hi*(B+2*D)+A)
    # Compute Q
    Q = ho*(T4-To)

    return (Q,T1,T2,T3,T4)
###################################################################################
def compute_heat_flux3(samples):
    """
    Computes heat flux, and temperature of four window pane surfaces, neglecting conduction
    in air gap between panes, and radiative heat transfer. Uses a solved system of 5 equations.

    Input: Array of samples (scalars) of the uncertain, Gaussian parameters.
    Output: Heat flux Q, inside window temperature T1, temprature of the interior of first window pane temperature T2,
    temperature of the interior of second window pane T3, and outer window temperature T4.
    """

    Ti=samples[0]
    To=samples[1]
    dw=samples[2]
    da=samples[3]
    kw=samples[4]
    ka=samples[5]
    hi=samples[6]
    ho=samples[7]

    #Consolidation of variables
    A = ho*ka*kw
    B = ho*da*kw
    C = hi*ka*kw
    D = ho*ka*dw

    # Compute T1
    T1 = (A*To+((C+hi*B+2*D*hi)*Ti))/(C+hi*B+2*D*hi+A)
    # Compute T2
    T2 = (A*To+D*hi*To+((C+B*hi+D*hi)*Ti))/(C+hi*B+hi*2*D+A)
    # Compute T3
    T3 = (To*(B*hi+D*hi+A)+(C+hi*D)*Ti)/(C+hi*(B+2*D)+A)
    # Compute T4
    T4 = (((hi*B+2*D*hi+A)*To)+C*Ti)/(C+hi*(B+2*D)+A)
    # Compute Q
    Q = ho*(T4-To)

    return (Q,T1,T2,T3,T4)
###################################################################################
def r_heat_flux(samples, estimates):
    """
    Function to compute Q,T1,and T2 assuming radiative heat transfer occurs.
    Assumes radiative heat transfer to atmosphere and requires solving a nonlinear system of equations

    Input: An array of samples of the 7 uncertain, Gaussian parameters
           estimates: For the required estimates of Q,T1, and T2 needed to solve the nonlinear system,
           we use the values obtained by solving the system assuming no radiative heat transfer

    Output: Heat Flux Q

    """
    Ti=samples[0]
    To=samples[1]
    dw=samples[2]
    kw=samples[3]
    hi=samples[4]
    ho=samples[5]
    TA=samples[6]

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
################################################################################
def r_heat_flux2(samples,estimates):
    """
    Function to compute Q,T1,T2,T3 and T4 taking into account convective heat transfer in window air gap as
    well as contributions from radiation.
    Since the convective coefficent for the air gap is unknown, the combined convective and conductive
    heat transfer for this region is calculated by in terms of known parameters for air. A nonlinear system
    of five equation is solved to compute Q.

    Input: An array of samples of the uncertain, Gaussian parameters
           estimates: For the required estimates of Q,T1,T2,T3 and T4 needed to solve the nonlinear system,
           we use the values obtained by solving the system assuming no convective heat transfer for air gap.

    Output: Heat Flux Q
    """
    def equations(variables):
        Q,T1,T2,T3,T4 = variables
        g = 9.8  # Acceleration due to gravity (m/s^2)
        e = 0.95 # Emissivity of uncoated glass, unitless value
        SBC = 5.67037321e-8 # Stefan-Boltzmann constant (W/m^2*K^4)
        mu = 1.86e-5 # Viscosity of air (kg/m*s)
        rho = 1.29 # Density of air (kg/m^3)
        beta = 3.67e-3 # Coefficient of thermal exspansion (1/K)
        dw = 0.005 # Width of the glass pane (m)
        da = 0.01  # Width of the gap between the panes (m)
        Ti = 293.0 # Room temperature (K)
        To = 273.0 # Outside temperature (K)
        f1 = hi*(Ti-T1)-Q
        f2 = (kw/dw)*(T1-T2)-Q
        f3 = ((0.41*(ka/da)*(((g*beta*rho**2*da**3)/mu**2)*abs(T2-T3))**0.16)*(T2-T3))-Q
        f4 = (kw/dw)*(T3-T4)-Q
        f5 = ho*(T4-To)+(e*SBC*(T4**4-Ts**4))-Q
        return(f1,f2,f3,f4,f5)

    Ts=samples[0]
    hi=samples[1]
    ho=samples[2]
    kw=samples[3]
    ka=samples[4]

    # Solve the nonlinear system of 5 equations using the estimates
    Q,T1,T2,T3,T4 = optimize.fsolve(equations,estimates)
    return(Q)
###################################################################################
def r_heat_flux3(samples,estimates):
    """
    Function to compute Q,T1,T2,T3 and T4 taking into account convective heat transfer in window air gap as
    well as contributions from radiation.
    Since the convective coefficent for the air gap is unknown, the combined convective and conductive
    heat transfer for this region is calculated by in terms of known parameters for air. A nonlinear system
    of five equation is solved to compute Q.

    Input: Array of samples of the 8 uncertain, Gaussian parameters
           estimates: For the required estimates of Q,T1,T2,T3 and T4 needed to solve the nonlinear system,
           we use the values obtained by solving the system assuming no convective heat transfer for air gap.

    Output: Heat Flux Q
    """
    def equations(variables):
        Q,T1,T2,T3,T4 = variables
        g = 9.8  # Acceleration due to gravity (m/s^2)
        e = 0.95 # Emissivity of uncoated glass, unitless value
        SBC = 5.67037321e-8 # Stefan-Boltzmann constant (W/m^2*K^4)
        f1 = hi*(Ti-T1)-Q
        f2 = (kw/dw)*(T1-T2)-Q
        f3 = ((0.41*(ka/da)*(((g*beta*rho**2*da**3)/mu**2)*abs(T2-T3))**0.16)*(T2-T3))-Q
        f4 = (kw/dw)*(T3-T4)-Q
        f5 = ho*(T4-To)+(e*SBC*(T4**4-Ts**4))-Q
        return(f1,f2,f3,f4,f5)

    Ti=samples[0]
    To=samples[1]
    dw=samples[2]
    da=samples[3]
    kw=samples[4]
    ka=samples[5]
    hi=samples[6]
    ho=samples[7]
    Ts=samples[8]
    mu=samples[9]
    rho=samples[10]
    beta=samples[11]

    # Solve the nonlinear system of 5 equations using the estimates
    Q,T1,T2,T3,T4 = optimize.fsolve(equations,estimates)
    return(Q)
###################################################################################
def fwd_model(samples, model, compute_rad, sub_verbose=0):
    """
    Evaluates the forward model
    Input:
        samples: array of uncertain parameters, Gaussian inputs
        sub_verbose: verbosity level [default = 0]
    Output:
        Q_evals: numpy array of evaluations of the forward model

    """
    # Determine number of samples
    n_samples=len(samples[0])

    # List to store values of Q calculated by the random samples of the parameters
    Q_samples=[]

    # Each entry of samples is an array for a single parameter
    # Transposing makes each entry have one value for each parameter
    samples=np.transpose(samples)

    # Calculate values of heat flux Q (assuming no radiative heat transfer)
    # for the different sample values and append to the list
    if sub_verbose > 0:
        print("\nComputing heat transfer for",n_samples,"parameter samples:" )

    # Calculate Q for each set of inputs
    for i_s in range(n_samples):
        if sub_verbose > 0 and i_s % (n_samples/10) == 0:
            print("  Sample",i_s,"out of",n_samples)

        if (model==1): # single-pane window
            (Q,T1,T2)=compute_heat_flux(samples[i_s])
            if (compute_rad): # adds radiation component
                # uses a non-radiation evaluation as a starting estimate
                Q=r_heat_flux(samples[i_s], (Q,T1,T2))

        if (model==2): # double-pane window
            (Q,T1,T2,T3,T4)=compute_heat_flux2(samples[i_s])
            if (compute_rad): #adds radiation component
                #uses a non-radiation evaluation as a starting estimate
                Q=r_heat_flux2(samples[i_s],(Q,T1,T2,T3,T4))

        if (model==3):
            (Q,T1,T2,T3,T4)=compute_heat_flux3(samples[i_s])
            if (compute_rad): #adds radiation component
                #uses a non-radiation evaluation as a starting estimate
                Q=r_heat_flux3(samples[i_s],(Q,T1,T2,T3,T4))

        Q_samples.append(Q)

    #Convert Q_samples to numpy array
    Q_evals = np.array(Q_samples)
    return Q_evals
###################################################################################
