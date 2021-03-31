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

import time

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
try:
    import matplotlib.pyplot as plt
    from matplotlib import rc
    rc('mathtext', default='regular')
except ImportError:
    print("Matplotlib not found")

#################################################################
def compute_heat_flux(kw,ka,hi,ho):
    """
    Computes heat flux, and temperature of four window pane surfaces, neglecting conduction
    in air gap between panes, and radiative heat transfer. Uses a solved system of 5 equations.

    Input: Samples (scalars) of the 8 uncertain, Gaussian parameters.
    Output: Heat flux Q, inside window temperature T1, temprature of the interior of first window pane temperature T2,
    temperature of the interior of second window pane T3, and outer window temperature T4.
    """
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

def r_heat_flux(Ts,kw,ka,hi,ho,estimates):
    """
    Function to compute Q,T1,T2,T3 and T4 taking into account convective heat transfer in window air gap as
    well as contributions from radiation.
    Since the convective coefficent for the air gap is unknown, the combined convective and conductive
    heat transfer for this region is calculated by in terms of known parameters for air. A nonlinear system
    of five equation is solved to compute Q.

    Input: Samples of the 8 uncertain, Gaussian parameters
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
    # Solve the nonlinear system of 5 equations using the estimates
    Q,T1,T2,T3,T4 = optimize.fsolve(equations,estimates)
    return(Q)


def GalerkinProjection(pc_model,f_evaluations):
    """
    Obtain PC coefficients by Galerkin Projection
    Input:
        pc_model : PC object with info about basis to project on
        f_evaluations: 1D numpy array (vector) with function to be projected,
                       evaluated at the quadrature points
    Output:
        Numpy array with PC coefficients
    """

    # Get parameters
    if len(f_evaluations.shape) > 1:
        print("This function can only project single variables for now")
        exit(1)
    # Number of quadrature points
    npce = pc_model.GetNumberPCTerms()

    nqp = f_evaluations.shape[0]

    # UQTk array for PC coefficients for one variable
    c_k_1d_uqtk = uqtkarray.dblArray1D(npce,0.0)

    # UQTk array for function evaluations at quadrature points for that variable
    f_uqtk = uqtkarray.dblArray1D(nqp,0.0)
    for ipt in range(nqp):
        f_uqtk[ipt]=f_evaluations[ipt]

    # Galerkin Projection
    pc_model.GalerkProjection(f_uqtk,c_k_1d_uqtk)

    # Put PC coefficients in numpy array
    c_k = np.zeros(npce)
    for ip in range(npce):
        c_k[ip] = c_k_1d_uqtk[ip]

    # Return numpy array of PC coefficients
    return c_k

def evaluate_pce(pc_model,pc_coeffs,germ_samples):
    """
    Evaluate PCE at a set of samples of the germ of this PCE
    Input:
        pc_model: PC object with info about PCE
        pc_coeffs: 1D numpy array with PC coefficients of the RVs to be evaluated.
                   Each column corresponds to one RV.
        germ_samples: numpy array with samples of the PCE germ at which the RVs
                      are to be evaluated. Each line is one sample. The number
                      of columns is the number of RVs.

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

    # Evaluate the PCEs for each input RV at those random samples
    pc_model.EvalPCAtCustPoints(rv_from_pce_uqtk,std_samples_uqtk,c_k_1d_uqtk)

    # Put evaluated samples in numpy array
    for isamp in range(n_test_samples):
        rvs_sampled[isamp] = rv_from_pce_uqtk[isamp]

    # Return numpy array of PCE evaluations
    return rvs_sampled

def get_quadpts(pc_model,ndim):
    """
    Generates quadrature points
    Input:
        pc_model: PC object with info about PCE
        ndim: number of dimensions of the PCE
    Output:
        qdpts: numpy array of quadrature points
    """
    # Get the quadrature points
    qdpts_uqtk = uqtkarray.dblArray2D()
    pc_model.GetQuadPoints(qdpts_uqtk)
    totquat = pc_model.GetNQuadPoints() # Total number of quadrature points

    # Convert quad points to a numpy array
    qdpts = np.zeros((totquat,ndim))
    qdpts_uqtk.getnpdblArray(qdpts)
    return qdpts, totquat

def fwd_model(Ts_samples,kw_samples,ka_samples,hi_samples,ho_samples):
    """
    Evaluates the forward model
    Input:
        Samples of the 8 uncertain, Gaussian inputs
    Output:
        Q_evals: numpy array of evaluations of the forward model

    """
    #Determine number of samples (totquat)
    totquat=len(kw_samples)
    # List to store values of Q (assuming no radiative heat transfer) calculated from
    # the random samples of the parameters
    Q_samples_4PCE=[]
    # List to store values of Q assuming radiative heat transfer occurs
    Q_r_samples_4PCE=[]
    # Calculate values of heat flux Q (assuming no radiative heat transfer)
    # for the different sample values and append to the list
    for i in range(totquat):
        (Q,T1,T2,T3,T4)=compute_heat_flux(kw_samples[i],ka_samples[i], hi_samples[i], ho_samples[i])
        Q_samples_4PCE.append(Q)

        # Calculate values of heat flux Q assuming radiative heat transfer to atmosphere and append to list
        # For the required estimates of Q,T1, and T2 needed to solve the nonlinear system,
        # we use the values obtained by solving the system assuming no radiative heat transfer
        Q2=r_heat_flux(Ts_samples[i],kw_samples[i],ka_samples[i], hi_samples[i], ho_samples[i],(Q,T1,T2,T3,T4))
        Q_r_samples_4PCE.append(Q2)

    # Convert Q_r_samples_4PCE to numpy array
    Q_evals = np.array(Q_r_samples_4PCE)
    return Q_evals


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
    xpts_pce=np.linspace(fcn_evals.min(),fcn_evals.max(),200)
    # Evaluate the estimated PDF at these points
    PDF_data_pce=kern_pce(xpts_pce)
    return xpts_pce, PDF_data_pce

def get_multi_index(pc_model,ndim):
    """
    Function that returns a 2D array of the PC multiindex.
    Input:
       pc_model, ndim.
    Output:
       2D array of the PC multiindex
    """
    # Get number of PC terms
    totpc = pc_model.GetNumberPCTerms()
    # Create  2D int UQTk array with width of ndim and height of totpc
    mi_uqtk = uqtkarray.intArray2D(totpc,ndim)
    # Populate UQTk array with PC multiindex
    pc_model.GetMultiIndex(mi_uqtk)
    # Convert UQTk array to numpy array
    mi = np.zeros((totpc,ndim))
    mi_uqtk.getnpdblArray(mi)
    return mi

def plot_mi_dims(pc_model,c_k,ndim):
    """
    Function that creates a plot of the behavior of the absolute value of the
    PC coefficient for each order.
    Input:
      pc_model, ndim(number of parameters), and c_k(array of pc coefficients).
    Output:
      Matplotlib plot.
    """
    # Get array of PC multiindicies
    mi = get_multi_index(pc_model,ndim)
    # Get the order of the PC coefficient by taking the sum of the multiindex row
    # that corresponds to that value
    misum = np.sum(mi, axis=1)
    # Create an numpy array of the log of the absolute value of the PC coefficients
    cklen = len(c_k)
    ac_k = np.absolute(c_k)
    ac_k = np.log10(ac_k)
    #Create an array to represent the PC coefficient number
    x = np.arange(1,cklen+1)
    #x = np.log10(x)
    #Set the plot size
    plt.figure(figsize=(16,10))
    # Set Plot min, max
    xmin = np.amin(x)
    xmax = np.amax(x)
    ymin = np.amin(ac_k)
    ymax = np.amax(ac_k)
    plt.ylim(ymin,ymax)
    plt.xlim(xmin-2,xmax)
    # Create axis and title labels
    plt.xlabel("Coefficient Number", size=25)
    plt.ylabel("PC Coefficient Magnitude", size=25)
    plt.suptitle("Spectral Decay of the PC Coefficients", size=25)

    y=[-10,-8,-6,-4,-2,0]
    labels = [r'$10^{-10}$',r'$10^{-8}$',r'$10^{-6}$',r'$10^{-4}$',r'$10^{-2}$',r'$0$']
    plt.yticks(y,labels,size=20)
    plt.xticks(size=20)
    # Create verticle lines seperating orders
    dotted_line = plt.Line2D((x[1],x[1]),(ymin,ymax),lw=1,c = 'r')
    dotted_line2 = plt.Line2D((x[7],x[7]),(ymin,ymax),lw=1,c = 'r')
    dotted_line3 = plt.Line2D((x[23],x[23]),(ymin,ymax),lw=1,c = 'r')
    dotted_line4 = plt.Line2D((x[57],x[57]),(ymin,ymax),lw=1,c = 'r')
    plt.gca().add_line(dotted_line)
    plt.gca().add_line(dotted_line2)
    plt.gca().add_line(dotted_line3)
    plt.gca().add_line(dotted_line4)
    # Annotate plot
    plt.annotate(r'$0^{th}$', xy=(x[0],ac_k[100]),xytext=(x[0]-1.7,ac_k[100]),size = 16)
    plt.annotate(r'$1^{st}$', xy=(x[1],ac_k[100]),xytext=(((x[1]+x[7])/2)-1,ac_k[100]),size=20)
    plt.annotate(r'$2^{nd}$', xy=(x[7],ac_k[100]),xytext=((x[7]+x[22])/2,ac_k[100]),size=20)
    plt.annotate(r'$3^{rd}$', xy=(x[23],ac_k[100]),xytext=((x[23]+x[57])/2,ac_k[100]),size=20)
    plt.annotate(r'$4^{th}$', xy=(x[57],ac_k[100]),xytext=((x[57]+x[cklen-1])/2,ac_k[100]),size=20)
    # Plot figure
    plt.plot(x,ac_k,linewidth=2,color='b',)
    # Set plot name, and save as PDF
    fig_name="Multi_Index_Dim.pdf"
    plt.savefig(fig_name)
    print("\nMulti_Index_Dim.pdf has been saved.")
    plt.show()
