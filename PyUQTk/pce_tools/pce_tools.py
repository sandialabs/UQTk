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
import sys
sys.path.append('../pyuqtkarray/')
sys.path.append('../pce/')
sys.path.append('../quad/')
sys.path.append('../tools/')
sys.path.append('../pce_tools/')

# Import Modules
try:
    import pyuqtkarray as uqtkarray
    import _quad as uqtkquad
    import _pce as uqtkpce
    import tools as uqtktools
except ImportError:
    import PyUQTk.pyuqtkarray as uqtkarray
    import PyUQTk._quad as uqtkquad
    import PyUQTk._pce as uqtkpce
    import PyUQTk.tools as uqtktools
except ImportError:
    print("PyUQTk array, quad, pce, or tools modules not found")

try:
    import numpy as np
except ImportError:
    print("Numpy module could not be found")

################################################################################
def UQTkMap2PCE(pc_model,rvs_in,verbose=0):
    """Obtain PC representation for the random variables that are described by samples.
    Employ a Rosenblatt transformation to build a map between the input RVs and the space
    of the PC germ.
    Input:
        pc_model: object with properties of the PCE to be constructed
        rvs_in    : numpy array with input RV samples. Each line is a sample. The columns
                  represent the dimensions of the input RV.
        verbose : verbosity level (more output for higher values)

    Output:
        Numpy array with PC coefficients for each RV in the original rvs_in input
    """

    # Dimensionality and number of samples of input RVs
    ndim = rvs_in.shape[1]
    nsamp = rvs_in.shape[0]

    # Algorithm parameters
    bw = -1 # KDE bandwidth for Rosenblatt (on interval 0.1)
    iiout = 50 # interval for output to screen

    # Number of PCE terms
    npce = pc_model.GetNumberPCTerms()

    # Get the default quadrature points
    qdpts = uqtkarray.dblArray2D()
    pc_model.GetQuadPoints(qdpts)

    totquat = pc_model.GetNQuadPoints()
    print("Total number of quadrature points =",totquat)

    # Set up transpose of input data for the inverse Rosenblatt transformation in a UQTk array
    ydata_t = uqtkarray.dblArray2D(ndim,nsamp)
    ydata_t.setnpdblArray(np.asfortranarray(rvs_in.T))

    # Set up numpy array for mapped quadrature points
    invRosData = np.zeros((totquat,ndim))

    # Map all quadrature points in chosen PC set to the distribution given by the data
    # using the inverse Rosenblatt transformation
    for ipt in range(totquat):

        # print("Converting quadrature point #",ipt)

        # Set up working arrays
        quadunif = uqtkarray.dblArray1D(ndim,0.0)
        invRosData_1s = uqtkarray.dblArray1D(ndim,0.0)

        # First map each point to uniform[0,1]
        # PCtoPC maps to [-1,1], which then gets remapped to [0,1]
        for idim in range(ndim):
            quadunif[idim] = (uqtktools.PCtoPC(qdpts[ipt,idim],pc_model.GetPCType(),pc_model.GetAlpha(),pc_model.GetBeta(),"LU",0.0,0.0)+1.0)/2.0

        # Map each point from uniform[0,1] to the distribution given by the original samples via inverse Rosenblatt
        if bw > 0:
            uqtktools.invRos(quadunif,ydata_t,invRosData_1s,bw)
        else:
            uqtktools.invRos(quadunif,ydata_t,invRosData_1s)

        # Store results
        for idim in range(ndim):
            invRosData[ipt,idim] = invRosData_1s[idim]

        # Screen diagnostic output
        if ((ipt+1)%iiout == 0) or ipt==0 or (ipt+1)==totquat:
            print("Inverse Rosenblatt for Galerkin projection:",(ipt+1),"/",totquat,"=",(ipt+1)*100/totquat,"% completed")

    # Get PC coefficients by Galerkin projection
    # Set up numpy array for PC coefficients (one column for each transformed random variable)
    c_k = np.zeros((npce,ndim))

    # Project each random variable one by one
    # Could replace some of this with the UQTkGalerkinProjection function below
    for idim in range(ndim):

        # UQTk array for PC coefficients for one variable
        c_k_1d = uqtkarray.dblArray1D(npce,0.0)

        # UQTk array for map evaluations at quadrature points for that variable
        invRosData_1d = uqtkarray.dblArray1D(totquat,0.0)
        # invRosData_1d.setnpdblArray(np.asfortranarray(invRosData[:,idim])
        for ipt in range(totquat):
            invRosData_1d[ipt]=invRosData[ipt,idim]

        # Galerkin Projection
        pc_model.GalerkProjection(invRosData_1d,c_k_1d)

        # Put coefficients in full array
        for ip in range(npce):
            c_k[ip,idim] = c_k_1d[ip]

    # Return numpy array of PC coefficients
    return c_k
################################################################################
def UQTkEvalPC(pce_model,pce_coeffs,germ_sample):
    """
    Use UQTkEvaluatePCE instead
    """
    print("Use UQTkEvaluatePCE instead of UQTkEvalPC.")
    exit(1)
################################################################################
def UQTkDrawSamplesPCE(pc_model,pc_coeffs,n_samples):
    """
    Draw samples of the germ underneath the pc_model and evaluates one PCE
    for those samples.

    Input:
        pc_model:   PC object with into about PCE
        pc_coeffs:  1D numpy array with PC coefficients of the RV to be evaluated. [n_dim]
        n_samples:    number of samples to be drawn
    Output:
        1D Numpy array with PCE evaluations

    """

    #need a 1d array passed into pc_coeffs
    if(len(pc_coeffs.shape) != 1):
        print("UQTkEvaluatePCE only takes one PCE. pc_coeff needs to be 1 dimension.")
        exit(1)

    #get number of nTerms
    npce = pc_model.GetNumberPCTerms()

    # Create and fill UQTk array for PC coefficients
    p = uqtkarray.dblArray1D(npce,0.0)
    for ip in range(npce):
        p[ip] = pc_coeffs[ip]

    #create UQTk array to store outputs in
    samples = uqtkarray.dblArray1D(n_samples,0.0)

    #draw the samples
    pc_model.DrawSampleSet(p, samples)

    #convert samples to a numpy array
    pce_samples = np.zeros(n_samples)
    for isamp in range(n_samples):
        pce_samples[isamp] = samples[isamp]

    #return samples in numpy array
    return pce_samples
################################################################################
def UQTkEvaluatePCE(pc_model,pc_coeffs,samples):
    """
    Evaluate PCE at a set of samples of this PCE
    Input:
        pc_model:   PC object with into about PCE
        pc_coeffs:  1D numpy array with PC coefficients of the RV to be evaluated. [npce]
        samples:    2D numpy array with samples of the PCE at which the RV
                    are to be evaluated. Each line is one sample. [n_samples, ndim]
    Output:
        1D Numpy array with PCE evaluations
    """

    #need a 1d array passed into pc_coeffs
    if(len(pc_coeffs.shape) != 1):
        print("UQTkEvaluatePCE only takes one PCE. pc_coeff needs to be 1 dimension.")
        exit(1)

    # Get data set dimensions etc.
    n_test_samples = samples.shape[0]
    ndim = samples.shape[1]
    npce = pc_model.GetNumberPCTerms()

    # Put PC samples in a UQTk array
    std_samples_uqtk = uqtkarray.dblArray2D(n_test_samples, ndim)
    std_samples_uqtk.setnpdblArray(np.asfortranarray(samples))

    # Create and fill UQTk array for PC coefficients
    c_k_1d_uqtk = uqtkarray.dblArray1D(npce,0.0)
    for ip in range(npce):
        c_k_1d_uqtk[ip] = pc_coeffs[ip]

    # Create UQTk array to store outputs in
    rv_from_pce_uqtk = uqtkarray.dblArray1D(n_test_samples,0.0)

    # Evaluate the PCEs for reach input RV at those random samples
    pc_model.EvalPCAtCustPoints(rv_from_pce_uqtk,std_samples_uqtk,c_k_1d_uqtk)

    # Numpy array to store all RVs evaluated from sampled PCEs
    rvs_sampled = np.zeros((n_test_samples,))

    # Put evaluated samples in full 2D numpy array
    for isamp in range(n_test_samples):
        rvs_sampled[isamp] = rv_from_pce_uqtk[isamp]

    # return numpy array of PCE evaluations
    return rvs_sampled
################################################################################
def UQTkGalerkinProjection(pc_model,f_evaluations):
    """
    Obtain PC coefficients by Galerkin Projection via UQTk

    Note: need to generalize this to allow projecting multiple variables at the time

    Input:
        pc_model : PC object with info about basis to project on
        f_evaluations: 1D numpy array (vector) with function to be projected,
                       evaluated at the quadrature points
    Output:
        1D Numpy array with PC coefficients
    """

    # Get parameters
    if len(f_evaluations.shape) > 1:
        print("This function can only project single variables for now")
        exit(1)

    npce = pc_model.GetNumberPCTerms()
    nqp = f_evaluations.shape[0]        # Number of quadrature points

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
################################################################################
def UQTkGetQuadPoints(pc_model):
    """
    Generates quadrature points through UQTk and returns them in numpy array
    Input:
        pc_model: PC object with info about PCE
    Output:
        qdpts: numpy array of quadrature points [totquat,n_dim]
        totquat: total number of quadrature points
    """

    # Info about PCE
    n_dim = pc_model.GetNDim()

    # Get the quadrature points
    qdpts_uqtk = uqtkarray.dblArray2D()
    pc_model.GetQuadPoints(qdpts_uqtk)
    totquat = pc_model.GetNQuadPoints() # Total number of quadrature points

    # Convert quad points to a numpy array
    qdpts = np.zeros((totquat,n_dim))
    qdpts_uqtk.getnpdblArray(qdpts)
    return qdpts, totquat
################################################################################
def UQTkStDv(pc_model,pc_coeffs):
    """
    Compute Standard Deviation of a PCE through UQTk
    Input:
        pc_model: PC object with info about PCE
        pc_coeffs: 1D numpy array with PCE coefficients
    Output:
        pc_stdv: Standard Deviation of the PCE
    """

    # Info about PCE
    n_pce = pc_model.GetNumberPCTerms()

    # Create and fill UQTk array for PC coefficients
    c_k_1d_uqtk = uqtkarray.dblArray1D(n_pce,0.0)
    for ip in range(n_pce):
        c_k_1d_uqtk[ip] = pc_coeffs[ip]

    pc_stdv = pc_model.StDv(c_k_1d_uqtk)

    return pc_stdv
################################################################################
