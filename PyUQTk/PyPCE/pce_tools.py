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

# This is only necessary for the user to pass the ctest without the need for install
import sys
sys.path.append('../pyuqtkarray/')
sys.path.append('../pce/')
sys.path.append('../quad/')
sys.path.append('../tools/')
sys.path.append('../..')


# Import Modules
try:
    import uqtkarray
    import quad as uqtkquad
    import pce as uqtkpce
    import tools as uqtktools
    import bcs as bcs
    import PyUQTk.utils.multiindex as uqtkmi
except ImportError:
    import PyUQTk.uqtkarray as uqtkarray
    import PyUQTk.quad as uqtkquad
    import PyUQTk.pce as uqtkpce
    import PyUQTk.tools as uqtktools
except ImportError:
    print("PyUQTk array, quad, pce, tools, or bcs modules not found")

try:
    import numpy as np
except ImportError:
    print("Numpy module could not be found")

try:
    from scipy import stats
    import math
except ImportError:
    print("Scipy stats or math module could not be found")

try:
    import matplotlib.pyplot as plt
    from matplotlib import rc
    rc('mathtext', default='regular')
except ImportError:
    print("Matplotlib not found")
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
    if (verbose>0):
        print("Total number of quadrature points =",totquat)

    # Set up transpose of input data for the inverse Rosenblatt transformation in a UQTk array
    ydata_t = uqtkarray.dblArray2D(ndim,nsamp)
    #ydata_t.setnpdblArray(np.asfortranarray(rvs_in.T))
    y_data_t=uqtkarray.numpy2uqtk(np.asfortranarray(rvs_in.T))

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
        if (verbose>0):
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
        #p[ip] = pc_coeffs[ip]
        p.assign(ip,pc_coeffs[ip])

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
        1D Numpy array with PCE evaluations [n_test_samples,]
    """

    #need a 1d array passed into pc_coeffs
    if(len(pc_coeffs.shape) != 1):
        print("UQTkEvaluatePCE only takes one PCE. pc_coeff needs to be 1 dimension.")
        exit(1)

    # Get data set dimensions etc.
    n_test_samples = samples.shape[0]
    npce = pc_model.GetNumberPCTerms()

    # Put PC samples in a UQTk array
    if len(samples.shape)>1:
        ndim = samples.shape[1]
        std_samples_uqtk = uqtkarray.dblArray2D(n_test_samples, ndim)
        std_samples_uqtk = uqtkarray.numpy2uqtk(np.asfortranarray(samples))
    else:
        std_samples_uqtk=uqtkarray.dblArray2D(n_test_samples,1) #UQTk array for samples - [nsam, ndim]
        for i in range(n_test_samples):
            std_samples_uqtk.assign(i, 0, samples[i])

    # Create and fill UQTk array for PC coefficients
    c_k_1d_uqtk = uqtkarray.dblArray1D(npce,0.0)
    for ip in range(npce):
        c_k_1d_uqtk.assign(ip,pc_coeffs[ip])

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

    Note: Need to generalize this to allow projecting multiple variables at a time

    Input:
        pc_model : PC object with info about basis to project on
        f_evaluations: 1D numpy array (vector) with function to be projected,
                       evaluated at the quadrature points [npq,]
    Output:
        1D Numpy array with PC coefficients [npce,]
    """

    # Sends error message if y-values are multi-dimensional
    if len(f_evaluations.shape) > 1:
        print("This function can only project single variables for now")
        exit(1)

    # Get parameters
    npce = pc_model.GetNumberPCTerms()  # Number of PC terms
    nqp = f_evaluations.shape[0]        # Number of quadrature points

    # UQTk array for PC coefficients for one variable
    c_k_1d_uqtk = uqtkarray.dblArray1D(npce,0.0)

    # UQTk array for function evaluations at quadrature points for that variable
    f_uqtk = uqtkarray.dblArray1D(nqp,0.0)
    for ipt in range(nqp):
        f_uqtk.assign(ipt,f_evaluations[ipt])

    # Galerkin Projection
    pc_model.GalerkProjection(f_uqtk,c_k_1d_uqtk)

    # Put PC coefficients in numpy array
    c_k = np.zeros(npce)
    for ip in range(npce):
        c_k[ip] = c_k_1d_uqtk[ip]

    # Return numpy array of PC coefficients
    return c_k
################################################################################
def UQTkRegression(pc_model,f_evaluations, samplepts):
    """
    Obtain PC coefficients by regression

    Note: Need to generalize this to allow projecting multiple variables at a time

    Input:
        pc_model :     PC object with info about basis
        f_evaluations: 1D NumPy array (vector) with function evaluated at the
                            sample points [nsam,]
        samplepts:     n-dimensional NumPy array with sample points
                            [nsam, ndim]
    Output:
        1D Numpy array with PC coefficients for each PC term [npce,]
    """

    # Sends error message if y-values are multi-dimensional
    if len(f_evaluations.shape) > 1:
        print("This function can only project single variables for now")
        exit(1)

    # Get parameters
    npce = pc_model.GetNumberPCTerms()  # Number of PC terms
    nsam = f_evaluations.shape[0]       # Number of sample points

    # Create UQTk array for sample points - [nsam, ndim]
    # if dim>1
    if len(samplepts.shape)>1:
        ndim=samplepts.shape[1]         # Number of dimensions
        sam_uqtk=uqtkarray.numpy2uqtk(np.asfortranarray(samplepts))
    # if dim = 1
    else:
        sam_uqtk=uqtkarray.dblArray2D(nsam,1)
        for i in range(nsam):
            sam_uqtk.assign(i, 0, samplepts[i])

    # UQTk array for the basis terms evaluated at the sample points
    psi_uqtk = uqtkarray.dblArray2D()
    pc_model.EvalBasisAtCustPts(sam_uqtk, psi_uqtk)

    # NumPy array for basis terms evaluated at the sample points - [nsam, npce]
    psi_np = uqtkarray.uqtk2numpy(psi_uqtk)

    # Regression
    c_k, resids, rank, s = np.linalg.lstsq(psi_np,f_evaluations,rcond=None)

    # Return numpy array of PC coefficients
    return c_k
################################################################################
def UQTkBCS(pc_model, f_evaluations, samplepts, sigma, eta, nfolds=5, upit=0,\
    conserve=False, verbose=False):
    """
    Obtain PC coefficients by Bayesian compressive sensing

    Note: Need to generalize this to allow multiple variables at a time
    ToDo: add documentation in UQTk manual on what BCS is and the basis growth schemes
    ToDo: generalize to weighted BCS

    Input:
        pc_model :     PC object with information about the starting basis
        f_evaluations: 1D numpy array (vector) with function, evaluated at the
                            sample points [#samples,]
        samplepts:     N-dimensional NumPy array with sample points [#samples,
                            #dimensions]
        sigma:         Inital noise variance we assume is in the data
        eta:           NumPy array, list, or float with the threshold for
                            stopping the algorithm. Smaller values
                            retain more nonzero coefficients. If eta is an array/list,
                            the optimum value of the array is chosen. If a float,
                            the given value is used.
        nfolds:        Number of folds to use for eta cross-validation; default is 5
        upit:          Number of up-iterations; default is 0
        conserve:      Whether to use conservative basis growth; default is a
                            non-conservative approach

    Output:
        c_k:      1D Numpy array with PC coefficients for each term of the final
                       model [#terms_in_final_basis,]
        pc_model: PC object with basis expanded by the up-iterations (if upit>0)

    """

    # Sends error message if y-values are multi-dimensional
    if len(f_evaluations.shape) > 1:
        print("This function can only project single variables for now.")
        exit(1)

    # Choose whether to optimize eta
    if (type(eta)==np.float64 or type(eta)==float):
        eta_opt = eta
    elif (type(eta)==np.ndarray or type(eta)==list):
        eta_opt = UQTkOptimizeEta(pc_model, f_evaluations, samplepts, sigma,\
         eta, conserve, upit, nfolds)
        if verbose:
            print("Optimal eta is", eta_opt)
    else:
        print("Invalid input for eta.")

    # UQTk array for sigma - [1,]
    sig_np=np.array([sigma])
    sig_uqtk=uqtkarray.numpy2uqtk(np.asfortranarray(sig_np))

    #UQTk array for samples - [#samples, #dimensions]
    sam_uqtk=uqtkarray.numpy2uqtk(np.asfortranarray(samplepts))

    # UQTk array for function f_evaluations - [#evaluations,]
    y = uqtkarray.numpy2uqtk(np.asfortranarray(f_evaluations))

    # Initial run of BCS
    weights, used = UQTkEvalBCS(pc_model, y, sam_uqtk, sig_uqtk, eta_opt, verbose)

    # Loop through up-iterations, growing basis with higher order terms
    if(upit>=0):
        for it in range(upit):
            if verbose:
                print("Up-iteration", it+1)
            weights, used, pc_model = UQTkUpItBCS(pc_model, used, sam_uqtk, y,\
                sig_uqtk, eta_opt, conserve, verbose)
    else:
        print("Invalid value for upit.")

    # Create NumPy array of PC coefficients
    c_k=np.zeros(pc_model.GetNumberPCTerms())
    for i in range(used.XSize()):
        pos = used[i]
        c_k[pos]=weights[i]

    # Return numpy array of PC coefficients and new pc model
    return c_k, pc_model
################################################################################
def UQTkOptimizeEta(pc_start, y, x, sigma, etas, conserve, upit, nfolds):
    """
    Choose the opimum eta for Bayesian compressive sensing
    Helper function for UQTkBCS

    Input:
        pc_start :     PC object with information about the starting basis
        y:             1D numpy array (vector) with function, evaluated at the
                            sample points [#samples,]
        x:             N-dimensional NumPy array with sample points [#samples,
                            #dimensions]
        sigma:         Inital noise variance we assume is in the data
        etas:          NumPy array or list with the threshold for stopping the
                            algorithm. Smaller values retain more nonzero
                            coefficients
        conserve:      Whether to use conservative basis growth
        upit:          Number of up-iterations
        nfolds:        Number of folds to use for eta cross-validation

    Output:
        eta_opt:      Optimum eta

    """
    # split data in k folds
    k=bcs.kfoldCV(x, y, nfolds)
    error=np.zeros(nfolds) # minimum error per fold
    e_k=[] # list of optimum etas per fold

    # loop through each fold
    for i in range(nfolds):
        # retrieve training and validation data
        x_tr=k[i]['xtrain']
        y_tr=k[i]['ytrain']
        x_val=k[i]['xval']
        y_val=k[i]['yval']
        RMSE_per_eta=[]

        # loop through each eta
        for eta in etas:

            # Obtain coefficients through BCS
            c_k, pc_final = UQTkBCS(pc_start, y_tr, x_tr, sigma, eta, 1, upit, conserve)

            # Evaluate the PCE at the validation points
            pce_evals = UQTkEvaluatePCE(pc_final, c_k, x_val)

            # Calculate error metric
            MSE = np.square(np.subtract(y_val, pce_evals)).mean()
            RMSE = math.sqrt(MSE)
            RMSE_per_eta.append(RMSE)

        # pick the lowest error and the corresponding eta
        min_error = min(RMSE_per_eta)
        e_k.append(etas[RMSE_per_eta.index(min_error)])

    # return the mean of the optimum etas for each fold
    eta_opt = np.array(e_k).mean()
    return eta_opt
################################################################################
def UQTkEvalBCS(pc_model, y, sam_uqtk, sig_uqtk, eta, verbose):
    """
    Perform one iteration of Bayesian compressive sensing
    Helper function for UQTkBCS

    Input:
        pc_model:  PC object with information about the basis
        y:         1D UQTk array of function evaluations [#samples,]
        sam_uqtk:  N-dimensional UQTk of array of samples [#samples, #dimensions]
        sig_uqtk:  1D UQTk array with the inital noise variance we assume is in
                        the data [1,]
        eta:       Threshold for stopping the algorithm. Smaller values
                            retain more nonzero coefficients.

    Output:
        weights: 1D UQTk array with PC coefficients at indices in used [#used,]
        used:    1D UQTk array with inidices of the sparse weights
    """
    # Configure BCS parameters to defaults
    lambda_init=np.array([]) # Parameter of the Laplace distribution and
                             # coefficient in the l_1 regularization term;
                             # if assigned an empty array, it will be computed,
                             # otherwise lambda will be fixed to the given value.
    adaptive = 0 # Flag for adaptive CS, using a generative basis, set to 0 or 1
    optimal = 1  # Flag for optimal implementation of adaptive CS, set to 0 or 1
    scale = 0.1  # Diagonal loading parameter; relevant only in adaptive,
                    # non-optimal implementation
    bcs_verbose = 0 # silence print statements

    # UQTk array for lambda_init - []
    lam_uqtk=uqtkarray.numpy2uqtk(np.asfortranarray(lambda_init))

    #UQTk array for the basis terms evaluated at the sample points
    psi_uqtk = uqtkarray.dblArray2D()
    pc_model.EvalBasisAtCustPts(sam_uqtk, psi_uqtk)

    # UQTk arrays for outputs
    weights = uqtkarray.dblArray1D()  # sparse weights
    used = uqtkarray.intArray1D()     # position of the sparse weights;
                                          #indices of selected basis terms
    errbars = uqtkarray.dblArray1D()  # 1 standard dev around sparse weights
    basis = uqtkarray.dblArray1D()    # if adaptive==1, basis = next projection
                                          #vector
    alpha = uqtkarray.dblArray1D()    # sparse hyperparameters (1/gamma)
    _lambda = uqtkarray.dblArray1D(1,0.0) # parameter controlling the sparsity
                                          # on output

    # Run BCS through the c++ implementation
    bcs.BCS(psi_uqtk, y, sig_uqtk, eta, lam_uqtk, adaptive, optimal, scale,\
     bcs_verbose, weights, used, errbars, basis, alpha, _lambda)

    # Print result of the BCS iteration
    if (verbose):
        print("BCS has selected", used.XSize(), "basis terms out of",\
            pc_model.GetNumberPCTerms())

    # Return coefficients and their locations with respect to the basis terms
    return weights, used
################################################################################
def UQTkUpItBCS(pc_model, used, sam_uqtk, y, sig_uqtk, eta, conserve, verbose):
    """
    Perform one up-iteration of Bayesian compressive sensing, growing the basis
        with higher order terms and determining coefficients
    Helper function for UQTkBCS

    Input:
        pc_model:  PC object with info about the basis
        used:      UQTk array of the indices of basis terms selected by BCS
        sam_uqtk:  N-dimensional UQTk of array of samples [#samples, #dimensions]
        y:         1D UQTk array of function evaluations [#samples,]
        sig_uqtk:  1D UQTk array with the inital noise variance we assume is in
                        the data [1,]
        eta:       Threshold for stopping the algorithm. Smaller values
                            retain more nonzero coefficients.
        conserve:  Whether to use a conservative approach to basis growth
        verbose:   Boolean flag for print statements

    Output:
        weights: 1D UQTk array with PC coefficients at indices in used [#used,]
        used:    1D UQTk array with inidices of the sparse weights
        pc_model: PC object with a basis expanded by the up-iteration

    """

    # UQTk array for multiindex - [#terms, #dimensions]
    mi_uqtk = uqtkarray.intArray2D(pc_model.GetNumberPCTerms(), sam_uqtk.YSize())

    # Retrieve multiindex
    pc_model.GetMultiIndex(mi_uqtk)

    # UQTk array new multiindex
    used_mi_uqtk=uqtkarray.intArray2D()

    # Calculate new multiindex from the rows that are used
    uqtkarray.subMatrix_row_int(mi_uqtk, used, used_mi_uqtk)

    # Convert to NumPy array
    used_mi_np=uqtkarray.uqtk2numpy(used_mi_uqtk)

    # Grow the multiindex in pertinent directions
    # using a conservative approach
    if (conserve):
        mi_next = np.array(uqtkmi.mi_addfront_cons(used_mi_np), dtype=object)[0]
    # using a non-conservative approach
    else:
        mi_next = np.array(uqtkmi.mi_addfront(used_mi_np),dtype=object)[0]

    if verbose:
        print(mi_next.shape[0]-used.XSize(), "terms added; new multiindex has", mi_next.shape[0], "terms")

    # Convert to UQTk array
    mi_next_uqtk = uqtkarray.intArray2D(mi_next.shape[0], mi_next.shape[1])
    for i in range(mi_next.shape[0]):
        for j in range(mi_next.shape[1]):
            mi_next_uqtk.assign(i, j, mi_next[i][j])

    # Create a PC object with the new multiindex
    pc_model=uqtkpce.PCSet("NISPnoq", mi_next_uqtk, pc_model.GetPCType(),\
            pc_model.GetAlpha(), pc_model.GetBeta())

    # Run an iteraion of BCS using the new basis
    weights, used = UQTkEvalBCS(pc_model, y, sam_uqtk, sig_uqtk, eta, verbose)

    # Return coefficients, their locations with respect to the the basis terms,
        # and a PC object with the updated basis
    return weights, used, pc_model
################################################################################
def UQTkGetQuadPoints(pc_model):
    """
    Generates quadrature points through UQTk and returns them in numpy array
    Input:
        pc_model: PC object with info about PCE
    Output:
        qdpts: numpy array of quadrature points [totquat, n_dim]
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
    #qdpts_uqtk.getnpdblArray(qdpts)
    qdpts = uqtkarray.uqtk2numpy(qdpts_uqtk)
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
def UQTkGSA(pc_model, pc_coeffs):
    """
    Computes Sobol' sensivity indices

    ToDo: refer to documentation in the UQTk manual

    Input:
        pc_model: PC object with information about the basis
        pc_coeffs: NumPy array of PC coefficients [#PCTerms,]
    Output:
        mainsens:  1D NumPy array of the main sensitivities for each dimension [#dim,]
        totsens:   1D NumPy array of the total sensivities for each dimension [#dim,]
        jointsens: 2D NumPy array of joint sensitivities for each pair of dimensions [#dim, #dim]
    """
    # coefficients in a uqtk array
    coef_uqtk = uqtkarray.numpy2uqtk(pc_coeffs)

    # Compute main sensitivities
    mainsens_uqtk=uqtkarray.dblArray1D()
    pc_model.ComputeMainSens(coef_uqtk,mainsens_uqtk)
    mainsens = uqtkarray.uqtk2numpy(mainsens_uqtk)

    # Compute total sensitivities
    totsens_uqtk=uqtkarray.dblArray1D()
    pc_model.ComputeTotSens(coef_uqtk,totsens_uqtk)
    totsens = uqtkarray.uqtk2numpy(totsens_uqtk)

    # Compute joint sensitivities
    jointsens_uqtk = uqtkarray.dblArray2D()
    pc_model.ComputeJointSens(coef_uqtk,jointsens_uqtk)
    for id in range(pc_model.GetNDim()):
        jointsens_uqtk.assign(id, id, mainsens_uqtk[id])
    jointsens = uqtkarray.uqtk2numpy(jointsens_uqtk)

    return mainsens, totsens, jointsens
################################################################################
def UQTkKDE(fcn_evals):
    """
    Performs kernel density estimation
    Input:
        fcn_evals: numpy array of evaluations of the forward model [n_samples,]
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
################################################################################
def UQTkGetMultiIndex(pc_model,ndim):
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
    mi = uqtkarray.uqtk2numpy(mi_uqtk)
    #mi_uqtk.getnpdblArray(mi)
    return mi
################################################################################
def UQTkPlotMiDims(pc_model,c_k,ndim, nord, type):
    """
    Function that creates a plot of the behavior of the absolute value of the
    PC coefficient for each order.
    Input:
        pc_model, ndim(number of parameters), and c_k(array of pc coefficients).
        Order of the PC
        string indicating the type of plot for labeling purposes
    Output:
        Matplotlib plot.
    """
    # Get array of PC multiindicies
    mi = UQTkGetMultiIndex(pc_model,ndim)

    # Get the order of the PC coefficient by taking the sum of the multiindex row
    # that corresponds to that value
    misum = np.sum(mi, axis=1)

    #find values that separate the orders
    sep=[]
    for i in range(misum.shape[0]):
        if misum[i]!=misum[i-1]:
            sep.append(i)

    # Create an numpy array of the log of the absolute value of the PC coefficients
    cklen = len(c_k)
    ac_k = np.absolute(c_k)
    ac_k = np.log10(ac_k)

    # Create an array to represent the PC coefficient number
    x = np.arange(1,cklen+1)

    # Set the plot size
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
    sup="Spectral Decay of the PC Coefficients for "+type+" Quadrature"
    plt.suptitle(sup, size=25)

    # Get the correct number of y-labels
    y=[]
    val=0
    while (val>ymin-2):
        y.insert(0, val)
        val=val-2
    labels=[]
    for val in y:
        new=r'$10^{'+str(val)+'}$'
        labels.append(new)

    plt.yticks(y,labels,size=20)
    plt.xticks(size=20)


    # Create verticle lines seperating orders
    for i in range(nord+1):
        if (i==1):
            label=r'$'+str(i)+'^{st}$'
        elif (i==2):
            label=r'$'+str(i)+'^{nd}$'
        elif (i==3):
            label=r'$'+str(i)+'^{rd}$'
        else:
            label=r'$'+str(i)+'^{th}$'

        if i>0:
            dotted_line=plt.Line2D((sep[i],sep[i]), (y[0],ymax), lw=1, c='r')
            plt.gca().add_line(dotted_line)

        if (sep[i]==sep[-1]):
            plt.annotate(label, xy=(sep[i],ac_k[2]),xytext=((sep[i]+cklen)/2,ac_k[2]),size = 16)
        else:
            plt.annotate(label, xy=(sep[i],ac_k[2]),xytext=((sep[i+1]+sep[i])/2,ac_k[2]),size = 16)

    # Plot figure
    plt.plot(x,ac_k,linewidth=2,color='b',)
    # Set plot name, and save as PDF
    #fig_name="Multi_Index_Dim_"+type+".pdf"
    #plt.savefig(fig_name)
    #print("\n"+fig_name+" has been saved.")
    plt.show()
