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

# Import Modules
try:
    import uqtkarray
    import quad as uqtkquad
    import pce as uqtkpce
    import tools as uqtktools
except ImportError:
    import PyUQTk.uqtkarray as uqtkarray
    import PyUQTk.quad as uqtkquad
    import PyUQTk.pce as uqtkpce
    import PyUQTk.tools as uqtktools
except ImportError:
    print("PyUQTk array, quad, pce, or tools modules not found")

try:
    import numpy as np
except ImportError:
    print("Numpy module could not be found")

try:
    from scipy import stats
except ImportError:
    print("Scipy stats module could not be found")

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
    pc_model.DrawSampleSet(p, samples) #Does the same thing as EvalPCAtCustPoints?

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
    ndim = samples.shape[1]
    npce = pc_model.GetNumberPCTerms()

    # Put PC samples in a UQTk array
    std_samples_uqtk = uqtkarray.dblArray2D(n_test_samples, ndim)
    std_samples_uqtk = uqtkarray.numpy2uqtk(np.asfortranarray(samples))

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

    Note: need to generalize this to allow projecting multiple variables at the time

    Input:
        pc_model : PC object with info about basis to project on
        f_evaluations: 1D numpy array (vector) with function to be projected,
                       evaluated at the quadrature points [npq,]
    Output:
        1D Numpy array with PC coefficients [npce,]
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
        #f_uqtk[ipt]=f_evaluations[ipt]
        f_uqtk.assign(ipt,f_evaluations[ipt])

    # Galerkin Projection
    pc_model.GalerkProjection(f_uqtk,c_k_1d_uqtk)

    # Put PC coefficients in numpy array
    c_k = np.zeros(npce)
    for ip in range(npce):
        c_k[ip] = c_k_1d_uqtk[ip]
        #c_k[ip] = c_k_1d_uqtk.at(ip)

    # Return numpy array of PC coefficients
    return c_k
################################################################################
def UQTkRegression(pc_model,f_evaluations, samplepts):
    """
    Obtain PC coefficients by regression via UQTk

    Note: need to generalize this to allow projecting multiple variables at the time

    Input:
        pc_model : PC object with info about basis
        f_evaluations: 1D numpy array (vector) with function,
                       evaluated at the quadrature points [nqp,]
    Output:
        1D Numpy array with PC coefficients for each PC term [npce,]
    """

    # Get parameters
    if len(f_evaluations.shape) > 1:
        print("This function can only project single variables for now")
        exit(1)

    npce = pc_model.GetNumberPCTerms()
    nqp = f_evaluations.shape[0]        # Number of sample points
    ndim=samplepts.shape[1]

    #UQTk array for samples - [nqp, ndim]
    sam_uqtk=uqtkarray.dblArray2D(nqp,ndim)
    sam_uqtk=uqtkarray.numpy2uqtk(np.asfortranarray(samplepts))

    # UQTk array for polynomials evaluated at the sample points
    psi_uqtk = uqtkarray.dblArray2D()
    pc_model.EvalBasisAtCustPts(sam_uqtk, psi_uqtk)

    # NumPy array for polynomials evaluated at the sample points
    psi_np = np.zeros( (nqp, npce) )
    psi_np = uqtkarray.uqtk2numpy(psi_uqtk)

    # Regression
    c_k, resids, rank, s = np.linalg.lstsq(psi_np,f_evaluations,rcond=None)

    # Return numpy array of PC coefficients
    return c_k
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
