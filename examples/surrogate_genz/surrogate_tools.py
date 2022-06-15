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
import sys

try:
    import numpy as np
except ImportError:
    print("Need numpy")

try:
    import PyUQTk.uqtkarray as uqtkarray
    import PyUQTk.quad as uqtkquad
    from PyUQTk.utils.func import *
except ImportError:
    print("PyUQTk array and quad module not found")


from math import *
from scipy.stats import qmc 

import PyUQTk.pce as uqtkpce
import PyUQTk.PyPCE.pce_tools as pce_tools

################################################################################
def surrogate(method, nord, ndim, pc_type, pc_alpha, pc_beta, model_genz, nSam, quad_type=None):
    """
    Create a PC surrogate for a Genz function
    
    Input:
        method: Method for surrogate construction; choose "galerkin" or "regression"
        nord: PC order
        ndim: number of dimensions
        pc_type: type of PCE; choose 'LU', 'HG', etc
        pc_alpha:
        pc_beta:
        param: number of quadrature points per dimension or level for sparse quadrature
        model_genz: which genz function to use
        nSam: number of samples to from the PCE
        quad_type: type of quadrature, default is None
        
    Output:
        1D Numpy array with PC coefficients [number PC terms,]
        1D Numpy array with actual outputs of the genz function at sample points [nSam,]
    """

    # Instantiate PC model and random number generator
    pc_model = uqtkpce.PCSet("NISPnoq", nord, ndim, pc_type, pc_alpha, pc_beta)
    rng = qmc.LatinHypercube(d=ndim, seed=42)
    
    # get training points
    if (method=='regression'):
        nTest=int((pc_model.GetNumberPCTerms())*1.1)
        #train_pts=np.random.normal(loc=0, scale=0.5, size=(nTest, ndim))
        train_pts=2*rng.random(n=(int(nTest)))-1
        
    if (method=='galerkin'):
        param=nord+1
        pc_model.SetQuadRule(pc_type, quad_type, param) # set full or sparse
        train_pts, totquat= pce_tools.UQTkGetQuadPoints(pc_model)
    
    # evaluate at training points
    f_evals=func(train_pts,model_genz,np.ones(ndim+1)) #genz-specific
    
    # Obtain Coefficients
    if (method=='regression'):
        c_k = pce_tools.UQTkRegression(pc_model, f_evals, train_pts)
    elif (method=='galerkin'):
        c_k = pce_tools.UQTkGalerkinProjection(pc_model,f_evals)
    
    #germ_samples=np.random.normal(0,1, (nSam,ndim))
    germ_samples=2*rng.random(n=nSam)-1
    pce_evals=pce_tools.UQTkEvaluatePCE(pc_model,c_k,germ_samples)
    f_actual=func(germ_samples,model_genz,np.ones(ndim+1))
   
    return pce_evals, f_actual
################################################################################