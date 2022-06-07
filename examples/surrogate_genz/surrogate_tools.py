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

import math  

import PyUQTk.pce as uqtkpce
import PyUQTk.PyPCE.pce_tools as pce_tools
from PyUQTk.utils.func import *

################################################################################
def surrogate(nord, ndim, pc_type, pc_alpha, pc_beta, param, quad_type, model_genz, coef_type, nSam):
    """
    Create a PC surrogate for a Genz function
    
    Input:
        nord:
        ndim:
        pc_type:
        pc_alpha:
        pc_beta:
        param:
        quad_type:
        model_genz:
        coef_type:
        nSam:
    Output:
        1D Numpy array with PC coefficients
    """

    # Instantiate PC model
    pc_model = uqtkpce.PCSet("NISPnoq", nord, ndim, pc_type, pc_alpha, pc_beta)
    pc_model.SetQuadRule(pc_type, quad_type, param) # set full or sparse
    
    # Get and evaluate quadrature points
    qdpts, totquat= pce_tools.UQTkGetQuadPoints(pc_model)
    f_evals=func(qdpts,model_genz,np.ones(ndim+1)) #genz-specific
    
    # Obtain Coefficients
    if (coef_type=='regression'):
        c_k = pce_tools.UQTkRegression(pc_model, f_evals)
    elif (coef_type=='galerkin'):
        c_k = pce_tools.UQTkGalerkinProjection(pc_model,f_evals)
    else:
        print("Invalid input. Galerkin projection used by default.")
        c_k = pce_tools.UQTkGalerkinProjection(pc_model,f_evals)
    
    germ_samples=np.random.normal(0,1, (nSam,ndim))
    pce_evals=pce_tools.UQTkEvaluatePCE(pc_model,c_k,germ_samples)
    f_actual=func(germ_samples,model_genz,np.ones(ndim+1))
   
    return pce_evals, f_actual
################################################################################