#!/usr/bin/env python
#=====================================================================================
#                     The UQ Toolkit (UQTk) version 3.0.4
#                     Copyright (2017) Sandia Corporation
#                     http://www.sandia.gov/UQToolkit/
#
#     Copyright (2017) Sandia Corporation. Under the terms of Contract DE-AC04-94AL85000
#     with Sandia Corporation, the U.S. Government retains certain rights in this software.
#
#     This file is part of The UQ Toolkit (UQTk)
#
#     UQTk is free software: you can redistribute it and/or modify
#     it under the terms of the GNU Lesser General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
#
#     UQTk is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU Lesser General Public License for more details.
#
#     You should have received a copy of the GNU Lesser General Public License
#     along with UQTk.  If not, see <http://www.gnu.org/licenses/>.
#
#     Questions? Contact Bert Debusschere <bjdebus@sandia.gov>
#     Sandia National Laboratories, Livermore, CA, USA
#===================================================================================== 

import os
import getopt
import shutil
import sys
import numpy as np
import math

sys.path.append(os.environ['UQTK_SRC'])
from PyUQTk.utils.func import *


#############################################################
#############################################################

def model(modelPar):
    #return model_example(modelPar)
    return model_genz(modelPar)
#############################################################

def model_genz(modelPar):
    npar=modelPar.shape[0]
    mdim=modelPar.shape[1]
    genzparams=1./np.arange(1,mdim+2)
    #genzparams[1:]=1./(genzparams[1:]+1)
    print genzparams
    modelnames=['genz_gaus','genz_exp']


    nout=len(modelnames)

    output=np.empty((npar,nout))
    for j in range(nout):
        output[:,j]=func(modelPar,modelnames[j],genzparams)

    return output

# #############################################################
     
def model_example(modelPar):

    npar=modelPar.shape[0]
    mdim=modelPar.shape[1]
    
    nout=7
    
    # Create the design parameters #TODO maybe should be done outside
    designPar=np.array(range(nout)).reshape(-1,1) #/float(nout-1)
    np.savetxt('designPar.dat',designPar)
    
    output=np.empty((npar,nout))
    for i in range(npar):
        print i+1, ": Running the model with parameter setting ", modelPar[i,:]
        for j in range(nout):
            aa=np.dot(modelPar[i,:]+modelPar[i,:]**2,pow(np.arange(1,mdim+1),-designPar[j]-1.))
            output[i,j]=aa *(sum(modelPar[i,:]))

    return output

#############################################################

def main(argv):
    modelPar_file=argv[0]
    output_file=argv[1]
    #designPar_file=argv[2]
    modelPar=np.loadtxt(modelPar_file,ndmin=2) 
    
    output=model(modelPar)
    np.savetxt(output_file,output)

if __name__ == "__main__":
    main(sys.argv[1:])






