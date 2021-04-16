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
    from numpy import *
except ImportError:
    print("Need numpy")

sys.path.append('../../PyUQTk/pyuqtkarray/')
sys.path.append('../../PyUQTk/pyuqtkarray_tools/')
sys.path.append('../../PyUQTk/quad/')

#try:
import pyuqtkarray as uqtkarray
import pyuqtkarray_tools as uqtkarray_tools
import _quad as uqtkquad
from PyUQTk.utils.func import *
#except ImportError:
    #print("PyUQTk array and quad module not found")


from math import *
import random as rnd
import itertools
import numpy as np
#####################################################

def generate_qw(ndim,param,sp='full',type='LU'):

    import pyuqtkarray as uqtkarray

    # get quad points and weights
    x = uqtkarray.dblArray2D()
    w = uqtkarray.dblArray1D()


    #print 'Create an instance of Quad class'
    q = uqtkquad.Quad(type,sp,ndim,param,0.0,1.0)

    #print 'Now set and get the quadrature rule...'
    q.SetRule()
    q.GetRule(x,w)

    # print out x and w
    #print 'Displaying the quadrature points and weights:\n'
    #print x
    #print w
    n = x.XSize()

    # get quad points
    x_np = zeros((n,ndim))
    x_np = uqtkarray_tools.uqtk2numpy(x)
    #x.getnpdblArray(x_np)

    # get quad weights
    w_np = zeros(n)
    w_np = uqtkarray_tools.uqtk2numpy(w)

    xpts=array((x_np))

    return xpts,w_np


def find_error(pts,ndim,model,integ_ex,func_params):
    #Initialize average error to zero
    avg_error=0
    #Do MC integration 10 times
    for i in range(1,11):
        #Generate random points to evaluate the function at
        x_mc=np.random.uniform(-1,1, size=(pts, ndim))
        #Evaluate function
        mc_ypts=func(x_mc,model,func_params)
        #Initialize mc_int (value of monte carlo integration) to zero
        mc_int=0
        #Average the function values
        for ypt in mc_ypts:
            mc_int+=ypt/float(pts)
        #Find error from exact integral
        mc_error=abs(integ_ex-mc_int)
        #Add error/10 to average error
        avg_error+=mc_error/10.
    return avg_error
