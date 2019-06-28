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
try:
    from numpy import *
except ImportError:
    "Need numpy"

try:
    import PyUQTk.uqtkarray as uqtkarray
    import PyUQTk.quad as uqtkquad
    from PyUQTk.utils.func import *
except ImportError:

    print "PyUQTk array and quad module not found"


from math import *
import random as rnd
import itertools
import numpy as np
#####################################################

def generate_qw(ndim,param,sp='full',type='LU'):
    
    # get quad points and weights
    x = uqtkarray.dblArray2D()
    w = uqtkarray.dblArray1D()
    

    #print 'Create an instance of Quad class'
    q = uqtkquad.Quad(type,sp,ndim,param)

    #print 'Now set and get the quadrature rule...'
    q.SetRule()
    q.GetRule(x,w)

    # print out x and w
    #print 'Displaying the quadrature points and weights:\n'
    #print x
    #print w
    n = len(x)

    # get quad points
    x_np = zeros((n,ndim))
    x.getnpdblArray(x_np)
    
    # get quad weights
    w_np = zeros(n)
    w.getnpdblArray(w_np)
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
    





