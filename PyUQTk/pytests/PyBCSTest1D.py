#=====================================================================================
#
#                      The UQ Toolkit (UQTk) version 3.1.3
#                          Copyright (2023) NTESS
#                        https://www.sandia.gov/UQToolkit/
#                        https://github.com/sandialabs/UQTk
#
#     Copyright 2023 National Technology & Engineering Solutions of Sandia, LLC (NTESS).
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
from __future__ import print_function # To make print() in Python 2 behave like in Python 3

import os
src = os.getenv('UQTK_SRC')

# include path for PyUQTk.
import sys
sys.path.append('../bcs/') # imports as build lib so installing not needed
sys.path.append('../pyuqtkarray/')
sys.path.append('../tools/')
sys.path.append('../pce/')
sys.path.append('../pyuqtkarray_tools/')
sys.path.append('../bcs_ext/')
sys.path.append('../quad')

import os
dir_path = os.path.dirname(os.path.realpath(__file__))
print(dir_path)

try:
	import numpy as np
except ImportError:
	print("Need numpy to test PyUQTk")
try:
	import matplotlib.pyplot as mpl
except ImportError:
	print("Need matplotlib to test PyUQTk")
try:
	import pdb
except ImportError:
	print("Need pdb module to test PyUQTk")

try:
    import _uqtkarray as uqtkarray
except ImportError:
	print("PyUQTk array module not found")
try:
    import pyuqtkarray_tools
except ImportError:
	print("pyuqtkarray_tools module not found")

try:
    import _quad as quad
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

try:
	import _bcs as bcs
	import bcs_ext as bcsTools
except ImportError:
    print("BCS Regression module not found")


nord = 8
ndim = 1
level = 16
q = quad.Quad("LU","full",ndim,level,0.0,1.0)

x = uqtkarray.dblArray2D()
w = uqtkarray.dblArray1D()
index = uqtkarray.intArray2D()
q.SetRule()
q.GetRule(x,w,index)

pcmodel = uqtkpce.PCSet("NISPnoq",nord,ndim,"LEG",0.0,1.0)

Phi = uqtkarray.dblArray2D()
pcmodel.EvalBasisAtCustPts(x,Phi)

y = uqtkarray.dblArray1D(x.XSize(),0.0)
for i in range(x.XSize()):
	value = 1.0/(1.0 + x.at(i,0)*x.at(i,0));
	y.assign(i,value)

#sigma = uqtkarray.dblArray1D(1,1e-8)
sigma = 1e-8
eta = 1e-8
lambda_init = uqtkarray.dblArray1D()
scale = 0.1

weights = uqtkarray.dblArray1D()
errbars = uqtkarray.dblArray1D()
basis = uqtkarray.dblArray1D()
alpha = uqtkarray.dblArray1D()
used = uqtkarray.intArray1D()
#_lambda = uqtkarray.dblArray1D(1,0.0)
_lambda=0.0

adaptive = 1
optimal = 1
verbose = 0

bcs.BCS(Phi,y,sigma,eta,lambda_init,adaptive,optimal,scale,verbose,weights,used,errbars,basis,alpha,_lambda)

uqtkarray.printarray(weights)
uqtkarray.printarray(used)
uqtkarray.printarray(errbars)

print(abs(0.785397 - weights[0]))
print(abs(-0.353987 - weights[1]))
print(abs(0.00316971 - weights[4]))
print(abs(used[1] - 2))
print(abs(used[2] - 4))
print(abs(used[4] - 8))

assert abs(0.785397 - weights[0]) < 1e-6
assert abs(-0.353987 - weights[1]) < 1e-6
assert abs(0.00316971 - weights[4]) < 1e-6
assert abs(used[1] - 2) < 1e-16
assert abs(used[2] - 4) < 1e-16
assert abs(used[4] - 8) < 1e-16
