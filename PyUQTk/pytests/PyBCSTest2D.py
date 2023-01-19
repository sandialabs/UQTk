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
    import _quad as uqtkquad
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

ck = uqtkarray.dblArray1D(10,0.0)
mindex = uqtkarray.intArray2D(10,2)

ck.assign(0,0.666666666666664)
ck.assign(1,1.600000000000499)
ck.assign(2,1.000000000000289)
ck.assign(5,-0.6666666666668039)
ck.assign(6,0.4000000000008473)

mindex.assign(0,0,0)
mindex.assign(0,1,0)
mindex.assign(1,0,1)
mindex.assign(1,1,0)
mindex.assign(2,0,0)
mindex.assign(2,1,1)
mindex.assign(3,0,2)
mindex.assign(3,1,0)
mindex.assign(4,0,1)
mindex.assign(4,1,1)
mindex.assign(5,0,0)
mindex.assign(5,1,2)
mindex.assign(6,0,3)
mindex.assign(6,1,0)
mindex.assign(7,0,2)
mindex.assign(7,1,1)
mindex.assign(8,0,1)
mindex.assign(8,1,2)
mindex.assign(9,0,0)
mindex.assign(9,1,3)

pcmodel = uqtkpce.PCSet("NISPnoq",mindex,"LEG",0.0,1.0)

q = uqtkquad.Quad("LU","sparse",2,5,0.0,1.0)
x = uqtkarray.dblArray2D()
w = uqtkarray.dblArray1D()
q.SetRule()
q.GetRule(x,w)

y = uqtkarray.dblArray1D(x.XSize(),0.0)
pcmodel.EvalPCAtCustPoints(y,x,ck)

Phi = uqtkarray.dblArray2D()
pcmodel.EvalBasisAtCustPts(x,Phi)

#sigma = uqtkarray.dblArray1D(1,1e-8)
sigma = 1e-8
eta = 1e-12
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
uqtkarray.printarray(ck)

assert used[0] == 1
assert used[1] == 2
assert used[2] == 0

assert abs(weights[0] - 1.6) < 1e-8
assert abs(weights[1] - 1) < 1e-8
assert abs(weights[2] - float(2/3)) < 1e-8
