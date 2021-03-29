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

import os
dir_path = os.path.dirname(os.path.realpath(__file__))
print(dir_path)

try:
	import numpy as np
	import matplotlib.pyplot as mpl
	import pdb
except ImportError:
	print("Need numpy and matplotlib to test PyUQTk")

try:
    import pyuqtkarray as uqtkarray
    import pyuqtkarray_tools
except ImportError:
	print("PyUQTk array module not found")

try:
    import _pce as uqtkpce
    import _tools as uqtktools
except ImportError:
    print("PyUQTk PCE module not found")

try:
	import _bcs as bcs
	import bcs_ext as bcsTools
except ImportError:
    print("BCS Regression module not found")

'''
This example uses BCS to fit

f(x,y) = 1 + x + .5(3y^2-1)

using 100 randomly generating training data points
and 20 test data points. Sensitivity analysis is also
performed post fitting.

'''

# set dimension
ndim = 2

# Create training data
rn = np.random.RandomState(145)
X = 2*rn.rand(100,ndim) - 1
x1,x2 = X.T[0],X.T[1]
f = lambda x1,x2: 1 + x1 + .5*(3*x2**2-1)
y = f(x1,x2)

# create test data
Xtest = 2*rn.rand(20,ndim) - 1
ytest = f(Xtest.T[0],Xtest.T[1])
testdata = {'X': Xtest, 'y': ytest}

# BCS hyperparameter definitions
sigsq=None
pcorder = 2
pctype = "LU"
tol=1e-12
upit=1

# setup, git and predict bcs model
regmodel = bcsTools.bcsreg(ndim=2,pcorder=pcorder,pctype="LU")
err, coeff, mindex = regmodel.fit(X,y,upit=upit,tol=tol)
ypred = regmodel.predict(Xtest)

# print mean squared prediction error
mse = np.mean((ypred - ytest)**2)
nmse = np.mean((ypred - ytest)**2)/np.mean(ytest)
print("\nMSE is {:.5g}".format(mse))
print("NMSE is {:.5g}".format(nmse))

# print sensitivities
print("\nSensitivities are ", regmodel.getsens())

prec = 1e-7
assert mse < prec, "BCS failed to recover the coefficients to desired precision :-("
