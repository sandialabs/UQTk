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
# include path for PyUQTk.
import sys
sys.path.append('../bcs/') # imports as build lib so installing not needed
sys.path.append('../uqtkarray/')
sys.path.append('../tools/')
sys.path.append('../pce/')

import os 
dir_path = os.path.dirname(os.path.realpath(__file__))
print dir_path

try:
	import numpy as np
	import matplotlib.pyplot as mpl
	import pdb
except ImportError:
	"Need numpy and matplotlib to test PyUQTk"

try:
	import uqtkarray as uqtkarray
	import pce as uqtkpce
	import tools as uqtktools
	from bcs import bcsreg
except ImportError:
	print "PyUQTk array and quad module not found"

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
regmodel = bcsreg(ndim=2,pcorder=pcorder,pctype="LU")
c, mindex = regmodel.fit(X,y,upit=upit,tol=tol)
ypred = regmodel.predict(Xtest)

# print mean squared prediction error
mse = np.mean((ypred - ytest)**2)
nmse = np.mean((ypred - ytest)**2)/np.mean(ytest)
print "\nMSE is %.5g" %mse
print "NMSE is %.5g" %nmse

# print sensitivities
print "\nSensitivities are ", regmodel.getsens()

prec = 1e-7
assert mse < prec, "BCS failed to recover the coefficients to desired precision :-("



