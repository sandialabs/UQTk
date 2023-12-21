#=====================================================================================
#
#                      The UQ Toolkit (UQTk) version 3.1.4
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
#     Questions? Contact the UQTk Developers at https://github.com/sandialabs/UQTk/discussions
#     Sandia National Laboratories, Livermore, CA, USA
#=====================================================================================
from __future__ import print_function # To make print() in Python 2 behave like in Python 3

# include path to include PyUQTk
import sys
sys.path.append('../../')

try:
	from numpy import *
except ImportError:
	print("NumPy needed to test PyUQTk.")
try:
	from matplotlib.pyplot import *
except ImportError:
	print("Matplotlib needed to test PyUQTk")

try:
	import PyUQTk.array as uqtkarray
except ImportError:
	print("PyUQTk array module not found")
try:
	import PyUQTk.quad as uqtkquad
except ImportError:
	print("PyUQTk quad module not found")
try:
	import PyUQTk.pce as uqtkpce
except ImportError:
	print("PyUQTk pce module not found")


# get quad points and weights
x = uqtkarray.dblArray2D()
w = uqtkarray.dblArray1D()

print('Create an instance of Quad class')
ndim = 2
level = 8
q = uqtkquad.Quad('LU','full',ndim,level)
print('Now set and get the quadrature rule...')
q.SetRule()
q.GetRule(x,w)

# print out x and w
print('Displaying the quadrature points and weights:\n')
print(x)
print(w)
n = len(x)
print('Number of quad points is ', n, '\n')
# conver to numpy arrays
x_np = zeros((n,2))
w_np = zeros(len(x))
x.getnpdblArray(x_np)
w.getnpdblArray(w_np)

# define function for evaluation over [-1,1]
f = lambda x: x[:,0]*x[:,1] + x[:,0]**2 + sqrt(abs(x[:,1]))
y = f(x_np)
y.shape = (len(y),) # make 1d array

# convert numpy y to 1d array
ydata = uqtkarray.dblArray1D(len(y),0)
ydata.setnpdblArray(y)

'''
Define PCSet object
'''
# Instantiate object
nord = 4
chaos_type = "LEG"
pcmodel = uqtkpce.PCSet('NISPnoq',nord,ndim,'LEG')

# set quad rule for pc model
pcmodel.SetQuadRule(q)
nup = pcmodel.GetNumberPCTerms()-1
totquad = pcmodel.GetNQuadPoints()

# Get the multiindex for postprocessing
mindex = uqtkarray.intArray2D();
pcmodel.GetMultiIndex(mindex);

# get the coefficients using the quadrature rule
# to calculate the projections
ck = uqtkarray.dblArray1D(nup+1,0.0)
pcmodel.GalerkProjection(ydata,ck);
c_np = zeros(len(ck))
ck.getnpdblArray(c_np)

# compute main sensitivities
mainsens = uqtkarray.dblArray1D(ndim,0)
pcmodel.ComputeMainSens(ck,mainsens)

# compute total sensitivity
totsens = uqtkarray.dblArray1D(ndim,0)
pcmodel.ComputeTotSens(ck,totsens)

#compute joint sensitivity
jointsens = uqtkarray.dblArray2D(ndim,ndim,0)
pcmodel.ComputeJointSens(ck,jointsens)

print(mainsens, totsens, jointsens)

# '''
# Evaluate PC Model at random points
# '''
# xall = linspace(-1,1,1000); xall.shape = (len(xall),1)
# yeval = uqtkarray.dblArray1D(len(xall),0.0)
# xeval = uqtkarray.dblArray2D(len(xall),1,0.0)
# xeval.setnpdblArray(asfortranarray(xall))
# pcmodel.EvalPCAtCustPoints(yeval,xeval,ck)

# y_exact = f(xall)
# y_pce = array(yeval.flatten())
# plot(xall,y_exact,'k',lw=2,alpha=.4)
# plot(xeval,y_pce,'--r',lw=1)

# '''
# Evaluate PC Model at quad points
# '''
# yevalq = uqtkarray.dblArray1D(len(x),0.0)
# pcmodel.EvalPCAtCustPoints(yevalq,x,ck)
# y_pceq = array(yevalq.flatten())
# plot(x,y_pceq,'or',lw=1)

# '''
# Plot the ImportError
# '''
# figure()
# plot(xeval,abs(y_pce - y_exact[:,0]),'k')
