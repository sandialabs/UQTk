#=====================================================================================
#                     The UQ Toolkit (UQTk) version @UQTKVERSION@
#                     Copyright (@UQTKYEAR@) Sandia Corporation
#                     http://www.sandia.gov/UQToolkit/
#
#     Copyright (@UQTKYEAR@) Sandia Corporation. Under the terms of Contract DE-AC04-94AL85000
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
from __future__ import print_function # To make print() in Python 2 behave like in Python 3

# include path to include PyUQTk
import sys
sys.path.append('../../')

try:
	from numpy import *
	from matplotlib.pyplot import *
except ImportError:
	print("Need numpy and matplotlib to test PyUQTk")

try:
	import PyUQTk.array as uqtkarray
	import PyUQTk.quad as uqtkquad
	import PyUQTk.pce as uqtkpce
except ImportError:
	print("PyUQTk array and quad module not found")


# get quad points and weights
x = uqtkarray.dblArray2D()
w = uqtkarray.dblArray1D()

print('Create an instance of Quad class')
ndim = 1
level = 8
q = uqtkquad.Quad('LU','full',ndim,level)
print('Now set and get the quadrature rule...')
q.SetRule()
q.GetRule(x,w)

# print out x and w
print('Displaying the quadrature points and weights:\n')
# print(x)
# print(w)
n = len(x)
print('Number of quad points is ', n, '\n')

# conver to numpy arrays
x_np = zeros((len(x),1))
w_np = zeros(len(x))
x.getnpdblArray(x_np)
w.getnpdblArray(w_np)

# define function for evaluation over [-1,1]
f = lambda x: 1./(1 + x**2)
y = f(x_np)
y.shape = (len(y),) # make 1d array

# convert numpy y to 1d array
ydata = uqtkarray.dblArray1D(len(y),0)
ydata.setnpdblArray(y)

'''
Define PCSet object
'''
# Instantiate object
nord = 8
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

'''
Evaluate PC Model at random points
'''
xall = linspace(-1,1,1000); xall.shape = (len(xall),1)
yeval = uqtkarray.dblArray1D(len(xall),0.0)
xeval = uqtkarray.dblArray2D(len(xall),1,0.0)
xeval.setnpdblArray(asfortranarray(xall))
pcmodel.EvalPCAtCustPoints(yeval,xeval,ck)

y_exact = f(xall)
y_pce = array(yeval.flatten())
plot(xall,y_exact,'k',lw=2,alpha=.4)
plot(xeval,y_pce,'--r',lw=1)

'''
Evaluate PC Model at quad points
'''
yevalq = uqtkarray.dblArray1D(len(x),0.0)
pcmodel.EvalPCAtCustPoints(yevalq,x,ck)
y_pceq = array(yevalq.flatten())
plot(x,y_pceq,'or',lw=1)

'''
Plot the ImportError
'''
figure()
plot(xeval,abs(y_pce - y_exact[:,0]),'k')
