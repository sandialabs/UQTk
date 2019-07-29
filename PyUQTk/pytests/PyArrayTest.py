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
# only necessary for cmake tests, so that user doesn have to "make install" to run
# python tests
import sys
sys.path.append('../uqtkarray/')

# try to import numpy and matplotlib
try:
	from numpy import *
	from matplotlib.pyplot import *
except ImportError:
	print("Need numpy and matplotlib to test PyUQTk")

# try to import uqtk array library and
# functions to convert between uqtk and numpy arrays
try:
	import uqtkarray
	from uqtkarray import numpy2uqtk
	from uqtkarray import uqtk2numpy
except ImportError:
	print("PyUQTk array module not found")
	print("If installing in a directory other than the build directory, make sure PYTHONPATH includes the install directory")

import unittest

''' Test converting 1d numpy array to 1d uqtk array '''
# create 1d array
N = 35
x = uqtkarray.dblArray1D(N,0)

# create 1d numpy array
x_np = random.randn(N)

# set uqtk array to numpy array
x.setnpdblArray(x_np)

# test to make sure array elements are the same
for i in range(N):
	assert x[i] == x_np[i]

''' Test converting 2d numpy array to 2d uqtk array '''
# create 2d array in uqtk
m = 100
n = 3
y = uqtkarray.dblArray2D(m,n,1)

# set 2d array to numpy array
# make sure to pass asfortranarray
y_np = random.randn(m,n)
y.setnpdblArray(asfortranarray(y_np))

for i in range(m):
	for j in range(n):
		assert y[i,j] == y_np[i,j]

''' alternative using uqtk2numpy and numpy2uqtk '''

# test conversion from 1d numpy array to 1d uqtk array
nn = 10
x1 = random.rand(nn)
y1 = numpy2uqtk(x1)
z1 = uqtk2numpy(y1)
for i in range(nn):
	assert x1[i] == y1[i]

# test conversion from 1d uqtk array to numpy
for i in range(nn):
	assert z1[i] == x1[i]

# test for conversion from 2d numpy to 2d uqtk
nn = 10
mm = 5
X1 = random.rand(mm,nn)
Y1 = numpy2uqtk(X1)
Z1 = uqtk2numpy(Y1)
for i in range(mm):
	for j in range(nn):
		assert X1[i,j] == Y1[i,j]

# test for conversion from 2d uqtk array to numpy array
for i in range(mm):
	for j in range(nn):
		assert Z1[i,j] == X1[i,j]
print('done')
