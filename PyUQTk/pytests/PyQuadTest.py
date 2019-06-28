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
# include path to include PyUQTk
import sys
sys.path.append('../uqtkarray/')
sys.path.append('../quad/')

try:
	from numpy import *
	from matplotlib.pyplot import *
except ImportError:
	"Need numpy and matplotlib to test PyUQTk"

try:
	import uqtkarray
	from uqtkarray import numpy2uqtk
	from uqtkarray import uqtk2numpy
	import quad as uqtkquad
except ImportError:
	print "PyUQTk array and quad module not found"

'''
This file tests the quadrature pyqutk routine
'''

# true quad points for sparse LU with ndim = 2 and level = 3
qpnts_ref = array([[-9.681602395076263079e-01, 0.000000000000000000e+00],
[-9.061798459386638527e-01, -7.745966692414832933e-01],
[-9.061798459386638527e-01, 0.000000000000000000e+00],
[-9.061798459386638527e-01, 7.745966692414834043e-01],
[-8.360311073266358806e-01, 0.000000000000000000e+00],
[-7.745966692414832933e-01, -9.061798459386638527e-01],
[-7.745966692414832933e-01, -7.745966692414832933e-01],
[-7.745966692414832933e-01, -5.384693101056832187e-01],
[-7.745966692414832933e-01, 0.000000000000000000e+00],
[-7.745966692414832933e-01, 5.384693101056829967e-01],
[-7.745966692414832933e-01, 7.745966692414834043e-01],
[-7.745966692414832933e-01, 9.061798459386638527e-01],
[-6.133714327005902467e-01, 0.000000000000000000e+00],
[-5.384693101056832187e-01, -7.745966692414832933e-01],
[-5.384693101056832187e-01, 0.000000000000000000e+00],
[-5.384693101056832187e-01, 7.745966692414834043e-01],
[-3.242534234038090268e-01, 0.000000000000000000e+00],
[0.000000000000000000e+00, -9.681602395076263079e-01],
[0.000000000000000000e+00, -9.061798459386638527e-01],
[0.000000000000000000e+00, -8.360311073266358806e-01],
[0.000000000000000000e+00, -7.745966692414832933e-01],
[0.000000000000000000e+00, -6.133714327005902467e-01],
[0.000000000000000000e+00, -5.384693101056832187e-01],
[0.000000000000000000e+00, -3.242534234038090268e-01],
[0.000000000000000000e+00, 0.000000000000000000e+00],
[0.000000000000000000e+00, 3.242534234038088048e-01],
[0.000000000000000000e+00, 5.384693101056829967e-01],
[0.000000000000000000e+00, 6.133714327005905798e-01],
[0.000000000000000000e+00, 7.745966692414834043e-01],
[0.000000000000000000e+00, 8.360311073266353254e-01],
[0.000000000000000000e+00, 9.061798459386638527e-01],
[0.000000000000000000e+00, 9.681602395076263079e-01],
[3.242534234038088048e-01, 0.000000000000000000e+00],
[5.384693101056829967e-01, -7.745966692414832933e-01],
[5.384693101056829967e-01, 0.000000000000000000e+00],
[5.384693101056829967e-01, 7.745966692414834043e-01],
[6.133714327005905798e-01, 0.000000000000000000e+00],
[7.745966692414834043e-01, -9.061798459386638527e-01],
[7.745966692414834043e-01, -7.745966692414832933e-01],
[7.745966692414834043e-01, -5.384693101056832187e-01],
[7.745966692414834043e-01, 0.000000000000000000e+00],
[7.745966692414834043e-01, 5.384693101056829967e-01],
[7.745966692414834043e-01, 7.745966692414834043e-01],
[7.745966692414834043e-01, 9.061798459386638527e-01],
[8.360311073266353254e-01, 0.000000000000000000e+00],
[9.061798459386638527e-01, -7.745966692414832933e-01],
[9.061798459386638527e-01, 0.000000000000000000e+00],
[9.061798459386638527e-01, 7.745966692414834043e-01],
[9.681602395076263079e-01, 0.000000000000000000e+00]])

# initiate uqtk arrays for quad points and weights
x = uqtkarray.dblArray2D()
w = uqtkarray.dblArray1D()

# create instance of quad class and output
# points and weights
print 'Create an instance of Quad class'
ndim = 2
level = 3
q = uqtkquad.Quad('LU','sparse',ndim,level,0,1)
print 'Now set and get the quadrature rule...'
q.SetRule()
q.GetRule(x,w)

# print out x and w
print 'Displaying the quadrature points and weights:\n'
x_np = uqtk2numpy(x)
print x_np
n = len(x)
print 'Number of quad points is ', n, '\n'

# plot the quadrature points
print 'Plotting the points (get points in column major order as a flattened vector)'
print 'need to use reshape with fortran ordering'
xpnts = zeros((n,ndim))
x.getnpdblArray(xpnts)
# plot(xpnts[:,0], xpnts[:,1],'ob',ms=10,alpha=.25)
# show()

# convert the quad weights to numpy arrays
w_np = zeros(n)
w.getnpdblArray(w_np)

# asserting the quadrature points are correct
m,n = x_np.shape
for i in range(m):
	for j in range(n):
		assert x_np[i,j] == qpnts_ref[i,j]

