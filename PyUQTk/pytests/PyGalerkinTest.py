#=====================================================================================
#
#                      The UQ Toolkit (UQTk) version 3.1.5
#                          Copyright (2024) NTESS
#                        https://www.sandia.gov/UQToolkit/
#                        https://github.com/sandialabs/UQTk
#
#     Copyright 2024 National Technology & Engineering Solutions of Sandia, LLC (NTESS).
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
sys.path.append('../pce/')
sys.path.append('../')
sys.path.append('../PyPCE/')

try:
    import _pce as uqtkpce
except:
	print('PyUQTk pce module not found')

try:
    import pce_tools
except ImportError:
    print("PyUQTk pce_tools module not found")

try:
	import numpy as np
except ImportError:
	print("Need numpy to test PyUQTk")

'''
This file tests the Galerkin projection pyqutk routine
'''

# create instance of the PCset class
print('Create a PC object')
nord=5
ndim=3
pc_type='LU'
param=nord+1
poly=uqtkpce.PCSet("NISPnoq", nord, ndim, pc_type, 0, 1)

# define true coefficients as a sequential series of numbers
npce=poly.GetNumberPCTerms()
coef=np.array([i+1 for i in range(npce)])

# obtain quadrature points
print('Obtain quadrature points')
poly.SetQuadRule(pc_type, 'full', param)
qdpts, totquat = pce_tools.UQTkGetQuadPoints(poly)

# evaluate the polynomial at the quadrature points
print('Evaluate the polynomial at quadrature points')
f_evals=pce_tools.UQTkEvaluatePCE(poly, coef, qdpts)

# use galerkin projection to find the coefficients
print('Find coefficients with Galerkin projection')
c_k = pce_tools.UQTkGalerkinProjection(poly,f_evals)

# print out coefficients
print('Displaying coefficients')
print(c_k)

# asserting the coefficients are correct
dec_place=5
reg=(np.round(c_k, dec_place)==coef)
assert (all(reg) and reg[0]==True)
