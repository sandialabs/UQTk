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
sys.path.append('../pyuqtkarray/')
sys.path.append('../pyuqtkarray_tools/')
sys.path.append('../')


try:
	import numpy as np
except ImportError:
	print("Need numpy to test PyUQTk")

try:
	import uqtkarray
except ImportError:
	print("PyUQTk array module not found")
	print("If installing in a directory other than the build directory, make sure PYTHONPATH includes the install directory")

'''
This file tests to make sure conversion from numpy -> uqtkarray does
not change the row-major (C contiguous) format of the original numpy array

Also, when converting form uqtkarray-> numpy we want to make sure that the
resulting numpy array is *only* row major (C contiguous)

'''

# create numpy matrix and show flags
a_np = np.array([[0, 2.00],[0.1, 1],[1, 5.0]])
print("flags for a_np to show whether C or F contiguous")
print(a_np.flags)


# get a uqtk array from a numpy array (memory is copied, not shared)
a_uqtk = uqtkarray.numpy2uqtk(a_np)
print("\nflags for original numpy array to make sure it hasn't changed to F continguous after converting")
# verify that the original numpy array is only C contiguous
assert a_np.flags['F_CONTIGUOUS'] == False
assert a_np.flags['C_CONTIGUOUS'] == True

print("\nConvert uqtk array back to numpy array and make sure C contiguous")
b_np = uqtkarray.uqtk2numpy(a_uqtk)
# test to make sure new numpy array is *only* C contiguous (row - major)
assert b_np.flags['F_CONTIGUOUS'] == False
assert b_np.flags['C_CONTIGUOUS'] == True

# test for the dot product
print("\ncompute dot product which should be [2,1.1,6] (Note that if F contigous, the dot product would be [.1,3,6]:")
dp = np.dot(b_np,np.ones(2))
assert np.all( dp ==  np.array([2.,1.1,6.]))
