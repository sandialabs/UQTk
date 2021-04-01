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

import numpy as np

try:
	import pyuqtkarray
except ImportError:
	print("PyUQTk array module not found")
	print("If installing in a directory other than the build directory, make sure PYTHONPATH includes the install directory")

def uqtk2numpy(x):
	if x.type() == 'int':
		s = x.shape()
		imin = np.argmin(s)
		if len(s) == 1:
			n = s[0]
			y = np.zeros(n,dtype='int64')
			pyuqtkarray.getnpintArray(x,y)
		if len(s) == 2 and np.amin(s) > 1:
			n = s[0]
			m = s[1]
			y = np.zeros((n,m),dtype='int64')
			pyuqtkarray.getnpintArray(x,y)
		if len(s) == 2 and np.amin(s) == 1:
			y = np.array(x.flatten())
		return y.copy()
	else:
		s = x.shape()
		imin = np.argmin(s)
		if len(s) == 1:
			n = s[0]
			y = np.zeros(n)
			pyuqtkarray.getnpdblArray(x,y)
		if len(s) == 2 and np.amin(s) > 1:
			n = s[0]
			m = s[1]
			y = np.zeros((n,m))
			pyuqtkarray.getnpdblArray(x,y)
		if len(s) == 2 and np.amin(s) == 1:
			y = np.array(x.flatten())
		return y.copy()

def numpy2uqtk(y):
	s = np.shape(y)
	if len(s) == 1:
		n = s[0]
		x = pyuqtkarray.dblArray1D(n)
	if len(s) == 2:
		n = s[0]
		m = s[1]
		x = pyuqtkarray.dblArray2D(n,m)
	#x.setnpdblArray(np.asfortranarray(y.copy()))
	pyuqtkarray.setnpdblArray(x,np.asfortranarray(y.copy()))
	return x
