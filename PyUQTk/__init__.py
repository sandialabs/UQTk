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
# swig interface modules (only compiled if PyUQTK=On)
try:
	import uqtkarray
except ImportError:
	print "PyUQTk SWIG array interface not created."

try:
	import quad
except ImportError:
	print "PyUQTk SWIG quad interface not created."

try:
	import tools
except ImportError:
	print "PyUQTk SWIG tools interface not created."

try:
	import kle
except ImportError:
	print "PyUQTk SWIG kle interface not created."

try:
	import pce
except ImportError:
	print "PyUQTk SWIG pce interface not created."

try:
	import bcs
except ImportError:
	print "PyUQTk SWIG bcs interface not created."

try:
	import mcmc
except ImportError:
	print "PyUQTk SWIG mcmc interface not created."

try:
	import dfi
except:
	print "PyUQTk SWIG dfi interface not created."

# pure python tools (always included)
try:
	import inference
	import plotting
	import sens
except:
	print "Scipy and/or matplotlib may need to be installed"

import utils
import multirun
