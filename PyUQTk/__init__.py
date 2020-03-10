#=====================================================================================
#
#                      The UQ Toolkit (UQTk) version 3.1.0
#                          Copyright (2020) NTESS
#                        https://www.sandia.gov/UQToolkit/
#                        https://github.com/sandialabs/UQTk
#
#     Copyright 2020 National Technology & Engineering Solutions of Sandia, LLC (NTESS).
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
# swig interface modules (only compiled if PyUQTK=On)
try:
	from . import uqtkarray
except ImportError:
	print('PyUQTk SWIG array interface not created.')

try:
	from . import quad
except ImportError:
	print('PyUQTk SWIG quad interface not created.')

try:
	from . import tools
except ImportError:
	print('PyUQTk SWIG tools interface not created.')

try:
	from . import kle
except ImportError:
	print('PyUQTk SWIG kle interface not created.')

try:
	from . import pce
	from . import PyPCE
except ImportError:
	print('PyUQTk SWIG pce interface and PyPCE module not available.')

try:
	from . import bcs
except ImportError:
	print('PyUQTk SWIG bcs interface not created.')

try:
	from . import mcmc
except ImportError:
	print('PyUQTk SWIG mcmc interface not created.')

# pure python tools (always included)
try:
	from . import inference
	from . import plotting
	from . import sens
except:
	print('Scipy and/or matplotlib may need to be installed')


from . import utils
from . import multirun
