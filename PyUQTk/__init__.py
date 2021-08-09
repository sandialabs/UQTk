#=====================================================================================
#
#                      The UQ Toolkit (UQTk) version 3.1.1
#                          Copyright (2021) NTESS
#                        https://www.sandia.gov/UQToolkit/
#                        https://github.com/sandialabs/UQTk
#
#     Copyright 2021 National Technology & Engineering Solutions of Sandia, LLC (NTESS).
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
except ImportError as ie:
    print('PyUQTk SWIG uqtkarray interface not imported.', ie)

try:
    from . import quad
except ImportError as ie:
    print('PyUQTk SWIG quad interface not imported.', ie)

try:
    from . import tools
except ImportError as ie:
    print('PyUQTk SWIG tools interface not imported.', ie)

try:
    from . import kle
except ImportError as ie:
    print('PyUQTk SWIG kle interface not imported.', ie)

try:
    from . import pce
except ImportError as ie:
    print('PyUQTk SWIG pce interface not imported.', ie)

try:
    from . import PyPCE
except ImportError as ie:
    print('PyUQTk SWIG PyPCE module not imported.', ie)

try:
    from . import bcs
except ImportError as ie:
    print('PyUQTk SWIG bcs interface not imported.', ie)

try:
    from . import mcmc
except ImportError as ie:
    print('PyUQTk SWIG mcmc interface not imported.', ie)

# pure python tools (always included)
try:
    from . import inference
except ImportError as ie:
    print('Scipy and/or matplotlib may need to be installed', ie)

try:
    from . import plotting
except ImportError as ie:
    print('Scipy and/or matplotlib may need to be installed', ie)

try:
    from . import sens
except ImportError as ie:
    print('Scipy and/or matplotlib may need to be installed', ie)

from . import utils
from . import multirun
