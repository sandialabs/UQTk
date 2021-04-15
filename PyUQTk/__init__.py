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
# swig interface modules (only compiled if PyUQTK=On)
try:
    from pyuqtkarray import *
    from pyuqtkarray_tools import *
except ImportError:
    print('PyUQTk SWIG array interface not created.')

try:
    from quad import *
except ImportError:
    print('PyUQTk SWIG quad interface not created.')

try:
    from tools import *
except ImportError:
    print('PyUQTk SWIG tools interface not created.')

try:
    from kle import *
except ImportError:
    print('PyUQTk SWIG kle interface not created.')

try:
    from pce import *
    from adaptation_tools import *
    from pce_tools import *
except ImportError:
    print('PyUQTk SWIG pce interface and PyPCE module not available.')

try:
    from bcs import *
    from bcs_ext import *
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
