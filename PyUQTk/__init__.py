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
# Set up paths
import sys
import os
# Append the folder where this init script is located to the search path
# sys.path.append(os.environ['UQTK_INS'] + '/PyUQTk')
sys.path.append(os.path.realpath(os.path.dirname(__file__)))

# pure python tools (always included)
try:
    from . import plotting
except ImportError as ie:
    print('Some modules in the plotting package did not load')

try:
    from . import sens
except ImportError as ie:
    print('Some modules in the sensitivity package did not load')

try:
    from . import inference
except ImportError as ie:
    print('Some modules in the inference package did not load')

from . import utils
from . import multirun

try:
    import bcs
except:
    print('PyBCS module not imported')

try:
    import kle
except:
    print('PyKLE module not imported')

try:
    import pce
except:
    print('PyPCE module not imported')

try:
    import quad
except:
    print('PyQuad module not imported')

try:
    import tools
except:
    print('PyTools module not imported')

try:
    import uqtkarray
except:
    print('PyUQTkArray module not imported')
