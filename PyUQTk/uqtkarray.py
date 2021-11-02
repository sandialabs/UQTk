# Import the low-level C/C++ module
#if __package__ or "." in __name__:
#    from . import _uqtkarray
#else:
#    import _uqtkarray
import sys
sys.path.append('pyuqtkarray')

from _uqtkarray import *
from pyuqtkarray_tools import *
