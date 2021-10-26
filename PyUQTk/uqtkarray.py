# Import the low-level C/C++ module
if __package__ or "." in __name__:
    from . import _uqtkarray
else:
    import _uqtkarray

# Import the pyuqtkarray_tools to go from uqtk to numpy and vice versa
from pyuqtkarray_tools import *
