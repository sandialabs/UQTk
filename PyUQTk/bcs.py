# Import the low-level C/C++ module
if __package__ or "." in __name__:
    from . import _bcs
else:
    import _bcs
# Import Python Utility for PyBCS
from bcs_ext import *
