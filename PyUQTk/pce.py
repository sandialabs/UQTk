# Import the low-level C/C++ module
'''if __package__ or "." in __name__:
    from . import _pce
else:
    import _pce
'''
import sys
sys.path.append('pce')
from _pce import *
