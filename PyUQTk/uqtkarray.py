# Import the low-level C/C++ module
#if __package__ or "." in __name__:
#    from . import _uqtkarray
#else:
#    import _uqtkarray
import sys
sys.path.append('pyuqtkarray')

from _uqtkarray import *

# Import the _uqtkarray_tools to go from uqtk to numpy and vice versa
import numpy as np

def uqtk2numpy(x):
    if x.type() == 'int':
        s = x.shape()
        imin = np.argmin(s)
        if len(s) == 1:
            n = s[0]
            y = np.zeros(n,dtype='int64')
            y = getnpintArray(x)
        if len(s) == 2:
            #n = s[0]
            #m = s[1]
            #z = np.zeros((n,m),dtype='int64')
            #z = _uqtkarray.getnpintArray(x)
            #y = fixer(z)
            list = getnpintArray(x)
            m = x.XSize();
            n = x.YSize();
            y = np.full((m,n),0,dtype=int)
            counter = 0
            for i in range(m):
                for j in range(n):
                    y[i,j]=list[counter]
                    counter = counter + 1
#        if len(s) == 2 and np.amin(s) == 1:
#            y = np.array(x.flatten())
#            y = y[...,None]
        return y.copy()
    else:
        s = x.shape()
        imin = np.argmin(s)
        if len(s) == 1:
            n = s[0]
            y = np.zeros(n)
            y = getnpdblArray(x)
        if len(s) == 2:
            list = getnpdblArray(x)
            m = x.XSize();
            n = x.YSize();
            y = np.full((m,n),0,dtype=float)
            counter = 0
            for i in range(m):
                for j in range(n):
                    y[i,j]=list[counter]
                    counter = counter + 1
#        if len(s) == 2 and np.amin(s) == 1:
#            y = np.array(x.flatten())
#            y = y[...,None]
        return y.copy()
