#!/usr/bin/env python

import sys
import numpy as npy

ndim=5
a=[0.001,0.01,0.01,0.1,1.0]

if ( len(sys.argv) != ndim+1 ):
  quit()

x=[float(sys.argv[i]) for i in range(1,ndim+1)]

ax2 = [a[i]*x[i]*x[i] for i in range(ndim)]
rval = npy.exp(-sum(ax2))

f = open("rval.txt", "w")
f.write( str(rval)  ) 
f.close()

quit()

