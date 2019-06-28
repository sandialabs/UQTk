#!/usr/bin/env python
#=====================================================================================
#                     The UQ Toolkit (UQTk) version 3.0.4
#                     Copyright (2017) Sandia Corporation
#                     http://www.sandia.gov/UQToolkit/
#
#     Copyright (2017) Sandia Corporation. Under the terms of Contract DE-AC04-94AL85000
#     with Sandia Corporation, the U.S. Government retains certain rights in this software.
#
#     This file is part of The UQ Toolkit (UQTk)
#
#     UQTk is free software: you can redistribute it and/or modify
#     it under the terms of the GNU Lesser General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
#
#     UQTk is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU Lesser General Public License for more details.
#
#     You should have received a copy of the GNU Lesser General Public License
#     along with UQTk.  If not, see <http://www.gnu.org/licenses/>.
#
#     Questions? Contact Bert Debusschere <bjdebus@sandia.gov>
#     Sandia National Laboratories, Livermore, CA, USA
#===================================================================================== 

import os
import shutil
import sys
try:
    import numpy as np
except ImportError:
    print "Numpy was not found. "

try:
    import matplotlib
except ImportError:
    print "Matplotlib was not found. "
import math
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.mlab as mlab


sys.path.append(os.environ['UQTK_INS'])
import PyUQTk.plotting.surrogate as ss
import PyUQTk.plotting.inout as io

from pylab import *

import cPickle as pick

rc('legend',loc='best', fontsize=22)
rc('lines', linewidth=2, color='r')
rc('axes',linewidth=3,grid=True,labelsize=28)
rc('xtick',labelsize=20)
rc('ytick',labelsize=20)

    
#############################################################
#############################################################


plotid=sys.argv[1]

# Parallel coordinates
if(plotid=='pcoord'): 
    
    # Parse the arguments
    ndim=int(sys.argv[2])


    if os.path.exists('pnames.txt'):
        with open('pnames.txt') as f:
            pnames = f.read().splitlines()
            assert(len(pnames)==ndim)
    else:
        pnames=['Param # '+str(i) for i in range(1,ndim+1)]


    ndcut=int(sys.argv[3])
    if ndcut==0:
        ndcut=ndim

    ndg=int((ndim-1)/ndcut)+1

    inputfiles=sys.argv[4:]


    inputs=np.empty((0,ndim))
    labels=[]
    for inputfile in inputfiles:
        this_input=np.loadtxt(inputfile)
        inputs=np.vstack((inputs,this_input))
        labels.extend([inputfile]*this_input.shape[0])

    
    #labels= (inputs[:,iout]>0.0)#(out_ss==0.0) #(out<0.004)

        
    
    for i in range(ndg):
        print "Plotting %d / %d " % (i+1,ndg)
        names=range(1+i*ndcut,min(1+(i+1)*ndcut,ndim+1))
        
        values=inputs[:,i*ndcut:min((i+1)*ndcut,ndim)].T
        io.parallel_coordinates(names, values, labels,'pcoord_'+str(i+1)+'.eps')

        #labels_only=labels[labels==True]
        #values_only=values[:,labels==True]
        #io.parallel_coordinates(names, values_only, labels_only)

# Plot one input versus another
elif(plotid=='xx'):

    ndim=int(sys.argv[2])
    d1=int(sys.argv[3])
    d2=int(sys.argv[4])


    if os.path.exists('pnames.txt'):
        with open('pnames.txt') as f:
            pnames = f.read().splitlines()
            pname1=pnames[d1]
            pname2=pnames[d2]
    else:
        pname1,pname2=['Param # '+str(i) for i in (d1,d2)]

    print "Plotting %s vs %s" % (pname1,pname2)

    inputfiles=sys.argv[5:]

    inputs=np.empty((0,ndim))
    labels=[]
    for inputfile in inputfiles:
        this_input=np.loadtxt(inputfile)
        inputs=np.vstack((inputs,this_input))
        labels.extend([inputfile]*this_input.shape[0])

    io.plot_xx(d1,d2,pnames, inputs, labels,'xx_'+pname1+'_'+pname2+'.eps')

# Plot output vs input
elif(plotid=='xy'):

    d=int(sys.argv[2])
    o=int(sys.argv[3])


    if os.path.exists('pnames.txt'):
        with open('pnames.txt') as f:
            pnames = f.read().splitlines()
            pname=pnames[d]
    else:
        pname='Param # '+str(d) 

    if os.path.exists('outnames.txt'):
        with open('outnames.txt') as f:
            outnames = f.read().splitlines()
            outname=outnames[o]
    else:
        outname='Output # '+str(o) 



    print "Plotting %s vs %s" % (pname,outname)

    inputfile=sys.argv[4]
    outputfile=sys.argv[5]

    input=np.loadtxt(inputfile,ndmin=2)
    output=np.loadtxt(outputfile,ndmin=2)

    io.plot_xy(input[:,d],output[:,o],pname, outname,'xy_'+pname+'_'+outname+'.eps')
    
# Plot output vs input1,input2
elif(plotid=='xxy'):

    d1=int(sys.argv[2])
    d2=int(sys.argv[3])
    o=int(sys.argv[4])


    if os.path.exists('pnames.txt'):
        with open('pnames.txt') as f:
            pnames = f.read().splitlines()
    else:
        pnames=['Param # '+str(d) for d in [d1,d2]]

    if os.path.exists('outnames.txt'):
        with open('outnames.txt') as f:
            outnames = f.read().splitlines()
            outname=outnames[o]
    else:
        outname='Output # '+str(o) 



    print "Plotting %s,%s vs %s" % (pnames[0],pnames[1],outname)

    inputfile=sys.argv[5]
    outputfile=sys.argv[6]

    input=np.loadtxt(inputfile,ndmin=2)
    output=np.loadtxt(outputfile,ndmin=2)

    io.plot_xxy(input[:,d1],input[:,d2],output[:,o],[pnames[d1],pnames[d2]], outname,savefig='xxy_'+pnames[0]+'_'+pnames[1]+'_'+outname+'.eps')
    

else:
    print "plotid not recognized. Exiting."
    sys.exit()
