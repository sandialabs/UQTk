#!/usr/bin/env python
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

import argparse
import os
import sys
try:
    import numpy as np
except ImportError:
    print("Numpy was not found. ")

try:
    import matplotlib
except ImportError:
    print("Matplotlib was not found. ")
import math

sys.path.append(os.environ['UQTK_INS'])
import PyUQTk.plotting.surrogate as ss
import PyUQTk.plotting.inout as io

from pylab import *


rc('legend',loc='best', fontsize=22)
rc('lines', linewidth=2, color='r')
rc('axes',linewidth=3,grid=True,labelsize=28)
rc('xtick',labelsize=20)
rc('ytick',labelsize=20)


#############################################################
#############################################################


## Parsing the inputs
usage_str='           Makes exploratory visualizations,usually before the surrogates are constructed.\n\
           First argument defines the plot type. \n\
           Try "plot_prep.py <plot_type> -h" for help on additional arguments under each plot_type. \n\
           When relevant, parameter names and output names can be read in pnames.txt and outnames.txt (if files not present, will use generic names).\n\
           Many options are quite experimental yet, and may not be optimal visually.'

parser = argparse.ArgumentParser(description=usage_str,formatter_class=argparse.RawTextHelpFormatter)
arg1_choices=['pcoord','xx','xy','xxy']
#parser.add_argument('plot_type', type=str,nargs=1,help="Plot type", choices=arg1_choices)
subparsers=parser.add_subparsers()
for arg1 in arg1_choices:
    if arg1=='pcoord' :
        sbparser = subparsers.add_parser(arg1,formatter_class=argparse.RawTextHelpFormatter,\
            description='Plots input parameters in parallel coordinates.', \
            epilog='Example: \n"plot_prep.py pcoord 2 0 qtrain.dat qval.dat"')
        sbparser.add_argument('ndim', type=int,nargs='?',help="The dimensionality of inputs")
        sbparser.add_argument('ndcut', type=int,nargs='?',help="Chunk size, if one wants to plot a parameter chunk at a time, useful for very high-d (if 0, plots all in one shot)")
        sbparser.add_argument('inputfiles', type=str,nargs='*',help="Parameter file names, can include more than one.")
    elif arg1=='xx' :
        sbparser = subparsers.add_parser(arg1,formatter_class=argparse.RawTextHelpFormatter,\
            description='Plots input[d1] versus input[d2], i.e. 2d projection of input samples.', \
            epilog='Example: \n"plot_prep.py xx 2 0 1 qtrain.dat qval.dat"')
        sbparser.add_argument('ndim', type=int,nargs='?',help="The dimensionality of inputs")
        sbparser.add_argument('d1', type=int,nargs='?',help="First dimension (count from 0)")
        sbparser.add_argument('d2', type=int,nargs='?',help="Second dimension (count from 0)")
        sbparser.add_argument('inputfiles', type=str,nargs='*',help="Parameter file names, can include more than one.")
    elif arg1=='xy' :
        sbparser = subparsers.add_parser(arg1,formatter_class=argparse.RawTextHelpFormatter,\
            description='Plots one output versus one input.', \
            epilog='Examples: \n"plot_prep.py xy 1 1 ptrain.dat ytrain.dat", \n"plot_prep.py xy 1 1 qval.dat yval.dat"')
        sbparser.add_argument('d', type=int,nargs='?',help="Input dimension (count from 0)")
        sbparser.add_argument('o', type=int,nargs='?',help="Output number (count from 0)")
        sbparser.add_argument('inputfile', type=str,nargs='?',help="Input parameter file name, matrix of size Nx(ndim).")
        sbparser.add_argument('outputfile', type=str,nargs='?',help="Output file name, matrix of size Nx(nout).")
    elif arg1=='xxy' :
        sbparser = subparsers.add_parser(arg1,formatter_class=argparse.RawTextHelpFormatter,\
            description='Plots one output versus two inputs.', \
            epilog='Examples: \n"plot_prep.py xxy 0 1 1 ptrain.dat ytrain.dat", \n"plot_prep.py xxy 0 1 1 qval.dat yval.dat"')
        sbparser.add_argument('d1', type=int,nargs='?',help="First input dimension (count from 0)")
        sbparser.add_argument('d2', type=int,nargs='?',help="Second input dimension (count from 0)")
        sbparser.add_argument('o', type=int,nargs='?',help="Output number (count from 0)")
        sbparser.add_argument('inputfile', type=str,nargs='?',help="Input parameter file name, matrix of size Nx(ndim).")
        sbparser.add_argument('outputfile', type=str,nargs='?',help="Output file name, matrix of size Nx(nout).")
    else:
        print("The code should not get here. Please double check the arguments. Exiting.")
        sys.exit()

args = parser.parse_args()
plot_type=sys.argv[1]


## Parallel coordinates
if(plot_type=='pcoord'):

    # Parse the arguments
    ndim=args.ndim #int(sys.argv[2])


    if os.path.exists('pnames.txt'):
        with open('pnames.txt') as f:
            pnames = f.read().splitlines()
            assert(len(pnames)==ndim)
    else:
        pnames=['Param # '+str(i) for i in range(1,ndim+1)]


    ndcut=args.ndcut #int(sys.argv[3])
    if ndcut==0:
        ndcut=ndim

    ndg=int((ndim-1)/ndcut)+1

    inputfiles=args.inputfiles #sys.argv[4:]


    inputs=np.empty((0,ndim))
    labels=[]
    for inputfile in inputfiles:
        this_input=np.loadtxt(inputfile)
        inputs=np.vstack((inputs,this_input))
        labels.extend([inputfile]*this_input.shape[0])


    #labels= (inputs[:,iout]>0.0)#(out_ss==0.0) #(out<0.004)



    for i in range(ndg):
        print("Plotting %d / %d " % (i+1,ndg))
        names=range(1+i*ndcut,min(1+(i+1)*ndcut,ndim+1))

        values=inputs[:,i*ndcut:min((i+1)*ndcut,ndim)].T
        io.parallel_coordinates(names, values, labels,'pcoord_'+str(i+1)+'.eps')

        #labels_only=labels[labels==True]
        #values_only=values[:,labels==True]
        #io.parallel_coordinates(names, values_only, labels_only)

## Plot one input versus another
elif(plot_type=='xx'):

    ndim=args.ndim #int(sys.argv[2])
    d1=args.d1 #int(sys.argv[3])
    d2=args.d2 #int(sys.argv[4])


    if os.path.exists('pnames.txt'):
        with open('pnames.txt') as f:
            pnames = f.read().splitlines()
            pname1=pnames[d1]
            pname2=pnames[d2]
    else:
        pnames=['Param # '+str(i) for i in range(1,ndim+1)]
        pname1,pname2=['Param # '+str(i) for i in (d1,d2)]

    print("Plotting %s vs %s" % (pname1,pname2))

    inputfiles=args.inputfiles #sys.argv[5:]

    inputs=np.empty((0,ndim))
    labels=[]
    for inputfile in inputfiles:
        this_input=np.loadtxt(inputfile)
        inputs=np.vstack((inputs,this_input))
        labels.extend([inputfile]*this_input.shape[0])

    io.plot_xx(d1,d2,pnames, inputs, labels,'xx_'+pname1+'_'+pname2+'.eps')

## Plot output vs input
elif(plot_type=='xy'):

    d=args.d #int(sys.argv[2])
    o=args.o #int(sys.argv[3])


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



    print("Plotting %s vs %s" % (pname,outname))

    inputfile=args.inputfile #sys.argv[4]
    outputfile=args.outputfile #sys.argv[5]

    input=np.loadtxt(inputfile,ndmin=2)
    output=np.loadtxt(outputfile,ndmin=2)

    io.plot_xy(input[:,d],output[:,o],pname, outname,'xy_'+pname+'_'+outname+'.eps')

# Plot output vs input1,input2
elif(plot_type=='xxy'):

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



    print("Plotting %s,%s vs %s" % (pnames[0],pnames[1],outname))

    inputfile=sys.argv[5]
    outputfile=sys.argv[6]

    input=np.loadtxt(inputfile,ndmin=2)
    output=np.loadtxt(outputfile,ndmin=2)

    io.plot_xxy(input[:,d1],input[:,d2],output[:,o],[pnames[d1],pnames[d2]], outname)


else:
    print("plot_type not recognized. Exiting.")
    sys.exit()
