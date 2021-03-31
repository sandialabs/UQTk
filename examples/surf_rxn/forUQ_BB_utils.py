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

import numpy as npy
import scipy as spy
from scipy import stats, mgrid, c_, reshape, random, rot90
import matplotlib.pyplot as plt
import os
import sys
import shutil
import fileinput

def funcBB(inpfile,outfile,**kwargs):
    
    # get var arguments
    xmltpl = kwargs.get("xmltpl", "surf_rxn.in.xml.tp3")   # input xml template
    xmlin  = kwargs.get("xmlin",  "surf_rxn.in.xml")         # input xml file
    spid   = kwargs.get("spid",  1)                          # species id to analyze
    fstore = kwargs.get("fstore", False)

    # Load input samples
    indata=npy.loadtxt(inpfile)
    if len(indata.shape)==1:
        indata=npy.array(npy.transpose([indata]))

    nruns,ndim=indata.shape    

    open("logSurfRxn.dat", 'w').close()
    ave=[]
    for ip in range(nruns):

        # Create the xml files
        shutil.copyfile(xmltpl,xmlin)
        for idim in range(1,ndim+1):
            for line in fileinput.input("surf_rxn.in.xml", inplace = 1):
                print(line.replace('PAR_'+str(idim), str(indata[ip,idim-1])), end=' ')

        #print "Running the ODE solver for parameter %d/%d"%(ip+1,nruns),
        os.system("./SurfRxnDet.x >> logSurfRxn.dat")

        # Store the solution file
        if fstore:
            shutil.copyfile('solution.dat','solution_'+str(ip)+'.dat')

        # Load solution and get average of the second half of the time series 
        cursol=npy.loadtxt("solution.dat")
        tail=cursol.shape[0]//2
        ave.append(npy.average(cursol[-tail:,spid]))

        #print "Average: ", npy.array(ave)
    npy.savetxt(outfile, ave, fmt="%.18e", delimiter=" ", newline="\n")

def plotKdeBin(splfile,figout,**kwargs):
    # get var arguments
    npdf  = kwargs.get("npdf", 100)   # no. of points for pdf evaluation
    nbins = kwargs.get("nbins", 50)   # no. of bins

    spls=npy.genfromtxt(splfile)
    # construct bins
    bh,be=npy.histogram(spls,nbins,density=True)
    bc=[]
    for i in range(nbins):
        bc.append(0.5*(be[i]+be[i+1]))
    xlin=npy.linspace(spls.min(),spls.max(),npdf) ;
    kernlin=stats.kde.gaussian_kde(spls);
    pdflin1=kernlin(xlin);
    # check for scipy version; if 11 or above then play with bandwidths
    spv=[]
    spver=spy.__version__
    for i,c in enumerate(spver):
        if c==".": spv.append(i)
    spver=int(spver[spv[0]+1:spv[1]])
    print("scipy version: ",spver)
    if spver > 10:
        kernlin.set_bandwidth(bw_method=kernlin.factor/2.0)
        pdflin2=kernlin(xlin);
        kernlin.set_bandwidth(bw_method=kernlin.factor*4.0)
        pdflin3=kernlin(xlin);
        fig = plt.figure(figsize=(8,6))
        ax=fig.add_axes([0.12,0.10,0.83,0.85])
        pl1,=plt.plot(xlin,pdflin1,linewidth=2,label="optimal")
        pl2,=plt.plot(xlin,pdflin2,linewidth=2,label="(1/2 * optimal)")
        pl3,=plt.plot(xlin,pdflin3,linewidth=2,label="(2 * optimal)")
    else:
        fig = plt.figure(figsize=(8,6))
        ax=fig.add_axes([0.12,0.10,0.83,0.85])
        pl1,=plt.plot(xlin,pdflin1,linewidth=2,label="optimal")
    plb,=plt.plot(bc,bh,linewidth=2,label="binning")
    ax.set_xlabel("PC surrogate",fontsize=24)
    ax.set_ylabel("PDF",fontsize=24)
    plt.legend(loc="upper center")
    plt.savefig(figout)

def modelMC(mu,stdfrac,nspl,outfile,**kwargs):

    # get var arguments
    xmltpl = kwargs.get("xmltpl", "surf_rxn.in.xml.tp3")   # input xml template
    xmlin  = kwargs.get("xmlin",  "surf_rxn.in.xml")       # input xml file
    spid   = kwargs.get("spid",  1)                        # species id to analyze

    # other parameters
    ndim = 1
    
    # Generate samples from a normal distribution with mu and sigma
    sigma=stdfrac*mu
    spls=npy.random.normal(mu, sigma, nspl)

    open("logSurfRxn.dat", 'w').close()
    ave=[]
    for ip in range(nspl):

        if ((ip+1)%(nspl/100) == 0):
            perc=(ip+1)/(nspl/100)
            print("\rSampling model: %d percent done"%(perc), end=' ')
            sys.stdout.flush()
            npy.savetxt(outfile, ave, fmt="%.18e", delimiter=" ", newline="\n")
        
        # Create the xml files
        shutil.copyfile(xmltpl,xmlin)
        for line in fileinput.input("surf_rxn.in.xml", inplace = 1):
            print(line.replace('PAR_1', str(spls[ip])), end=' ')

        #print "Running the ODE solver for parameter %d/%d"%(ip+1,nruns),
        os.system("./SurfRxnDet.x >> logSurfRxn.dat")

        # Load solution and get average of the second half of the time series 
        cursol=npy.loadtxt("solution.dat")
        tail=cursol.shape[0]/2
        ave.append(npy.average(cursol[-tail:,spid]))

        #print "Average: ", npy.array(ave)
        
    npy.savetxt(outfile, ave, fmt="%.18e", delimiter=" ", newline="\n")

def plotFunQdpts(ifun,ofun,iqdpts,oqdpts,figout,**kwargs):

    # get var arguments
    npdf  = kwargs.get("npdf", 100)   # no. of points for pdf evaluation
    nbins = kwargs.get("nbins", 50)   # no. of bins

    xfun=npy.genfromtxt(ifun)
    yfun=npy.genfromtxt(ofun)
    xqdpts=npy.genfromtxt(iqdpts)
    yqdpts=npy.genfromtxt(oqdpts)

    fig = plt.figure(figsize=(8,6))
    ax=fig.add_axes([0.12,0.10,0.83,0.85])
    pl1,=plt.plot(xfun,yfun,linewidth=2,label="model")
    pl2,=plt.plot(xqdpts,yqdpts,linestyle='none',marker="o",markeredgewidth=2,markersize=8,mfc='none',
                  label="quadrature points")
    ax.set_xlabel(r"$b$",fontsize=24)
    ax.set_ylabel(r"$u_{ss}(b)$",fontsize=24)
    plt.legend(loc="upper right")
    plt.savefig(figout)

