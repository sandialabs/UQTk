#!/usr/bin/env python
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
"""
Plotting scripts for samples, covariance matrices, and KL modes for
1D and 2D configurations. It takes 2-4 command line arguments depending
on the value of the first argument
"""
import os
import sys
import numpy as np
from   scipy import stats
import matplotlib
import matplotlib.pylab as plt
import matplotlib.tri as tri
#from   pylab import *

from pyutils import readfile, column, checkPtInside

sigma = 5.0
Npl   = 5

rtype="anlcov"
if (len(sys.argv) > 1):
    rtype  = sys.argv[1]
    if rtype == "samples":
        clen  = sys.argv[2]
    if (rtype == "pltKLrecon1D") | (rtype == "pltKLrecon2D"):
        clen  = sys.argv[2]
        ids   = int(sys.argv[3])
        Nkl   = int(sys.argv[4])
        step  = int(sys.argv[5])
    if (rtype == "pltKLeig1D") | (rtype == "pltKLeig2D"):
        nspl  = sys.argv[2]
        if (len(sys.argv) > 3):
            nspl1  = sys.argv[3]
    if rtype == "anlcov":
        ctype = sys.argv[2]
        clen  = sys.argv[3]
    if rtype == "numcov":
        clen  = sys.argv[2]
        nreal = sys.argv[3]
    if rtype == "anlKLevec":
        ctype = sys.argv[2]
        clen  = sys.argv[3]
    if rtype == "numKLevec":
        clen  = sys.argv[2]
        nreal = sys.argv[3]
        dolg  = sys.argv[4]
    if rtype == "xidata1D":
        clen  = sys.argv[2]
        nreal = sys.argv[3]
    if (rtype == "samples2D") | (rtype == "samples2Du"):
        clen  = sys.argv[2]
        nreal = sys.argv[3]
        nspl  = int(sys.argv[4])
    if ( rtype == "anlcov2D" ) | ( rtype == "anlcov2Du" ):
        ctype = sys.argv[2]
        clen  = sys.argv[3]
    if ( rtype == "numcov2D" ) | ( rtype == "numcov2Du" ) :
        clen  = sys.argv[2]
        nreal = sys.argv[3]
    if ( rtype == "anlKLevec2D" ) | ( rtype == "anlKLevec2Du" ):
        ctype = sys.argv[2]
        clen  = sys.argv[3]
    if ( rtype == "numKLevec2D" ) | (rtype == "numKLevec2Du" ) :
        clen  = sys.argv[2]
        nreal = sys.argv[3]

if rtype == "samples":
    fname = "cvspl_"+clen+"_512/samples_"+clen+"_512.dat"
    print("Processing file ",fname)
    din,nliles=readfile(fname);
    # Plot samples
    fs1=18
    lw1=2
    xp=column(din,0);
    fig = plt.figure(figsize=(4,4))
    ax = fig.add_axes([0.19, 0.13, 0.76, 0.82]) 
    for i in range(1,Npl+1):
        plt.plot(xp,column(din,i),linewidth=lw1)
    plt.xlabel('$x$',fontsize=fs1)
    plt.ylabel(r'$F(x,\theta)$',fontsize=fs1)
    ax.set_ylim([-4*sigma,4*sigma])
    ax.set_yticks([-4*sigma,-2*sigma,0,2*sigma,4*sigma])
    plt.savefig("rf1D_"+clen+".pdf")

if rtype == "pltKLeig1D":
    #parameters
    lw1 = 2
    fs1 = 18
    clrs = ['b','g','r']
    lst  = ['--','-']
    # plot eigenvalue spectrum
    fig = plt.figure(figsize=(6,4))
    ax=fig.add_axes([0.15, 0.15, 0.80, 0.80]) 
    pleg=[]
    corl=["0.05","0.10","0.20"];
    i=-1
    for clen in corl:
        i = i+1
        fname = "cvspl_"+clen+"_"+nspl+"/eig_"+str(clen)+"_"+nspl+".dat"
        din,nliles=readfile(fname);
        eig=column(din,0);
        if nspl1 > 0:
            fname = "cvspl_"+clen+"_"+nspl1+"/eig_"+str(clen)+"_"+nspl1+".dat"
            din,nliles=readfile(fname);
            eig1=column(din,0);
            pleg.append(plt.plot(eig, linestyle='--',linewidth=lw1,color=clrs[i]))
            pleg.append(plt.plot(eig1,linestyle='-', linewidth=lw1,color=clrs[i]))
        else:
            pleg.append(plt.plot(eig,linewidth=lw1,color=clrs[i]))
    if nspl1>0:
        leg=plt.legend( (pleg[0][0], pleg[2][0], pleg[4][0]),
                    (r"$c_l="+corl[0]+"$",
                     r"$c_l="+corl[1]+"$",
                     r"$c_l="+corl[2]+"$"),'upper right' )
    else:
        leg=plt.legend( (pleg[0][0], pleg[1][0], pleg[2][0]),
                        (r"$c_l="+corl[0]+"$",
                        r"$c_l="+corl[1]+"$",
                        r"$c_l="+corl[2]+"$"),'upper right' )
    
    plt.xlabel("Eigenvalue #",        fontsize=fs1)
    plt.ylabel("Eigenvalue Magnitude",fontsize=fs1)
    ax.set_yscale("log")
    ax.set_ylim([1.e-6,100])
    ax.set_yticks([1.e-6,1.e-4,1.e-2,1.,1.e2])
    ax.set_xlim([0,50])
    ax.set_xticks([0,10,20,30,40,50])
    plt.savefig("eig_lc_"+nspl+".pdf")
    #plt.show()

if rtype == "pltKLeig2D":
    # parameters
    lw1 = 2
    fs1 = 16
    # plot eigenvalue spectrum
    fig = plt.figure(figsize=(4,4))
    ax=fig.add_axes([0.20, 0.15, 0.75, 0.75]) 
    pleg=[]
    corl=["0.1","0.2","0.5"];
    for clen in corl:
        fname = "cvspl2D_"+clen+"_"+nspl+"/eig2D_"+clen+"_"+nspl+".dat"
        din,nliles=readfile(fname);
        eig=column(din,0);
        pleg.append(plt.plot(eig,linewidth=lw1))
    plt.xlabel("Eigenvalue #",        fontsize=fs1)
    plt.ylabel("Eigenvalue Magnitude",fontsize=fs1)
    ax.set_yscale("log")
    ax.set_ylim([1.e-6,100])
    ax.set_yticks([1.e-6,1.e-4,1.e-2,1,1.e2])
    ax.set_xlim([0,128])
    ax.set_xticks([0,40,80,120])
    leg=plt.legend( (pleg[0][0], pleg[1][0], pleg[2][0]),
                    (r"$c_l="+corl[0]+"$",
                     r"$c_l="+corl[1]+"$",
                     r"$c_l="+corl[2]+"$"),'upper right' )
    plt.savefig("eig_lc2D_"+nspl+".pdf")
    #plt.show()

if rtype == "pltKLrecon1D":
    fspls = "cvspl_"+clen+"_512/samples_"+clen+"_512.dat"
    print("Processing file ",fspls)
    spls,nliles=readfile(fspls);
    feig = "cvspl_"+clen+"_512/eig_"+clen+"_512.dat"
    eig,nliles=readfile(feig)
    fKLm = "cvspl_"+clen+"_512/KLmodes_"+clen+"_512.dat"
    KLmodes,nliles=readfile(fKLm);
    fxi = "cvspl_"+clen+"_512/xidata_"+clen+"_512.dat"
    xi,nlines=readfile(fxi);
    #parameters
    lw1 = 2
    fs1 = 18
    xp=column(spls,0);
    npts=np.array(xp).size
    savg=[np.mean(spls[k][1::]) for k in range(npts)]
    fn=[]
    for n in range(Nkl):
        fn.append([0.0]*npts)
        for i in range(npts):
            if ( n == 0 ):
                fn[n][i] = savg[i]
            else:
                fn[n][i] = fn[n-1][i]+KLmodes[i][n]*xi[ids][n-1] 
    for n in range(0,Nkl,step):
        fig = plt.figure(figsize=(6,4))
        ax=fig.add_axes([0.15, 0.15, 0.80, 0.80]) 
        plt.plot(xp,column(spls,ids+1),color='black',linewidth=lw1)
        plt.plot(xp,fn[n],color='red',linewidth=lw1)
        plt.xlabel(r"$x$",     fontsize=fs1)
        plt.ylabel(r"$F_n(x)$",fontsize=fs1)
        ax.set_ylim([-15,15])
        ax.set_yticks([-10,0,10])
        ax.set_xlim([0,1])
        #ax.set_xticks([0,10,20,30,40])
        plt.text(0.72,17,str(n)+" terms",fontsize=18)
        plt.savefig("KLrecon_"+clen+"_"+str(n)+".pdf")

if rtype == "pltKLrecon2D":
    lw1 = 3
    fs1 = 22
    xp,nx=readfile("cvspl2D_"+clen+"_4096/xgrid.dat");
    yp,ny=readfile("cvspl2D_"+clen+"_4096/ygrid.dat");
    x=column(xp,0);
    y=column(yp,0);
    X,Y=np.meshgrid(x,y)
    fspls = "cvspl2D_"+clen+"_4096/samples2D_"+clen+"_4096.dat"
    samples,nlines=readfile(fspls);
    feig = "cvspl2D_"+clen+"_4096/eig2D_"+clen+"_4096.dat"
    eig,nlines=readfile(feig)
    fKLm = "cvspl2D_"+clen+"_4096/KLmodes2D_"+clen+"_4096.dat"
    KLmodes,nlines=readfile(fKLm);
    fxi = "cvspl2D_"+clen+"_4096/xidata2D_"+clen+"_4096.dat"
    xi,nlines=readfile(fxi);
    savg=[np.mean(samples[k][:]) for k in range(nx*ny)]
    fn=[]
    for n in range(Nkl):
        fn.append([0.0]*nx*ny)
        for i in range(nx*ny):
            if ( n == 0 ):
                fn[n][i] = savg[i]
            else:
                fn[n][i] = fn[n-1][i]+KLmodes[i][n+1]*xi[ids][n-1]
    for n in range(0,Nkl,step):
        Zex=np.array(column(samples,ids)).reshape(nx,ny)
        Z=np.array(fn[n]).reshape(nx,ny)
        fig = plt.figure(figsize=(4,4))
        ax=fig.add_axes([0.15, 0.15, 0.8, 0.8]) 
        plt.contour(X,Y,Zex,11,colors="black",linewidth=lw1)
        plt.contour(X,Y,Z,  11,colors="red",  linewidth=lw1)
        plt.xlabel('$x_1$',fontsize=fs1)
        plt.ylabel('$x_2$',fontsize=fs1)
        ax.set_xlim([0.0,1.0])
        ax.set_ylim([0.0,1.0])
        ax.set_aspect('equal')
        plt.savefig("KLrecon2D_"+clen+"_"+str(n)+".pdf")

if rtype == "anlcov":
    fname = "klcov_"+ctype+"_"+clen+"/cov_"+clen+"_"+ctype+"_anl.dat"
    print("Processing file ",fname)
    cov,nlines=readfile(fname);
    vmax = np.array(cov).max()
    fig=plt.figure(figsize=(4, 4))
    ax = fig.add_axes([0.05, 0.05, 0.85, 0.85]) 
    plt.imshow(cov, interpolation='nearest', vmin=-0.5*vmax, vmax=vmax,
              cmap=plt.cm.RdBu_r)
    plt.xticks(())
    plt.yticks(())
    plt.title("$c_l=$"+clen)
    plt.savefig("cov_"+ctype+"_"+clen+"_anl.eps")

if rtype == "numcov":
    fname = "cvspl_"+clen+"_"+nreal+"/cov_"+clen+"_"+nreal+".dat"
    print("Processing file ",fname)
    cov,nlines=readfile(fname);
    vmax = np.array(cov).max()
    fig=plt.figure(figsize=(4, 4))
    ax = fig.add_axes([0.02, 0.02, 0.96, 0.96]) 
    plt.imshow(cov, interpolation='nearest', vmin=-0.5*vmax, vmax=vmax,
              cmap=plt.cm.RdBu_r)
    plt.xticks(())
    plt.yticks(())
    #plt.title(r"$c_l="+clen+", N_{\Theta}="+nreal+"$")
    plt.savefig("cov_"+clen+"_"+nreal+"_num.pdf")

if rtype == "anlKLevec":
    # plot KL modes
    fname = "klcov_"+ctype+"_"+clen+"/KLmodes_"+clen+"_"+ctype+"_anl.dat"
    din,nliles=readfile(fname);
    #parameters
    lw=2
    fs=16
    Nkl=4
    xp=column(din,0);
    fig = plt.figure(figsize=(4,6))
    ax=fig.add_axes([0.17, 0.15, 0.75, 0.75]) 
    pleg=[]
    for i in range(1,Nkl+1):
        y=column(din,i);
        if ( clen == "0.10" ):
            if ( i == 1 ):
                y=(-1)*np.array(y)
        if ( clen == "0.20" ):
            if ( i == 2 ) | (i == 3) :
                y=(-1)*np.array(y)
        pleg.append(plt.plot(xp,y,linewidth=lw))
    plt.xlabel("x",fontsize=fs)
    plt.ylabel(r"$\sqrt{\lambda_k}f_k$",fontsize=fs)
    ax.set_ylim([-4,4])
    #ax.set_yticks([1.e-4,1.e-2,1.,1.e2])
    #ax.set_xlim([0,40])
    #ax.set_xticks([0,10,20,30,40])
    leg=plt.legend( (pleg[0][0], pleg[1][0], pleg[2][0], pleg[3][0]),
                    (r"$f_1$", r"$f_2$", r"$f_3$", r"$f_4$"),'lower right' )
    plt.title("$c_l=$"+clen)
    plt.savefig("KLmodes_"+ctype+"_"+clen+"_anl.eps")

if rtype == "numKLevec":
    # plot KL modes
    fname = "cvspl_"+clen+"_"+nreal+"/KLmodes_"+clen+"_"+nreal+".dat"
    din,nliles=readfile(fname);
    #parameters
    lw1 = 3
    fs1 = 22
    Nkl = 4
    xp=column(din,0);
    fig = plt.figure(figsize=(4,6))
    ax=fig.add_axes([0.18, 0.10, 0.77, 0.87]) 
    pleg=[]
    for i in range(1,Nkl+1):
        y=column(din,i);
        if ( clen == "0.20" ):
            if ( nreal == "8192" ):
                if ( i == 1 ):
                    y=(-1)*np.array(y)
                if ( i == 4 ):
                    y=(-1)*np.array(y)
            if ( nreal == "512" ):
                if ( i == 1 ):
                    y=(-1)*np.array(y)
        if ( clen == "0.05" ):
            if ( nreal == "512" ):
                if ( i == 1 ):
                    y = (-1)*np.array(y)
                if ( i == 2 ) | ( i == 3 ):
                    y = (-1)*np.array(y)
            if ( nreal == "8192" ):
                if ( i == 1 ):
                    y = (-1)*np.array(y)
        if ( clen == "0.10" ):
            if ( nreal == "512" ):
                if ( i == 1 ):
                    y=(-1)*np.array(y)
                if ( i == 3 ):
                    y=(-1)*np.array(y)
        pleg.append(plt.plot(xp,y,linewidth=lw1))
    plt.xlabel("$x$",fontsize=fs1)
    plt.ylabel(r"$\sqrt{\lambda_k}f_k$",fontsize=fs1,labelpad=-5)
    ax.set_ylim([-4,4])
    ax.set_yticks([-4,-2,0,2,4])
    #ax.set_xlim([0,40])
    #ax.set_xticks([0,10,20,30,40])
    if dolg == "on":
        leg=plt.legend( (pleg[0][0], pleg[1][0], pleg[2][0], pleg[3][0]),
                        (r"$f_1$", r"$f_2$", r"$f_3$", r"$f_4$"),'lower right' )
    #plt.title(r"$c_l="+clen+", N_{\Theta}="+nreal+"$")
    plt.savefig("KLmodes_"+clen+"_"+nreal+".pdf")

if rtype == "xidata1D":
    fname = "cvspl_"+clen+"_"+nreal+"/xidata_"+clen+"_"+nreal+".dat"
    din,nliles=readfile(fname);
    #parameters
    lw1 = 3
    fs1 = 18
    fig = plt.figure(figsize=(6,4))
    ax=fig.add_axes([0.12, 0.15, 0.85, 0.82]) 
    pleg = []
    Nxi  = 4
    for i in range(Nxi):
        x=np.linspace(np.array(column(din,i)).min(),np.array(column(din,i)).max(),100)
        kerns=stats.kde.gaussian_kde(column(din,i))
        pleg.append(plt.plot(x,kerns(x),linewidth=lw1))
    #pleg.append(plt.plot(x,np.exp(-x**2/(2*sigma*sigma))/(np.sqrt(2*np.pi)*sigma),linewidth=lw))
    plt.xlabel(r"$\xi(\theta)$",fontsize=fs1)
    plt.ylabel(r"$PDF(\xi)$",   fontsize=fs1)
    ax.set_xlim([-4,4])
    ax.set_ylim([0.0,0.45])
    ax.set_yticks([0.0, 0.10, 0.20, 0.30, 0.40])
    ax.set_xticks([-4,-2,0,2,4])
    leg=plt.legend( (pleg[0][0], pleg[1][0], pleg[2][0], pleg[3][0]),
                    (r"$\xi_1$", r"$\xi_2$", r"$\xi_3$", r"$\xi_4$"),'upper right' )
    #plt.title(r"$c_l="+clen+", N_{\Theta}="+nreal+"$")
    plt.savefig("xidata1D_"+clen+"_"+nreal+".pdf")


if rtype == "samples2D":
    #parameters
    fs1   = 22
    fs2   = 16
    fsize = 4
    #fs1   = 11
    #fs2   = 6
    #fsize = 2
    # plot KL modes
    fname = "cvspl2D_"+clen+"_"+nreal+"/samples2D_"+clen+"_"+nreal+".dat"
    din,nlines=readfile(fname);
    x,nx=readfile("cvspl2D_"+clen+"_"+nreal+"/xgrid.dat");
    y,ny=readfile("cvspl2D_"+clen+"_"+nreal+"/ygrid.dat");
    x=column(x,0);
    y=column(y,0);
    X,Y=np.meshgrid(x,y)
    for i in range(nspl):
        z=column(din,i);
        Z=[[z[k*ny+j] for j in range(ny)] for k in range(nx)]
        fig = plt.figure(figsize=(fsize,fsize))
        ax=fig.add_axes([0.15, 0.15, 0.8, 0.8]) 
        plt.contour(X,Y,Z,21,rasterized=True)
        plt.xlabel("$x_1$",fontsize=fs1)
        plt.ylabel("$x_2$",fontsize=fs1,labelpad=-2)
        ax.set_xlim([0.0,1.0])
        ax.set_ylim([0.0,1.0])
        for tick in ax.xaxis.get_major_ticks()+ax.yaxis.get_major_ticks():
            tick.label.set_fontsize(fs2)
        #ax.set_aspect('equal')
        plt.savefig("samples2D_"+clen+"_"+nreal+"_s"+str(i+1)+".pdf")

if rtype == "anlcov2D":
    fname = "klcov2D_"+ctype+"_"+clen+"/cov2D_"+clen+"_"+ctype+"_anl.dat"
    print("Processing file ",fname)
    cov,nlines=readfile(fname);
    vmax = np.array(cov).max()
    fig=plt.figure(figsize=(4, 4))
    ax = fig.add_axes([0.05, 0.05, 0.85, 0.85]) 
    plt.imshow(cov, interpolation='nearest', vmin=-0.5*vmax, vmax=vmax,
              cmap=plt.cm.RdBu_r)
    plt.xticks(())
    plt.yticks(())
    plt.title("$c_l=$"+clen)
    plt.savefig("cov2D_"+ctype+"_"+clen+"_anl.eps")

if rtype == "numcov2D":
    fname = "klsampl2D_"+clen+"_"+nreal+"/cov2D_"+clen+"_"+nreal+".dat"
    print("Processing file ",fname)
    cov,nlines=readfile(fname);
    vmax = np.array(cov).max()
    fig=plt.figure(figsize=(4, 4))
    ax = fig.add_axes([0.05, 0.05, 0.85, 0.85]) 
    plt.imshow(cov, interpolation='nearest', vmin=-0.5*vmax, vmax=vmax,
              cmap=plt.cm.RdBu_r)
    plt.xticks(())
    plt.yticks(())
    plt.title(r"$c_l="+clen+", N_{\Theta}="+nreal+"$")
    plt.savefig("cov2D_"+clen+"_"+nreal+"_num.eps")

if rtype == "anlKLevec2D":
    #parameters
    lw=2
    fs=16
    # plot KL modes
    Nkl=6
    fname = "klcov2D_"+ctype+"_"+clen+"/KLmodes2D_"+clen+"_"+ctype+"_anl.dat"
    din,nlines=readfile(fname);
    x=[din[i][0] for i in range(65)];
    y=[din[65*i][1] for i in range(65)];
    X,Y=np.meshgrid(x,y)
    for i in range(2,Nkl+2):
        z=column(din,i);
        Z=[[z[k*65+j] for j in range(65)] for k in range(65)]
        fig = plt.figure(figsize=(4,4))
        ax=fig.add_axes([0.15, 0.15, 0.75, 0.75]) 
        plt.contourf(X,Y,Z,101)
        plt.xlabel("x",fontsize=fs)
        plt.ylabel("y",fontsize=fs)
        ax.set_xlim([0.0,1.0])
        ax.set_ylim([0.0,1.0])
        ax.set_aspect('equal')
        plt.title("$c_l=$"+clen+", Mode "+str(i-1))
        plt.savefig("KLmodes2D_"+ctype+"_"+clen+"_m"+str(i-1)+".eps")

if rtype == "numKLevec2D":
    #parameters
    lw1 = 2
    fs1 = 22
    # plot KL modes
    Nkl = 8
    fname = "cvspl2D_"+clen+"_"+nreal+"/KLmodes2D_"+clen+"_"+nreal+".dat"
    din,nlines=readfile(fname)
    nx = ny = int(np.sqrt(nlines))
    x=[din[i][0] for i in range(nx)];
    y=[din[nx*i][1] for i in range(ny)]
    X,Y=np.meshgrid(x,y)
    for i in range(2,Nkl+2):
        z=column(din,i);
        Z=[[z[k*nx+j] for j in range(ny)] for k in range(nx)]
        fig = plt.figure(figsize=(4,4))
        ax=fig.add_axes([0.15, 0.15, 0.8, 0.8]) 
        plt.contourf(X,Y,Z,21)
        plt.xlabel("$x_1$",fontsize=fs1)
        plt.ylabel("$x_2$",fontsize=fs1)
        ax.set_xlim([0.0,1.0])
        ax.set_ylim([0.0,1.0])
        ax.set_aspect('equal')
        #plt.title("$c_l=$"+clen+", Mode "+str(i-1))
        plt.savefig("KLmodes2D_"+clen+"_"+nreal+"_m"+str(i-1)+".pdf")

if rtype == "samples2Du":
    #parameters
    lw1 = 2
    fs1 = 18
    # plot KL modes
    fname = "cvspl2Du_"+clen+"_"+nreal+"/samples2Du_"+clen+"_"+nreal+".dat"
    din,nlines=readfile(fname);
    xy,nxy=readfile("data/cali_grid.dat");
    x = np.array(column(xy,0))
    y = np.array(column(xy,1))
    # Create the Triangulation; no triangles so Delaunay triangulation created.
    triGrid = tri.Triangulation(x,y)
    # Mask off unwanted triangles.
    xyc,nl=readfile('data/cali.dat')
    xc = np.array(column(xyc,0))
    yc = np.array(column(xyc,1))
    xmid = x[triGrid.triangles].mean(axis=1)
    ymid = y[triGrid.triangles].mean(axis=1)
    mask = checkPtInside(xc,yc,xmid,ymid);
    #mask = np.where(xmid*xmid + ymid*ymid < min_radius*min_radius, 1, 0)
    triGrid.set_mask(mask)
    for i in range(nspl):
        fig = plt.figure(figsize=(4,4))
        ax=fig.add_axes([0.08, 0.08, 0.9, 0.9]) 
        ax.set_aspect('equal')
        plt.tricontour(triGrid,column(din,i),21)
        ax.set_xlabel("lon",fontsize=fs1)
        ax.set_ylabel("lat",fontsize=fs1)
        ax.set_xlim([x.min(),x.max()])
        ax.set_ylim([y.min(),y.max()])
        ax.set_xticks([])
        ax.set_yticks([])
        #plt.colorbar()
        #plt.title("$c_l=$"+clen)
        plt.savefig("samples2Du_"+clen+"_"+nreal+"_s"+str(i+1)+".pdf")

if rtype == "anlcov2Du":
    fname = "cvanl2Du_"+ctype+"_"+clen+"/cov2Du_"+ctype+"_"+clen+".dat"
    print("Processing file ",fname)
    cov,nlines=readfile(fname);
    vmax = np.array(cov).max()
    fig=plt.figure(figsize=(4, 4))
    ax = fig.add_axes([0.05, 0.05, 0.85, 0.85]) 
    plt.imshow(cov, interpolation='nearest', vmin=-vmax, vmax=vmax,
              cmap=plt.cm.RdBu_r)
    plt.xticks(())
    plt.yticks(())
    plt.title("$c_l=$"+clen)
    plt.savefig("cov2Du_"+ctype+"_"+clen+"_anl.pdf")

if rtype == "numcov2Du":
    fname = "cvspl2Du_"+clen+"_"+nreal+"/cov2Du_"+clen+"_"+nreal+".dat"
    print("Processing file ",fname)
    cov,nlines=readfile(fname);
    vmax = np.array(cov).max()
    fig=plt.figure(figsize=(4, 4))
    ax = fig.add_axes([0.05, 0.05, 0.85, 0.85]) 
    plt.imshow(cov, interpolation='nearest', vmin=-vmax, vmax=vmax,
              cmap=plt.cm.RdBu_r)
    plt.xticks(())
    plt.yticks(())
    #plt.title(r"$c_l="+clen+", N_{\Theta}="+nreal+"$")
    plt.savefig("cov2Du_"+clen+"_"+nreal+"_num.pdf")

if rtype == "anlKLevec2Du":
    #parameters
    fs1 = 18
    # plot KL modes
    Nkl = 16
    fname = "cvanl2Du_"+ctype+"_"+clen+"/KLmodes2Du_"+ctype+"_"+clen+".dat"
    din,nlines=readfile(fname);
    x = np.array(column(din,0))
    y = np.array(column(din,1))
    # Create the Triangulation; no triangles so Delaunay triangulation created.
    triGrid = tri.Triangulation(x,y)
    # Mask off unwanted triangles.
    xyc,nl=readfile('data/cali.dat')
    xc = np.array(column(xyc,0))
    yc = np.array(column(xyc,1))
    xmid = x[triGrid.triangles].mean(axis=1)
    ymid = y[triGrid.triangles].mean(axis=1)
    mask = checkPtInside(xc,yc,xmid,ymid);
    #mask = np.where(xmid*xmid + ymid*ymid < min_radius*min_radius, 1, 0)
    triGrid.set_mask(mask)
    for i in range(2,Nkl+2):
        fig = plt.figure(figsize=(4,4))
        ax=fig.add_axes([0.08, 0.08, 0.9, 0.9]) 
        ax.set_aspect('equal')
        plt.tricontourf(triGrid,column(din,i),21)
        ax.set_xlabel("lon",fontsize=fs1)
        ax.set_ylabel("lat",fontsize=fs1)
        ax.set_xlim([x.min(),x.max()])
        ax.set_ylim([y.min(),y.max()])
        ax.set_xticks([])
        ax.set_yticks([])
        #plt.colorbar()
        #plt.title("$c_l=$"+clen+", Mode "+str(i-1))
        plt.savefig("KLmodes2Du_"+ctype+"_"+clen+"_m"+str(i-1)+".pdf")


if rtype == "numKLevec2Du":
    # Plot parameters
    fs1 = 18
    # no. of KL modes
    Nkl = 16
    fname = "cvspl2Du_"+clen+"_"+nreal+"/KLmodes2Du_"+clen+"_"+nreal+".dat"
    din,nlines=readfile(fname);
    x = np.array(column(din,0))
    y = np.array(column(din,1))
    # Create the Triangulation; no triangles so Delaunay triangulation created.
    triGrid = tri.Triangulation(x,y)
    # Mask off unwanted triangles.
    xyc,nl=readfile('data/cali.dat')
    xc = np.array(column(xyc,0))
    yc = np.array(column(xyc,1))
    xmid = x[triGrid.triangles].mean(axis=1)
    ymid = y[triGrid.triangles].mean(axis=1)
    mask = checkPtInside(xc,yc,xmid,ymid);
    #mask = np.where(xmid*xmid + ymid*ymid < min_radius*min_radius, 1, 0)
    triGrid.set_mask(mask)
    for i in range(2,Nkl+2):
        fig = plt.figure(figsize=(4,4))
        ax=fig.add_axes([0.08, 0.08, 0.9, 0.9]) 
        ax.set_aspect('equal')
        plt.tricontourf(triGrid,column(din,i),21)
        ax.set_xlabel("lon",fontsize=fs1)
        ax.set_ylabel("lat",fontsize=fs1)
        ax.set_xlim([x.min(),x.max()])
        ax.set_ylim([y.min(),y.max()])
        ax.set_xticks([])
        ax.set_yticks([])
        #plt.colorbar()
        #plt.title("$c_l=$"+clen+", Mode "+str(i-1))
        plt.savefig("KLmodes2Du_"+clen+"_"+nreal+"_m"+str(i-1)+".pdf")





