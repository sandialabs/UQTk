#!/usr/bin/env python
#=====================================================================================
#
#                      The UQ Toolkit (UQTk) version 3.1.5
#                          Copyright (2024) NTESS
#                        https://www.sandia.gov/UQToolkit/
#                        https://github.com/sandialabs/UQTk
#
#     Copyright 2024 National Technology & Engineering Solutions of Sandia, LLC (NTESS).
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
#     Questions? Contact the UQTk Developers at https://github.com/sandialabs/UQTk/discussions
#     Sandia National Laboratories, Livermore, CA, USA
#=====================================================================================
import argparse
import matplotlib.pyplot as plt
import matplotlib.tri as tri
from matplotlib.lines import Line2D                        
import numpy as np
import math

THRSMALL = 1.e-6

def trArea(t0, t1, t2):
    '''
    Area of a 2D triangle
    '''
    return abs(t0[0]*(t1[1]-t2[1])+t1[0]*(t2[1]-t0[1])+t2[0]*(t0[1]-t1[1]))/2.0 

def checkPtInside(borders,x,y):
    nm=len(x);
    umask=np.zeros(nm)
    for i in range(nm):
        b=y[i]-x[i];
        nint = 0;
        for brd in borders:
            xc=brd[:,0];
            yc=brd[:,1];
            for j in range(xc.shape[0]-1):
                if abs(xc[j]-xc[j+1]) < 1.0e-30:
                    yi=xc[j]+b;
                    if (xc[j]<=x[i]) & ((yi-yc[j])*(yi-yc[j+1])<=0.0):
                        nint = nint+1;
                else:
                    a1=(yc[j+1]-yc[j]) / (xc[j+1]-xc[j]);
                    b1=yc[j]-a1*xc[j];
                    xi=(b1-b)/(1.0-a1);
                    if (xi<=x[i]) & ((xi-xc[j])*(xi-xc[j+1])<=0.0):
                        nint = nint+1;
        if ( nint%2 == 0 ):
            umask[i]=1;
    return umask;


def check_Point_In_Region(borders, x, y):
    b = y - x
    for region in borders:
        nint = 0
        xc = region[:,0]
        yc = region[:,1]
        for j in range(xc.shape[0]-1):
            if abs(xc[j]-xc[j+1]) < 1.0e-30:
                yi = xc[j] + b
                if (xc[j] <= x) & ((yi-yc[j]) * (yi-yc[j+1]) <= 0.0):
                    nint = nint+1
            else:
                a1 = (yc[j+1]-yc[j]) / (xc[j+1]-xc[j])
                b1 = yc[j]-a1*xc[j]
                xi = (b1-b) / (1.0-a1)
                if (xi<=x) & ((xi-xc[j])*(xi-xc[j+1])<=0.0):
                    nint = nint+1
        if ( nint%2 == 1 ):
            return True
    return False

def gen_samples_region(regions, npoints, seed = 2020):
    '''
    Generate samples inside a set of regions
    '''

    xmin,xmax = min([bnd[:,0].min() for bnd in regions]), max([bnd[:,0].max() for bnd in regions])
    ymin,ymax = min([bnd[:,1].min() for bnd in regions]), max([bnd[:,1].max() for bnd in regions])

    np.random.seed(seed)
    i = 0
    xy = np.zeros((npoints,2))
    while i<npoints:
        x = np.random.uniform(xmin,xmax)
        y = np.random.uniform(ymin,ymax)
        if check_Point_In_Region(regions, x, y):
            xy[i] = np.array([x,y])
            i=i+1

    for region in regions:
        xy = np.vstack((xy,region))

    return xy


def main(): 

    # Parse aeguments
    parser = argparse.ArgumentParser(prog='kl_prep',
                    description='Generate triangulation over a generic area')
    parser.add_argument('-r','--regions', dest='regions',  nargs='+', type=str, help='''list of regions''')
    parser.add_argument('-n','--npoints', dest='npoints', type=int, default = 256, help='''num. points''')
    parser.add_argument('-o','--output', dest='output', type=str, default = '', help='''output''')
    parser.add_argument('-s','--seed', dest='seed', type=int, default = 2020, help='''num. points''')
    args = parser.parse_args()

    print('Regions:', args.regions)
    print('No. of samples:', args.npoints)
    print('Output arguments:', args.output)

    # read boundaries
    borders = []
    for i in args.regions:
        borders.append(np.genfromtxt(i+'.dat'))

    # generate samples
    xy = gen_samples_region(borders, args.npoints, seed = args.seed)

    # Create triangulation
    triang = tri.Triangulation(xy[:,0],xy[:,1])
    xmid = xy[triang.triangles,0].mean(axis=1)
    ymid = xy[triang.triangles,1].mean(axis=1)

    # Mask off unwanted triangles
    mask = np.zeros(triang.triangles.shape[0])
    for ii, t in enumerate(triang.triangles):
        if not check_Point_In_Region(borders, xmid[ii], ymid[ii]):
            mask[ii] = 1

    triang.set_mask(mask)
    triangles_valid = triang.get_masked_triangles()

    np.savetxt('cali_grid'+args.output+'.dat', xy, fmt='%.10e %.10e', delimiter=' ', newline='\n')
    np.savetxt('cali_tria'+args.output+'.dat', triangles_valid, fmt='%d %d %d', delimiter=' ', newline='\n')

    # plot samples and triangles
    fig1 = plt.figure(figsize=[8,7])                                
    fig2 = plt.figure(figsize=[8,7])                                
    ax1 = fig1.add_axes([0.12,0.12,0.85,0.85])
    ax2 = fig2.add_axes([0.12,0.12,0.85,0.85])

    ax1.scatter(xy[:,0], xy[:,1], s=6)
    for k in range(len(triangles_valid)):
        x1, x2, x3 = xy[triangles_valid[k,:], 0]
        y1, y2, y3 = xy[triangles_valid[k,:], 1]
        aaa = ax2.add_line(Line2D([x1,x2],[y1,y2]))
        aaa = ax2.add_line(Line2D([x2,x3],[y2,y3]))
        aaa = ax2.add_line(Line2D([x3,x1],[y3,y1]))

    for brd in borders:
        ax1.plot(brd[:,0], brd[:,1], 'r', lw=3)
        ax2.plot(brd[:,0], brd[:,1], 'r', lw=3)

    for ax in [ax1, ax2]:
        ax.set_aspect('equal')
        ax.set_xlabel("longitude", fontsize=18)
        ax.set_ylabel("latitude",fontsize=18)
        for tick in ax.xaxis.get_major_ticks() + ax.yaxis.get_major_ticks():
            tick.label.set_fontsize(14) 

    fig1.savefig('cali_pts'+args.output+'.pdf')
    fig2.savefig('cali_tri'+args.output+'.pdf')
    plt.close(fig1)
    plt.close(fig2)
    #plt.show()


if __name__ == "__main__":
    main()

