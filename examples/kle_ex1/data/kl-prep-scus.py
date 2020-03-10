#!/usr/bin/env python
#=====================================================================================
#
#                      The UQ Toolkit (UQTk) version 3.1.0
#                          Copyright (2020) NTESS
#                        https://www.sandia.gov/UQToolkit/
#                        https://github.com/sandialabs/UQTk
#
#     Copyright 2020 National Technology & Engineering Solutions of Sandia, LLC (NTESS).
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

import matplotlib.pyplot as plt
import matplotlib.tri as tri
from matplotlib.lines import Line2D                        
import numpy as npy
import math

def checkPtInside(borders,x,y):
    nm=len(x);
    umask=npy.zeros(nm)
    for i in range(nm):
        b=y[i]-x[i];
        nint = 0;
        for brd in borders:
            xc=brd[:,0];
            yc=brd[:,1];
            for j in range(xc.shape[0]-1):
                if abs(xc[j]-xc[j+1])<1.0e-30:
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


def trArea(x,y):
    return (abs(x[0]*(y[1]-y[2])+x[1]*(y[2]-y[0])+x[2]*(y[0]-y[1]))/2.0) 

states = ['tx','ok','ks']
cvtfile = 'scus_cvt_4096.dat'

borders = []
for i in states:
    borders.append(npy.genfromtxt(i+'.dat'))

xy=npy.genfromtxt(cvtfile)
x = xy[:,0]
y = xy[:,1]

#nb=len(xc)
#nv=len(x)

# Create the Triangulation; no triangles so Delaunay triangulation created.
triang = tri.Triangulation(x, y)

# Mask off unwanted triangles.
xmid = x[triang.triangles].mean(axis=1)
ymid = y[triang.triangles].mean(axis=1)
mask = checkPtInside(borders,xmid,ymid);
triang.set_mask(mask)

# get valid traingles
forOut=triang.get_masked_triangles()
nt=len(forOut)

npy.savetxt("scus_grid.dat", xy, fmt='%.10e %.10e', delimiter=' ', newline='\n')
npy.savetxt("scus_tria.dat", forOut, fmt='%d %d %d', delimiter=' ', newline='\n')

# plot triangles
fig = plt.figure(figsize=[8,8])                                
ax = fig.add_axes([0.15,0.15,0.8,0.8])                             
ax.set_aspect('equal')

ax.set_xlim([x.min(),x.max()])
ax.set_ylim([y.min(),y.max()])
ax.set_xlabel("lon",fontsize=18)
ax.set_ylabel("lat",fontsize=18)
#plt.savefig("cali_pts.eps")
#plt.show()

for brd in borders:
    xc=brd[:,0];
    yc=brd[:,1];
    plt.plot(xc,yc,'r',lw=2)

for k in range(len(forOut)):
    x1=x[forOut[k,0]]
    x2=x[forOut[k,1]]
    x3=x[forOut[k,2]]
    y1=y[forOut[k,0]]
    y2=y[forOut[k,1]]
    y3=y[forOut[k,2]]
    aaa=ax.add_line(Line2D([x1,x2],[y1,y2]));
    aaa=ax.add_line(Line2D([x2,x3],[y2,y3]));
    aaa=ax.add_line(Line2D([x3,x1],[y3,y1]));

plt.show()
