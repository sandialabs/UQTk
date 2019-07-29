#=====================================================================================
#                     The UQ Toolkit (UQTk) version @UQTKVERSION@
#                     Copyright (@UQTKYEAR@) Sandia Corporation
#                     http://www.sandia.gov/UQToolkit/
#
#     Copyright (@UQTKYEAR@) Sandia Corporation. Under the terms of Contract DE-AC04-94AL85000
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
import numpy as np

def readfile(filename):
    d0  = []
    nlines = 0
    for line in file(filename):
        line = line.rstrip('\n')
        line_list = [float(x) for x in line.split()]
        d0.append(line_list)
        nlines = nlines+1
    return d0,nlines

def column(matrix, i):
    return [row[i] for row in matrix]

def checkPtInside(xc,yc,x,y):
    nc=len(xc);
    nm=len(x);
    umask=np.array(np.zeros(nm))
    for i in range(nm):
        b=y[i]-x[i];
        nint = 0;
        for j in range(nc-1):
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

