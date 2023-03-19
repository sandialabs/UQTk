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
import numpy as np

def readfile(filename):
    file1 = open(filename, 'r')
    Lines = file1.readlines()
    d0  = []
    nlines = 0
    for line in Lines:
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

