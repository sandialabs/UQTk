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

try: 
    import numpy as npy
except ImportError:
    print "Numpy was not found. "

def CRPSinteg(s1,s2):
    """ Computes integral of squared difference between two CDFs """
    Ns1 = s1.shape[0];
    Ns2 = s2.shape[0];
    ds1 = 1.0/Ns1;
    ds2 = 1.0/Ns2;
    # Combine samples and sort
    s12=npy.sort(npy.concatenate((s1,s2)));
    CRPS = 0.0;
    j1 = 0;
    j2 = 0;
    for i in range(Ns1+Ns2-1):
        if s12[i+1] <= s1[0]:
            fs1 = 0.0;
        elif s12[i] >= s1[Ns1-1]:
            fs1 = 1.0;
        else:
            j1 = j1+npy.argmax(s1[j1:]>s12[i])-1
            fs1 = (j1+1)*ds1;
        if s12[i+1] <= s2[0]:
            fs2 = 0.0;
        elif s12[i] >= s2[Ns2-1]:
            fs2 = 1.0;
        else:
            j2 = j2+npy.argmax(s2[j2:]>s12[i])-1
            fs2 = (j2+1)*ds2;
        CRPS = CRPS + (s12[i+1]-s12[i])*(fs1-fs2)**2;
    return CRPS

def CRPS(s1,s2):
    """ Computes CRPS score """
    nsamples = s1.shape[0]
    if nsamples != s2.shape[0]:
        print "The number of realizations in s1 and s2 is not the same:",nsamples,s2.shape[0]
        return (-1.0);
    crps = npy.zeros(nsamples)
    for i in range(nsamples):
        crps[i] =  CRPSinteg(s1[i],s2[i])
    return crps.mean()
