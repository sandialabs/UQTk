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

try:
    import numpy as npy
except ImportError:
    print('Numpy was not found.')

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
        print('The number of realizations in s1 and s2 is not the same: %d vs %d'%(nsamples,s2.shape[0]))
        return (-1.0);
    crps = npy.zeros(nsamples)
    for i in range(nsamples):
        crps[i] =  CRPSinteg(s1[i],s2[i])
    return crps.mean()
