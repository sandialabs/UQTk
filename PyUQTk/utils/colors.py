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

"""
Utilities for defining color lists.
Note that Python has standards as well; this is an alternative
that often produces enough variability in colors for eye-pleasing results.
"""

try:
    import numpy as np
except ImportError:
    print('Numpy was not found.')



def set_colors(npar):
    """ Sets a list of different colors of requested length, as rgb triples"""
    colors = []
    pp=1+int(npar/6)
    for i in range(npar):
        c=1-float(int((i/6))/pp)
        b=np.empty((3))
        for jj in range(3):
            b[jj]=c*int(i%3==jj)
        a=int(int(i%6)/3)
        colors.append(((1-a)*b[2]+a*(c-b[2]),(1-a)*b[1]+a*(c-b[1]),(1-a)*b[0]+a*(c-b[0])))

    return colors
