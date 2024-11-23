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

import os
import sys 
import numpy as np
from   scipy import stats
import matplotlib
import matplotlib.pylab as plt

def ReadDataFile(dataFileName):
    """Read records in file with name dataFileName into
    a numpy matrix of floats"""

    # Open text file with all data
    dataFile = open(dataFileName,"r")

    # Convert all data lines into floats and put them in a list
    dataList = []
    for line in dataFile.readlines():
     records = line.split()
     numRecords = [float(s) for s in records]
     dataList.append(numRecords)

    dataFile.close()

    # Convert list to 2D numpy matrix object
    data = np.array(dataList)

    return data

def SetPlotParams():
    """ Set line widths and font sizes.
    Can be expanded down the road to take
    an argument like paper or slide. 
    Returns linewidth and fontsize. """

    linewidth = 2;
    fontsize = 16;

    return linewidth, fontsize
