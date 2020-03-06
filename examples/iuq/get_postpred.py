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
#=====================================================================================

import argparse
import numpy as np
from PyUQTk.utils.func import func

##########################################################################
##########################################################################

# Parse input arguments
usage_str = 'Script to evaluate surrogates given samples of the input'
parser = argparse.ArgumentParser(description=usage_str)
parser.add_argument("-p", "--input", dest="input_file",
                    type=str, default='input.dat', help="Parameter input file of size Nxd")
parser.add_argument("-f", "--output", dest="output_file",
                    type=str, default='output.dat', help="Output file name, the size will be NxL")
parser.add_argument("-n", "--nout", dest="nout",
                    type=int, default=1, help="Number of outputs, L")

args = parser.parse_args()

input_file = args.input_file
output_file = args.output_file
nout = args.nout

# Parameter input file of size Nxd
pinput = np.loadtxt(input_file, ndmin=2)
ninput = pinput.shape[0]
dim = pinput.shape[1]

# Hardwired, but this is surrogate type which is always LU
# Do not confuse this with the embedded PC type
pctype = 'LU'

# Output container
outputs = np.zeros((ninput, nout))

# For each output, read corresponding multiindex and
# PC coefficients and evaluate the PC
for j in range(nout):
    mindex = np.loadtxt('mindexp.' + str(j) + '_pred.dat',
                        dtype=int).reshape(-1, dim)
    pcf = np.loadtxt('pccfp.' + str(j) + '_pred.dat').reshape(-1, 1)
    outputs[:, j] = func(pinput, 'PCmi', [mindex, pcf, pctype])

# Save the output in the requested file
np.savetxt(output_file, outputs)
