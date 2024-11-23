#!/bin/bash -e
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
#     Questions? Contact the UQTk Developers at <uqtk-developers@software.sandia.gov>
#     Sandia National Laboratories, Livermore, CA, USA
#=====================================================================================

#=====================================================================================

# Script location
export UQPC=${UQTK_INS}/examples/uqpc

###########################################################
## Assume input.dat (Nxd) and output.dat (Nxr) are given ##
###########################################################

# Use all inputs as training, and no validation
cp ${UQPC}/input.dat ptrain.dat
cp ${UQPC}/output.dat ytrain.dat

# Get ranges of inputs, with 10% 'cushion' from the dimension-wise extreme samples
${UQPC}/getrange.x ptrain.dat 0.1 > param_range.txt

# Scale the inputs
${UQPC}/scale.x ptrain.dat from param_range.txt qtrain.dat

# Get the number of training samples
NSAM=`echo | awk 'END{print NR}' ptrain.dat`

# Build surrogates
${UQPC}/uq_pc.py -r offline_post -p param_range.txt -m lsq -s rand -n $NSAM -v 0 -t 3

rm -f pnames.txt outnames.txt
# Plot model-vs-surrogate
${UQPC}/plot.py dm training
# Plot total sensitivities
${UQPC}/plot.py sens total

# Enjoy the .eps files!
