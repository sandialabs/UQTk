#!/bin/bash -e
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

#=====================================================================================

# Script location
export UQPC=${UQTK_INS}/examples/uqpc

# Parameter range file
echo "0 5" > param_range.txt
echo "0 2.5" >> param_range.txt
echo "-2 5" >> param_range.txt

# Prepare inputs
${UQPC}/uq_pc.py -r offline_prep -p param_range.txt -s rand -n 111 -v 33

# ptrain.dat is NTRAIN x DIM matrix, each row is a DIM-variate parameter sample
# qtrain.dat is the same scaled to [-1,1]
# wtrain.dat are quadrature weights only if sampling method is quadrature
# ytrain.dat is NTRAIN x NOUT vector of outputs
# pval.dat is NVAL x DIM matrix, each row is a DIM-variate parameter sample
# qval.dat is the same scaled to [-1,1]
# yval.dat is NVAL x NOUT vector of outputs

# Run the black-box model
${UQPC}/model.py -i ptrain.dat -o ytrain.dat
${UQPC}/model.py -i pval.dat -o yval.dat

# Build surrogates
${UQPC}/uq_pc.py -r offline_post -p param_range.txt -m lsq -s rand -n 111 -v 3 -t 3

rm -f pnames.txt outnames.txt
# Plot model-vs-surrogate
${UQPC}/plot.py dm training  validation
# Plot total sensitivities
${UQPC}/plot.py sens total

# Enjoy the .eps files!
