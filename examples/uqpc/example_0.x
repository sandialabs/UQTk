#!/bin/bash -e
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
#=====================================================================================

# Script location
export UQPC=${UQTK_INS}/examples/uqpc

###########################################################
## Assume input.dat (Nxd) and output.dat (Nxr) are given ##
###########################################################

# Use all inputs as training, and no validation
cp input.dat ptrain.dat
cp output.dat ytrain.dat

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