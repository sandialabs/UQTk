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


# 3 options for input parameter PC construction. Uncomment one of them.

## (a) Given PC explicitely
# INPC_ORDER=1
# echo "2 -1 3" > param_pcf.txt
# echo "1 0 0" >> param_pcf.txt
# echo "0 1 0" >> param_pcf.txt
# echo "0 0 1" >> param_pcf.txt

## (b) Given marginal, independent PCs per dimension
# echo "1 0.3 " > param_margpc.txt
# echo "3 0.1 0.2 0.1" >> param_margpc.txt
# echo "1 0.5 0.2" >> param_margpc.txt
# INPC_ORDER=3
# ${UQPC}/prepare_inpc.py marg param_margpc.txt $INPC_ORDER

## (c) Given samples of inputs (e.g. from a prior calibration study)
${UQPC}/generate_inputsamples.py
INPC_ORDER=3
${UQPC}/prepare_inpc.py sam param_sam.txt $INPC_ORDER


DIM=`awk 'NR==1{print NF}' param_pcf.txt`

# Prepare inputs
${UQPC}/uq_pc.py -r offline_prep -c param_pcf.txt -x HG -d $DIM -o $INPC_ORDER -m lsq -s rand -n 100 -v 0

# Run the black-box model,
# (a) Test function with 3 inputs and 5 outputs
awk '{print exp($1), $3*log($1**2)+$1+($2**2)*(1-exp(-$3**2)), $1+$3**2, $2+$3, $1*$3}' ptrain.dat > ytrain.dat

# (b) entry-wise exponential!
#awk '{for (i=1; i<=NF; i++) { printf("%lg ",exp($i)) }; printf("\n")}' ptrain.dat > ytrain.dat

# Build surrogates
${UQPC}/uq_pc.py -r offline_post -c param_pcf.txt -x HG -d $DIM -o $INPC_ORDER -m lsq -s rand -n 100 -v 0 -t 5


rm -f pnames.txt outnames.txt
# Plot output pdfs
${UQPC}/plot.py pdf
# Plot main sensitivities (barplot)
${UQPC}/plot.py sens main

# Plot total sensitivities (matrix-plot)
${UQPC}/plot.py sensmat total

# Plot sensitivities (circular plot)
${UQPC}/plot.py senscirc

# Enjoy the .eps files!