#!/bin/bash -e
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
#=====================================================================================

# Script-example for 1d model ensemble generation

# Directory forward UQ scripts
UQPC=${UQTK_INS}/examples/uqpc

# Create parameter ranges file
echo "-1.0 2.0" > prange.dat
echo "-1.0 4.0" >> prange.dat
echo "-2.0 1.0" >> prange.dat

NTRAIN=100  # N
NVAL=33     # V
SAMPLING=rand

# Prepare inputs
${UQPC}/uq_pc.py -r offline_prep -p prange.dat -s $SAMPLING -n $NTRAIN -v $NVAL


###################################################################################
###################################################################################

# Test function with 3 inputs and 4 outputs
awk '{print $3*log($1**2)+$1+($2**2)*(1-exp(-$3**2)), $1+$3**2, $2+$3, $1*$3+$2}' ptrain.dat > ytrain.dat
awk '{print $3*log($1**2)+$1+($2**2)*(1-exp(-$3**2)), $1+$3**2, $2+$3, $1*$3+$2}' pval.dat   > yval.dat

# Create input param names (optional)
awk '{print "par"NR}' prange.dat > pnames.txt
# Create output names (optional)
echo "Qoi1" > outnames.txt
echo "Qoi2" >> outnames.txt
echo "Qoi3" >> outnames.txt
echo "Qoi4" >> outnames.txt

NOUT=`awk 'NR==1{print NF}' ytrain.dat`

# Design conditions x
# In this case xdata_all.txt simply IDs the conditions, say, from 1 to 4
awk '{print NR}' outnames.txt > xdata_all.txt
