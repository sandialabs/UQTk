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

# Script-example for 2d model ensemble generation

# Directory forward UQ scripts
UQPC=${UQTK_INS}/examples/uqpc
UQBIN=${UQTK_INS}/bin

# Create parameter ranges file
echo "-1.0 0.5" > prange.dat
echo "1.0 2.0" >> prange.dat
echo "-2.0 1.0" >> prange.dat

NTRAIN=111  # N
NVAL=33     # V
SAMPLING=rand

# Prepare inputs
${UQPC}/uq_pc.py -r offline_prep -p prange.dat -s $SAMPLING -n $NTRAIN -v $NVAL


###################################################################################
###################################################################################

# Create 10x7 grid of points as x-conditions
${UQBIN}/generate_quad -d 2 -g NC -x full -p 10
head -n90 qdpts.dat | tail -n70 > xdata_all.txt

# Test function
${UQPC}/model.py -i ptrain.dat -o ytrain.dat -x xdata_all.txt -m ex_xp
${UQPC}/model.py -i pval.dat -o yval.dat -x xdata_all.txt -m ex_xp


# Create input param names (optional)
awk '{print "p"NR}' prange.dat > pnames.txt
# Create output names (optional)
awk '{print "out"NR}' xdata_all.txt > outnames.txt
