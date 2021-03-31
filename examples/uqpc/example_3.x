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
echo "-1 2" > param_range.txt
echo "0 3" >> param_range.txt
echo "2 5" >> param_range.txt
DIM=`awk 'END{print NR}' param_range.txt`

# Prepare inputs
${UQPC}/uq_pc.py -r offline_prep -p param_range.txt -s rand -n 222 -v 0

# Black-box model evaluator
NGRID=100


echo -n"" > ytrain.dat
for ((tgrid=1;tgrid<=$NGRID;tgrid++)); do
    awk '{print $3*$1**2+$2*exp((tgr/ngr)*$3)}' tgr=$tgrid ngr=$NGRID ptrain.dat > tmp
    paste ytrain.dat tmp > tmpp; mv tmpp ytrain.dat; rm -f tmp
done

# Build surrogates
${UQPC}/uq_pc.py -r offline_post -p param_range.txt -m lsq -s rand -n 222 -v 0 -t 3


rm -f pnames.txt outnames.txt
for ((i=1;i<=$NGRID;i++)); do
    echo "$i" >> outnames.txt
done

# Plot total sensitivities
${UQPC}/plot.py sens total

# Enjoy the .eps files!
