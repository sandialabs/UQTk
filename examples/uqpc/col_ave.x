#!/bin/sh -e
#=====================================================================================
#
#                      The UQ Toolkit (UQTk) version 3.1.2
#                          Copyright (2022) NTESS
#                        https://www.sandia.gov/UQToolkit/
#                        https://github.com/sandialabs/UQTk
#
#     Copyright 2022 National Technology & Engineering Solutions of Sandia, LLC (NTESS).
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
# Finds average of each column and the std

alias awk="awk -v OFMT='%.15e'"

#Input
INFILE=$1
OUTFILE=$2

# Dim check
N=`awk 'END{print NR}' $INFILE`
M=`awk 'END{print NF}' $INFILE`

echo "Make sure the data file has rectangular matrix form"

DIRF=`which $0`
DIR=`dirname $DIRF`


awk '{for (i=1; i<=NF; i++) { sum[i]+=$i; ssum[i]+=$i*$i}}; END{for (i=1; i<=NF; i++) { print sum[i]/n, sqrt(ssum[i]/n-sum[i]*sum[i]/(n*n)) }}' n=$N $INFILE | $DIR/transpose_file.x > $OUTFILE
