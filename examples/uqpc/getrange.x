#!/bin/bash
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

# Given a Nxd matrix of samples, compute dx2 ranges
# Example: getrange.x samples.dat > ranges.dat

shopt -s expand_aliases
alias awk="awk -v OFMT='%.15e'"

if [ $# -lt 1 ]; then
    echo "Number of arguments can not be less than 1"
    echo "Syntax: $0 <filename> [<cushion_fraction>]"
    exit
elif [ $# -gt 2 ]; then
    echo "Number of arguments can not be greater than 2"
    echo "Syntax: $0 <filename> [<cushion_fraction>]"
    exit
elif [ $# -eq 1 ]; then
    fr=0.0
elif [ $# -eq 2 ]; then
    fr=$2
fi

filename=$1

DIM=`awk 'END{print NF}' $filename`

for (( COL=1; COL<=$DIM ; COL++ )); do

awk 'BEGIN {
valmin=1e+100;
valmax=-1e+100;
line=1;
}
{
if( $(col) < valmin )
{
  valmin=$(col);
}
if( $(col) > valmax )
{
  valmax=$(col);
}
}
END{
print valmin-fr*(valmax-valmin), valmax+fr*(valmax-valmin)
}' col=$COL $filename
done
