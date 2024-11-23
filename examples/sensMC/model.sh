#!/bin/bash
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
#     Questions? Contact the UQTk Developers at https://github.com/sandialabs/UQTk/discussions
#     Sandia National Laboratories, Livermore, CA, USA
#=====================================================================================
# This script evaluates the Ishigami function:
# sin(X1) + a sin^2(X2) + b X3^4 sin(X1)
# where X1,X2,X3 are independent RVs in [-pi,pi]

# Default parameters
a="7.0"
b="0.05"
fin="spls.dat"
fout="meval.dat"

# parse command-line arguments
while [[ $# -ge 1 ]]
do
key="$1"
case $key in
    -a)
    a=$2;    shift
    ;;
    -b)
    b=$2;    shift
    ;;
    -i)
    fin=$2;  shift
    ;;
    -o)
    fout=$2; shift
    ;;
    *)
    # unknown option
    echo "Unknown key: ${key}"
    exit
    ;;
esac
shift
done


echo "Evaluating Ishigami function with a=$a and b=$b"
echo "  - input  : $fin"
echo "  - output : $fout"

awk -v a1=$a -v b1=$b '{printf("%12.8e\n",sin($1)+a1*(sin($2))^2+b1*($3)^4*sin($1))}' $fin > $fout 

