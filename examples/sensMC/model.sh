#!/bin/bash
# =====================================================================================
#                      The UQ Toolkit (UQTk) version 3.0.4
#                     Copyright (2017) Sandia Corporation
#                      http://www.sandia.gov/UQToolkit/
#
#     Copyright (2017) Sandia Corporation. Under the terms of Contract DE-AC04-94AL85000
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
# =====================================================================================
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

