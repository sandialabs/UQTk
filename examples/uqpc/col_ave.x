#!/bin/sh -e
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
