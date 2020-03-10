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

if [ "$#" -ne 4 ]; then
  echo "usage: $0 nspls iseed fid npar                  "
  echo "        nspls : no. of samples                  "
  echo "        iseed : seed for random number generator"
  echo "        fid   : file ID, typically 1 or 2       "
  echo "        npar  : no of parameters                "
  exit
fi

nspls=$1   # no. of samples
iseed=$2   # seed for random number generation 
fileid=$3  # file ID, typically 1 or 2 for the 2 matrices 
npar=$4    # no. of parameters

# Location of UQTk executables
bindir=../../bin

# Parameter types
# TODO these should be read from files
nrmpar=(1 2) # List of normal     parameters 
lnrpar=()    # List of log-normal parameters
unipar=(3)   # List of uniform    parameters

# (truncated) normal parameters
nrmMu=(0.0 0.2)             # mean values
nrmSg=(0.6 0.5)             # standard deviations
nrmxmi=(-3.14159 -3.14159)  # lower bounds (set -Large number for non-tructated computations)
nrmxma=( 3.14159  3.14159)  # upper bounds (set +Large number for non-tructated computations)

# Log-Normal parameters
lnrMu=()
lnrSg=()
lnrxmi=()
lnrxma=()

# Uniform parameters
unixmi=(-3.14159 -3.14159) # lower bounds
unixma=( 3.14159  3.14159) # upper bounds

# Generate (truncated) normal samples
j=0
for i in ${nrmpar[@]}; do
  if [ $i -le $npar ]; then
    echo "...Parameter ${i} has a truncated normal prior"
    ${bindir}/trdSpls -a ${nrmxmi[j]} -b ${nrmxma[j]} -m ${nrmMu[j]} -s ${nrmSg[j]} -n $nspls -i $iseed -t n
    mv samples.dat par${i}.dat
    j=$(( j + 1 ))
    let "iseed += 13"
  fi
done

# Generate (truncated) log-normal samples
j=0
for i in ${lnrpar[@]}; do
  if [ $i -le $npar ]; then
    echo "...Parameter ${i} has a truncated log-normal prior"
    ${bindir}/trdSpls -a ${lnrxmi[j]} -b ${lnrxma[j]} -m ${lnrMu[j]} -s ${lnrSg[j]} -n $nspls -i $iseed -t ln
    mv samples.dat par${i}.dat
    j=$(( j + 1 ))
    let "iseed += 13"
  fi
done

# Generate uniform samples
j=0
for i in ${unipar[@]}; do
  if [ $i -le $npar ]; then
    echo "...Parameter ${i} has an uniform prior"
    ${bindir}/trdSpls -a ${unixmi[j]} -b ${unixma[j]} -n $nspls -i $iseed -t u
    mv samples.dat par${i}.dat
    j=$(( j + 1 ))
    let "iseed += 13"
  fi
done

# assemble samples in a matrix, each column represents one parameter
/bin/cp par1.dat sgsa.dat; /bin/rm -rf par1.dat
for i in `seq 2 $npar` 
do
  paste sgsa.dat par${i}.dat > tmp.dat; /bin/rm -rf par${i}.dat
  /bin/mv tmp.dat sgsa.dat
done
/bin/mv sgsa.dat sgsa${fileid}.dat

