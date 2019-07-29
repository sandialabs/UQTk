#!/bin/bash
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

usage ()
{
  echo "Usage : $0 -c <gnu/intel> -p <ON/OFF> -h"
  exit
}

pyintf="OFF"
ctype="gnu"

while getopts ":p:c:h" opt; do
  case $opt in
    p) pyintf="$OPTARG"
    ;;
    c) ctype="$OPTARG"
    ;;
    h) usage
    ;;
    \?) echo "Invalid option -$OPTARG" >&2; usage
    ;;
  esac
done

echo "============================================"
echo "Compiling UQTk with:"
echo " - $ctype compilers"
echo " - python interface $pyintf"
echo "============================================"


if [ "${ctype}" == "gnu" ]; then
cmake -DCMAKE_INSTALL_PREFIX:PATH=$PWD/../UQTk-install     \
      -DCMAKE_Fortran_COMPILER=gfortran \
      -DCMAKE_C_COMPILER=gcc            \
      -DCMAKE_CXX_COMPILER=g++          \
      -DPyUQTk=${pyintf} \
      ../UQTk
elif [ "${ctype}" == "intel" ]; then
cmake -DCMAKE_INSTALL_PREFIX:PATH=$PWD/../UQTk-install \
      -DCMAKE_Fortran_COMPILER=ifort                   \
      -DCMAKE_C_COMPILER=icc                           \
      -DCMAKE_CXX_COMPILER=icpc                        \
      -DPyUQTk=${pyintf}                               \
      ../UQTk
else
  echo "Unknown compiler: ${ctype}"
fi
