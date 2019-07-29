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
  echo "Usage : $0 -p <ON/OFF> -h"
  exit
}

pyintf="OFF"

while getopts ":p:h" opt; do
  case $opt in
    p) pyintf="$OPTARG"
    ;;
    h) usage
    ;;
    \?) echo "Invalid option -$OPTARG" >&2; usage
    ;;
  esac
done

echo "============================================"
echo "Compiling UQTk with:"
echo " - python interface $pyintf"
echo "============================================"

cmake -DCMAKE_INSTALL_PREFIX:PATH=$PWD/../UQTk-install     \
      -DCMAKE_Fortran_COMPILER=/usr/share/Modules/apps/gcc/4.8.2/bin/gfortran \
      -DCMAKE_C_COMPILER=/usr/share/Modules/apps/gcc/4.8.2/bin/gcc \
      -DCMAKE_CXX_COMPILER=/usr/share/Modules/apps/gcc/4.8.2/bin/g++ \
      -DPyUQTk=${pyintf} \
      ../UQTk
