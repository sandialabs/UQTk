#!/bin/bash
#=====================================================================================
#                     The UQ Toolkit (UQTk) version 3.0
#                     Copyright (2015) Sandia Corporation
#                     http://www.sandia.gov/UQToolkit/
#
#     Copyright (2015) Sandia Corporation. Under the terms of Contract DE-AC04-94AL85000
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
  echo "Usage : $0 -d <UQTk_DIR> -c <gnu53/gnu62/gnu63/gnu71> -p <ON/OFF> -h"
  exit
}

uqtksrc="${HOME}/Projects/UQTk/3.0/gitdir"
uqtkins="${HOME}/local/uqtk"
pyintf="OFF"
ctype="gnu62"

while getopts ":p:c:d:h" opt; do
  case $opt in
    p) pyintf="$OPTARG"
    ;;
    c) ctype="$OPTARG"
    ;;
    d) uqtksrc="$OPTARG"
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

#PATH2MUQ=${HOME}/Projects/muq-install
GNUROOT=/opt/local

if [ "${ctype}" == "gnu53" ]; then
cmake -DCMAKE_INSTALL_PREFIX:PATH=${uqtkins}/gnu53 \
      -DCMAKE_Fortran_COMPILER=/usr/local/opt/gcc53/bin/gfortran-5.3.0 \
      -DCMAKE_C_COMPILER=/usr/local/opt/gcc53/bin/gcc-5.3.0 \
      -DCMAKE_CXX_COMPILER=/usr/local/opt/gcc53/bin/g++-5.3.0 \
      -DPATH2MUQ=${PATH2MUQ} \
      -DPyUQTk=${pyintf} \
      ${uqtksrc}
elif [ "${ctype}" == "gnu62" ]; then
cmake -DCMAKE_INSTALL_PREFIX:PATH=${uqtkins}/gnu62 \
      -DCMAKE_Fortran_COMPILER=$GNUROOT/gcc62/bin/gfortran-6.2.0 \
      -DCMAKE_C_COMPILER=$GNUROOT/gcc62/bin/gcc-6.2.0 \
      -DCMAKE_CXX_COMPILER=$GNUROOT/gcc62/bin/g++-6.2.0 \
      -DPATH2MUQ=${PATH2MUQ} \
      -DPyUQTk=${pyintf} \
      ${uqtksrc}
elif [ "${ctype}" == "gnu63" ]; then
cmake -DCMAKE_INSTALL_PREFIX:PATH=${uqtkins}/gnu63 \
      -DCMAKE_Fortran_COMPILER=$GNUROOT/gcc63/bin/gfortran-6.3.0 \
      -DCMAKE_C_COMPILER=$GNUROOT/gcc63/bin/gcc-6.3.0 \
      -DCMAKE_CXX_COMPILER=$GNUROOT/gcc63/bin/g++-6.3.0 \
      -DPATH2MUQ=${PATH2MUQ} \
      -DPyUQTk=${pyintf} \
      ${uqtksrc}
elif [ "${ctype}" == "gnu71" ]; then
cmake -DCMAKE_INSTALL_PREFIX:PATH=${uqtkins}/gnu71 \
      -DCMAKE_Fortran_COMPILER=$GNUROOT/gcc71/bin/gfortran-7.1.0 \
      -DCMAKE_C_COMPILER=$GNUROOT/gcc71/bin/gcc-7.1.0 \
      -DCMAKE_CXX_COMPILER=$GNUROOT/gcc71/bin/g++-7.1.0 \
      -DPATH2MUQ=${PATH2MUQ} \
      -DPyUQTk=${pyintf} \
      ${uqtksrc}
else
  echo "Unknown compiler: ${ctype}"
fi

