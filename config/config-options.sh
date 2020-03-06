#!/bin/bash
#=====================================================================================
#
#                      The UQ Toolkit (UQTk) version @UQTKVERSION@
#                          Copyright (@UQTKYEAR@) NTESS
#                        https://www.sandia.gov/UQToolkit/
#                        https://github.com/sandialabs/UQTk
#
#     Copyright @UQTKYEAR@ National Technology & Engineering Solutions of Sandia, LLC (NTESS).
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

usage ()
{
  echo "Usage : $0 -d <UQTk_DIR> -c <gnu53/gnu61/gnu61m/gnu71/intel/clang> -p <ON/OFF> -h"
  exit
}

uqtksrc="${HOME}/Projects/UQTk/3.0/gitdir"
pyintf="OFF"
ctype="gnu71"

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
cmake -DCMAKE_INSTALL_PREFIX:PATH=$PWD/../install_5.3 \
      -DCMAKE_Fortran_COMPILER=/usr/local/opt/gcc53/bin/gfortran-5.3.0 \
      -DCMAKE_C_COMPILER=/usr/local/opt/gcc53/bin/gcc-5.3.0 \
      -DCMAKE_CXX_COMPILER=/usr/local/opt/gcc53/bin/g++-5.3.0 \
      -DPATH2MUQ=${PATH2MUQ} \
      -DPyUQTk=${pyintf} \
      ${uqtksrc}
elif [ "${ctype}" == "gnu61" ]; then
cmake -DCMAKE_INSTALL_PREFIX:PATH=$PWD/../install_6.1 \
      -DCMAKE_Fortran_COMPILER=$GNUROOT/gcc61/bin/gfortran-6.1.0 \
      -DCMAKE_C_COMPILER=$GNUROOT/gcc61/bin/gcc-6.1.0 \
      -DCMAKE_CXX_COMPILER=$GNUROOT/gcc61/bin/g++-6.1.0 \
      -DPATH2MUQ=${PATH2MUQ} \
      -DPyUQTk=${pyintf} \
      ${uqtksrc}
elif [ "${ctype}" == "gnu61m" ]; then
cmake -DCMAKE_INSTALL_PREFIX:PATH=$PWD/../install_6.1m \
      -DCMAKE_Fortran_COMPILER=$GNUROOT/gcc61/bin/mpif90 \
      -DCMAKE_C_COMPILER=$GNUROOT/gcc61/bin/mpicc \
      -DCMAKE_CXX_COMPILER=$GNUROOT/gcc61/bin/mpic++ \
      -DPATH2MUQ=${PATH2MUQ} \
      -DPyUQTk=${pyintf} \
      ${uqtksrc}
elif [ "${ctype}" == "gnu71" ]; then
cmake -DCMAKE_INSTALL_PREFIX:PATH=$PWD/../install_7.1 \
      -DCMAKE_Fortran_COMPILER=$GNUROOT/gcc71/bin/gfortran-7.1.0 \
      -DCMAKE_C_COMPILER=$GNUROOT/gcc71/bin/gcc-7.1.0 \
      -DCMAKE_CXX_COMPILER=$GNUROOT/gcc71/bin/g++-7.1.0 \
      -DPATH2MUQ=${PATH2MUQ} \
      -DPyUQTk=${pyintf} \
      ${uqtksrc}
elif [ "${ctype}" == "intel" ]; then
cmake -DCMAKE_INSTALL_PREFIX:PATH=$PWD/../install_intel \
      -DCMAKE_Fortran_COMPILER=/opt/intel/composerxe/bin/ifort \
      -DCMAKE_C_COMPILER=/opt/intel/composerxe/bin/icc \
      -DCMAKE_CXX_COMPILER=/opt/intel/composerxe/bin/icpc \
      -DPyUQTk=${pyintf} \
      ${uqtksrc}
elif [ "${ctype}" == "clang" ]; then
cmake -DCMAKE_INSTALL_PREFIX:PATH=$PWD/../install_clang \
      -DCMAKE_Fortran_COMPILER=$GNUROOT/gcc61/bin/gfortran-6.1.0 \
      -DCMAKE_C_COMPILER=clang \
      -DCMAKE_CXX_COMPILER=clang++ \
      -DClangLibPath=/opt/local/gcc61/lib \
      -DPyUQTk=${pyintf} \
      ${uqtksrc}
else
  echo "Unknown compiler: ${ctype}"
fi
