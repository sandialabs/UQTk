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
  echo "No command-line parameters, run as is from a build directory"
  exit
}

while getopts ":h" opt; do
  case $opt in
    h) usage
    ;;
    \?) echo "Invalid option -$OPTARG" >&2; usage
    ;;
  esac
done

# If UQTk is installed in a directory other than "UQTk", then modify the line below accordingly
UQTK_SRC_DIR=$PWD/../UQTk
UQTK_INSTALL_DIR=$UQTK_SRC_DIR-install

echo "This script assumes the UQTk source code is in $UQTK_SRC_DIR"
echo "and will be installed in $UQTK_INSTALL_DIR"

cmake -DCMAKE_INSTALL_PREFIX:PATH=$UQTK_INSTALL_DIR    \
      -DCMAKE_Fortran_COMPILER=gfortran                \
      -DCMAKE_C_COMPILER=gcc                           \
      -DCMAKE_CXX_COMPILER=g++                         \
      $UQTK_SRC_DIR
