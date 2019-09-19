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
  echo "No command-line parameters, execute as is"
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

# If UQTk is installed in a directory with the version number after "UQTk", then define the suffix here
# (The test below may seem a bit convoluted, but basically it checks to see if the version number
# has a "." in it, which differentiates it from the placeholder string we have in place for
# our development branches of the code)
UQTK_VERSION=@UQTKVERSION@
if [[ "$UQTK_VERSION" == *.* ]]; then
  VERSION_SUFFIX=_v$UQTK_VERSION
fi

UQTK_SRC_DIR=$PWD/../UQTk$VERSION_SUFFIX
UQTK_INSTALL_DIR=$UQTK_SRC_DIR-install

echo "This script assumes the UQTk source code is in $UQTK_SRC_DIR"
echo "and will be installed in $UQTK_INSTALL_DIR"

cmake -DCMAKE_INSTALL_PREFIX:PATH=$UQTK_INSTALL_DIR    \
      -DCMAKE_Fortran_COMPILER=gfortran \
      -DCMAKE_C_COMPILER=gcc            \
      -DCMAKE_CXX_COMPILER=g++          \
      -DPYTHON_EXECUTABLE:FILEPATH=/opt/local/bin/python \
      -DPYTHON_LIBRARY:FILEPATH=/opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/libpython2.7.dylib \
      -DPyUQTk=ON \
      $UQTK_SRC_DIR
