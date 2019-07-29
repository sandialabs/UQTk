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
#     Need help with UQTk? Check out the resources on http://www.sandia.gov/UQToolkit/
#     or e-mail uqtk-users@software.sandia.gov
#     (subscription details listed at http://www.sandia.gov/UQToolkit/)
#     Other questions? Contact Bert Debusschere <bjdebus@sandia.gov>
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

UQTK_SRC_DIR=/Users/jlao/UQTk
UQTK_INSTALL_DIR=/Users/jlao/UQTk-install

echo "This script assumes the UQTk source code is in $UQTK_SRC_DIR"
echo "and will be installed in $UQTK_INSTALL_DIR"

cmake -DCMAKE_INSTALL_PREFIX:PATH=$UQTK_INSTALL_DIR    \
      -DCMAKE_Fortran_COMPILER=gfortran                \
      -DCMAKE_C_COMPILER=gcc                           \
      -DCMAKE_CXX_COMPILER=g++                         \
      -DHDF5LibPath=/opt/local/hdf5                    \
      $UQTK_SRC_DIR
