#!/bin/bash
#=====================================================================================
#
#                      The UQ Toolkit (UQTk) version 3.1.4
#                          Copyright (2023) NTESS
#                        https://www.sandia.gov/UQToolkit/
#                        https://github.com/sandialabs/UQTk
#
#     Copyright 2023 National Technology & Engineering Solutions of Sandia, LLC (NTESS).
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
#     Questions? Contact the UQTk Developers at https://github.com/sandialabs/UQTk/discussions
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

cmake -DCMAKE_INSTALL_PREFIX:PATH=$PWD/../UQTk-install     \
      -DCMAKE_Fortran_COMPILER=/opt/intel/fc/Compiler/11.1/080/bin/ia64/ifort \
      -DCMAKE_C_COMPILER=/opt/intel/cc/Compiler/11.1/080/bin/ia64/icc \
      -DCMAKE_CXX_COMPILER=/opt/intel/cc/Compiler/11.1/080/bin/ia64/icpc \
      -DIntelLibPath=/opt/intel/fc/Compiler/11.1/080/lib/ia64 \
      -DPyUQTk=OFF \
      ../UQTk
