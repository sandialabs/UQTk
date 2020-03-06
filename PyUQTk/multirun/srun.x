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
# This script is run automatically via multirun.py


# Get the script name (first entry of args.in)
SCRIPT=`cut -f 1 -d" " args.in`
# Get the output file name to dump the screen-output (second entry of args.in)
OUT=`cut -f 2 -d" " args.in`
# The rest of entries in args.in are parameters of the script
ARGUM=`cut -f 3- -d" " args.in`

# Informational print
THIS=`basename $PWD`
echo "Running $SCRIPT $ARGUM > $OUT in $THIS"

##echo $(< args.in)

# Running the script
cd ..
SCRIPT_ABS=`echo "$(cd "$(dirname "$SCRIPT")"; pwd)/$(basename "$SCRIPT")"`
cd -
ln -sf $SCRIPT_ABS linkToScript
./linkToScript $ARGUM > $OUT
