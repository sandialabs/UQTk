#!/bin/bash
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
