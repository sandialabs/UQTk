#!/bin/bash

BINPATH="../scripts"
OUTPUT=output.txt
INPUTXML=input.xml

echo "Getting data"

$BINPATH/get_data.py --ix $INPUTXML -g -e

echo "Fitting model to data"

$BINPATH/fit.py --ix $INPUTXML -w $OUTPUT

echo "Preforming Postprocessing"
$BINPATH/post.py --ix $INPUTXML --evidence -d -v 3
