#!/bin/bash -e
#=====================================================================================
#
#                      The UQ Toolkit (UQTk) version 3.1.1
#                          Copyright (2021) NTESS
#                        https://www.sandia.gov/UQToolkit/
#                        https://github.com/sandialabs/UQTk
#
#     Copyright 2021 National Technology & Engineering Solutions of Sandia, LLC (NTESS).
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
#=====================================================================================

# Script-example for 2d data generation


UQPC=${UQTK_INS}/examples/uqpc
UQBIN=${UQTK_INS}/bin

# Sample a random parameter vector and scale to range
${UQBIN}/pce_rv -w PCvar -x LEG -d 3 -p 3 -n 1 -s 11111
${UQPC}/scale.x rvar.dat to prange.dat params_true.dat

# Eval the model with some noise corruption
${UQPC}/model.py -i params_true.dat -o ydata_all.tmp -x xdata_all.txt -m ex_xp -r 0.05 -f 1.0
${UQPC}/transpose_file.x ydata_all.tmp > ydata_all.txt

# let's take surrogate errors as placeholder for data noise
cp surr_errors.dat dataerr_all.dat

# let's use data from all x-conditions for calibration
awk '{print 1}' xdata_all.txt > ind_calib.dat


# name the x-conditions for plotting convenience (optional),
# this should have as many rows as the number of columns in xdata_all.txt
# e.g. data is collected for various combinations of T and P
echo "T" > xcond_names.txt
echo "P" >> xcond_names.txt
