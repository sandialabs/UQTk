#!/bin/bash -e
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
#     Questions? Contact the UQTk Developers at <uqtk-developers@software.sandia.gov>
#     Sandia National Laboratories, Livermore, CA, USA
#=====================================================================================
#=====================================================================================

# Script-example for 1d data generation

# ydata_all.txt has 4 rows, one for each output,
# and k columns; in this case, k=1, i.e. one measurement per condition
echo "10" > ydata_all.txt
echo "3" >> ydata_all.txt
echo "4" >> ydata_all.txt
echo "3" >> ydata_all.txt

# let's take surrogate errors as placeholder for data noise
cp surr_errors.dat dataerr_all.dat

# let's use data from all x-conditions for calibration
awk '{print 1}' xdata_all.txt > ind_calib.dat

# generate_quad -d1 -g NC -x full -p 111
# mv qdpts.dat xdata_all.txt
# awk '{print 3+$1-$1*$1}' xdata_all.txt > ydata_all.txt


# name the x-conditions for plotting convenience,
# this should have as many rows as the number of columns in xdata_all.txt
echo "QoI Id" > xcond_names.txt
