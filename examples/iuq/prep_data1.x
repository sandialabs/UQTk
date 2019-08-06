#!/bin/bash -e
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
#     Questions? Contact Bert Debusschere <bjdebus@sandia.gov>
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
