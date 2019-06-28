#!/bin/bash
# =====================================================================================
#                      The UQ Toolkit (UQTk) version 3.0.4
#                     Copyright (2017) Sandia Corporation
#                      http://www.sandia.gov/UQToolkit/
#
#     Copyright (2017) Sandia Corporation. Under the terms of Contract DE-AC04-94AL85000
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
# =====================================================================================

nspl=100000    # no. of samples
ndim=3         # dimensionality (no. of parameters)
iseed1=2014    # seed for first set of random samples
iseed2=2015    # seed for second set of random samples
bindir=$PWD/../../bin # location of UQTk binaries

echo "Generate parameter matrices"
./genRndSpls.sh ${nspl} ${iseed1} 1 ${ndim}
./genRndSpls.sh ${nspl} ${iseed2} 2 ${ndim}

#
#    First-order sensitivities
#
echo "Generate MC samples for Si"
${bindir}/sens -a splFO -u sgsa1.dat -v sgsa2.dat -d ${ndim} -n ${nspl} >& genSpl_Si.log

echo "Evaluate model"
./model.sh -b 0.1 -i splFO.txt -o mevalFO.txt                           >& meval_Si.log

echo "Compute Si"
${bindir}/sens -a idxFO -x mevalFO.txt -d ${ndim} -n ${nspl}            >& comp_Si.log

#
#    Total effect sensitivities
#
echo "Generate MC samples for SiT"
${bindir}/sens -a splTO -u sgsa1.dat -v sgsa2.dat -d ${ndim} -n ${nspl} >& genSpl_SiT.log

echo "Evaluate model"
./model.sh  -b 0.1 -i splTO.txt -o mevalTO.txt                          >& genSpl_SiT.log

echo "Compute SiT"
${bindir}/sens -a idxTO -x mevalTO.txt -d ${ndim} -n ${nspl}            >& comp_SiT.log

#
#    Joint effect sensitivities
#
echo "Generate MC samples for Sij"
${bindir}/sens -a splJnt -u sgsa1.dat -v sgsa2.dat -d ${ndim} -n ${nspl} >& genSpl_Sij.log

echo "Evaluate model"
./model.sh  -b 0.1 -i splJ.txt -o mevalJ.txt                             >& genSpl_Sij.log

echo "Compute Sij"
${bindir}/sens -a idxJnt -x mevalJ.txt -d ${ndim} -n ${nspl}             >& comp_Sij.log
