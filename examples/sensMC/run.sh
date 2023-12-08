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
#     Questions? Contact the UQTk Developers at https://github.com/sandialabs/UQTk/discussions
#     Sandia National Laboratories, Livermore, CA, USA
#=====================================================================================

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
