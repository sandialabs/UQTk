#!/bin/bash
#=====================================================================================
#
#                      The UQ Toolkit (UQTk) version 3.1.0
#                          Copyright (2020) NTESS
#                        https://www.sandia.gov/UQToolkit/
#                        https://github.com/sandialabs/UQTk
#
#     Copyright 2020 National Technology & Engineering Solutions of Sandia, LLC (NTESS).
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
# Parallel evaluation of likelihood for TMCMC

nProc=$1 # number of processes
SLP1=0.05 # pause between checking existance of "done" files

# separate processes are run from corresponding subdirectories
for k in `seq 1 ${nProc}`
do
  rm -rf tmcmc_process_${k}
  mkdir tmcmc_process_${k}
  cd tmcmc_process_${k}
  cp ../model.x .
  cp ../mcmcstates_${k}.dat mcmcstates_local.dat
  ./model.x &
  cd ..
done

rm -f tmcmc_ll.dat
list=""
for k in `seq 1 ${nProc}`
do
	while [ ! -f tmcmc_process_${k}/done.txt ]
	do
		sleep ${SLP1}
	done
	list="$list tmcmc_process_${k}/tmcmc_ll.dat"
done

cat $list > tmcmc_ll.dat
