#!/bin/bash

nruns=$1
llfile="tmcmc_ll.dat"

# initialize nruns
touch ${llfile}; /bin/rm -rf ${llfile}; touch ${llfile}
for i in `seq 1 $nruns`; do
  awk '{print -(1-$1)**2-100.*($2-$1*$1)**2}' mcmcStates_$i.dat >> ${llfile}
done

