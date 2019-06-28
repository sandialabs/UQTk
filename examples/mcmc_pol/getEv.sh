#!/bin/bash

for i in `seq 2 7`; do
  grep wMean log_d${i} | awk -v j=$i 'BEGIN {prod=1} {prod=prod*$6} END {print j-1" "log(prod)}'
done

#for i in `seq 2 7`; do
#  grep wMean log_d${i}_500k | awk -v j=$i 'BEGIN {prod=1} {prod=prod*$6} END {print j-1" "log(prod)}'
#done

