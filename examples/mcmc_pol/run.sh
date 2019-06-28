#!/bin/bash

nspl=10000
dst=2
den=7

cmin=-10.0
cmax=10.0

echo "${cmin} ${cmax}" > rngs.dat
echo "n 0 1" > setup.dat
for i in `seq ${dst} ${den}`; do
  echo "${cmin} ${cmax}" >> rngs.dat
  echo "n 0 1" >> setup.dat
  ./mcmc_pol.x -d $i -p 1 -n $nspl >& log_d${i}_test
  mkdir $i
  mv samples.dat.* $i
  mv logl* $i
  mv gamma* $i
  mv Evid* $i
done
