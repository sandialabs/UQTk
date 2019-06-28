#!/bin/bash

for i in `seq 1 1`; do
  iseed=$((2014+$i))
  nproc=1
  samples=5000
  iterations=5
 
  echo "TMCMC - $samples samples, $iterations stages, $iseed seed" 
 
  mkdir RunResults
  ./mcmc_manifold.x -p $nproc -n $samples -s $iseed -i $iterations >& log1
  rm TMCMCiter.dat
  rm samples.dat.0
  mv log1 RunResults/log1
  rm prevdelta.dat
  mv deltaRate.dat RunResults/deltaRate.dat
  cp rngs.dat RunResults/.
  cp setup.dat RunResults/.
  mv gamma.dat  RunResults/gamma.dat
  mv Evidence.dat RunResults/evidence.dat
  
done

