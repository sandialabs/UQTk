#!/bin/bash

fact=( 1)

for i in "${fact[@]}"; do
  iseed=$((2015+$i))
  nproc=1
  samples=25000
  iterations=1
  betaThres=0.5
  mala=0
  MFactor=${i}
  dirName="RunResults"
  echo "TMCMC - $samples samples, $iterations stages, $iseed seed"
  echo "MALA - $mala, betaThreshold - $betaThres, MFactor - ${MFactor}"
  mkdir $dirName
  ./mcmc_hdf5.x -M ${MFactor} -p $nproc -B $betaThres -m $mala -n $samples -s $iseed -i $iterations >& log1

  mv RunResults $dirName/.
  rm TMCMCiter.dat
  rm samples.dat.0
  mv log1 $dirName/log1
  rm prevdelta.dat
  mv deltaRate.dat $dirName/deltaRate.dat
  cp rngs.dat $dirName/.
  cp setup.dat $dirName/.
  mv gamma.dat  $dirName/gamma.dat
  mv Evidence.dat $dirName/evidence.dat
  
done
