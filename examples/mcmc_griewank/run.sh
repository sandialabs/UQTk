#!/bin/bash

for i in `seq 1 1`; do
  iseed=$((2014+$i))
  nproc=1
  nsamples=25000 
  spc=0
  
  echo "Running TMCMC with $nsamples samples, $nproc processors, $iseed seed"
  ./mcmc_griewank.x -p $nproc -n $nsamples -s $iseed -c $spc >& log1
  mv chain.dat chain_10k_$i.dat
  mkdir $iseed 
  mv samples.dat.* $iseed
  mv loglik.dat.* $iseed
  mv loglk.dat $iseed
  mv tmcmc_ll.dat $iseed
  mv tmcmc_lp.dat $iseed
  mv logpr.dat $iseed
  mv logprior.dat.* $iseed
  mv log1 $iseed
  mv gamma.dat $iseed
done

