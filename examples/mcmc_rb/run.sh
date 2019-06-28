#!/bin/bash

for i in `seq 1 10`; do
  iseed=$((2014+$i))
  echo "Running 10k with seed $iseed"
  ./mcmc_rb.x -p 1 -n 10000 -s $iseed >& log1
  mv chain.dat chain_10k_$i.dat
  mkdir $iseed 
  mv samples.dat.* $iseed
  mv loglik.dat.* $iseed
  mv loglk.dat $iseed
  mv tmcmc_ll.dat $iseed
  mv logpr.dat $iseed
  mv logprior.dat.* $iseed
done

for i in `seq 1 10`; do
  iseed=$((2014+$i))
  echo "Running 100k with seed $iseed"
  ./mcmc_rb.x -p 1 -n 100000 -s $iseed >& log1
  mv chain.dat chain_100k_$i.dat
  /bin/rm -rf samples.dat* loglik.dat.* loglk.dat tmcmc_ll.dat
done

for i in `seq 1 10`; do
  iseed=$((2014+$i))
  echo "Running 1000k with seed $iseed"
  ./mcmc_rb.x -p 1 -n 1000000 -s $iseed >& log1
  mv chain.dat chain_1000k_$i.dat
  /bin/rm -rf samples.dat* loglik.dat.* loglk.dat tmcmc_ll.dat
done
