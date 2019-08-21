#!/bin/bash

folderName=$1

# rm -rf $folderName
mkdir $folderName
mv samples.dat.* $folderName
mv loglik.dat.* $folderName
rm tmcmc_ll.dat $folderName
mv logprior.dat.* $folderName
rm tmcmc_lp.dat $folderName
rm mcmcstates_* $folderName
