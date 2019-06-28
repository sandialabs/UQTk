#!/bin/bash

folderName=$1
lastSample=$2

mkdir $folderName
mv samples.dat.* $folderName
cp $folderName/$lastSample samples.dat.0
mv loglik.dat.* $folderName
mv loglk.dat $folderName
mv tmcmc_ll.dat $folderName
mv tmcmc_lp.dat $folderName
mv logpr.dat $folderName
mv logprior.dat.* $folderName
mv mcmcstates_* $folderName
mv chain* $folderName
mv delta.dat $folderName
mv $folderName RunResults
