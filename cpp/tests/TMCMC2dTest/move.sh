#!/bin/bash

folderName=$1
lastSample=$2

mkdir $folderName
mv samples.dat.* $folderName
cp $folderName/$lastSample samples.dat.0
mv loglik.dat.* $folderName
mv tmcmc_ll.dat $folderName
mv mcmcstates_* $folderName
mv delta.dat $folderName
mv $folderName RunResults
