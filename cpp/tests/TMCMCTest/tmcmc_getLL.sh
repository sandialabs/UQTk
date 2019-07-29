#!/bin/bash

nProc=$1
SLP1=0.05    # pause between checking existance of "done" files

if [ ! -f "delta.dat" ]
then
  touch delta.dat
  echo 2 > delta.dat
fi

for k in `seq 1 ${nProc}`
do
  # echo "Likelihood computation: Launching process ${k} of ${nProc}"
  rm -rf run_${k}
  mkdir run_${k}
  cd run_${k}
  cp ../modelLLik.x .
  cp ../mcmcstates_${k}.dat mcmcstates_local.dat
  ./modelLLik.x & 
  cd ..
done

rm -f tmcmc_ll.dat
list=""
for k in `seq 1 ${nProc}`
do
	while [ ! -f run_${k}/done.txt ]
	do
		sleep ${SLP1}
	done
	list="$list run_${k}/tmcmc_ll.dat"
done

cat $list > tmcmc_ll.dat