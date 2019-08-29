#!/bin/bash
# Parallel evaluation of prior

nProc=$1
SLP1=0.05    # pause between checking existance of "done" files

for k in `seq 1 ${nProc}`
do
  # echo "Prior computation: Launching process ${k} of ${nProc}"
  rm -rf tmcmc_process_${k}
  mkdir tmcmc_process_${k}
  cd tmcmc_process_${k}
  cp ../bimodal.x .
  cp ../mcmcstates_${k}.dat mcmcstates_local.dat
  ./bimodal.x -p &
  cd ..
done

rm -f tmcmc_lp.dat
list=""
for k in `seq 1 ${nProc}`
do
	while [ ! -f tmcmc_process_${k}/done.txt ]
	do
		sleep ${SLP1}
	done
	list="$list tmcmc_process_${k}/tmcmc_lp.dat"
done

cat $list > tmcmc_lp.dat