#!/bin/bash

bfac=0.1
npts=500
nchn=10

spls=("10k" "100k" "1000k")

for j in "${spls[@]}"; do
  for i in `seq 1 ${nchn}`; do
    ../../bin/pdf_cl -i chain_${j}_${i}.dat -g ${npts} -f ${bfac}
    mv dens.dat dens_${j}_${i}.dat 
  done
done

