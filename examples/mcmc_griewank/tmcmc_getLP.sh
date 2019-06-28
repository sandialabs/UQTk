#!/bin/bash

iter=$1
lpfile="tmcmc_lp.dat";

# initialize nruns
touch ${lpfile}; /bin/rm -rf ${lpfile};
./calcLPGriewank.x -f setup.dat -s mcmcstates_1.dat;
