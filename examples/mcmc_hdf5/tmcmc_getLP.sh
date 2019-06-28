#!/bin/bash

iter=$1
lpfile="tmcmc_lp.dat";

tmcmcIter=$(head -n 1 TMCMCiter.dat)

if [ $tmcmcIter -eq 0 ]
then
  ./modelLP.x -f setup.dat -s mcmcstates_1.dat;
else
  curdelta=$(head -n 1 delta.dat)
  deltarate=$(head -n 1 deltaRate.dat)
  echo $(awk "BEGIN {print $curdelta / $deltarate}") > prevdelta.dat
  ./modelLL.x $(head -n 1 prevdelta.dat)
  # python3 PyKDE.py samples.dat.0 samples.dat.${iter};
fi
