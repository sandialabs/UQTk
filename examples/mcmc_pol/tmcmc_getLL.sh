#!/bin/bash

nruns=$1
llfile="tmcmc_ll.dat"

# initialize nruns
touch ${llfile}; /bin/rm -rf ${llfile}; touch ${llfile}
python model.py
