#!/bin/bash -e
#=====================================================================================
#
#                      The UQ Toolkit (UQTk) version 3.1.4
#                          Copyright (2023) NTESS
#                        https://www.sandia.gov/UQToolkit/
#                        https://github.com/sandialabs/UQTk
#
#     Copyright 2023 National Technology & Engineering Solutions of Sandia, LLC (NTESS).
#     Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government
#     retains certain rights in this software.
#
#     This file is part of The UQ Toolkit (UQTk)
#
#     UQTk is open source software: you can redistribute it and/or modify
#     it under the terms of BSD 3-Clause License
#
#     UQTk is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     BSD 3 Clause License for more details.
#
#     You should have received a copy of the BSD 3 Clause License
#     along with UQTk. If not, see https://choosealicense.com/licenses/bsd-3-clause/.
#
#     Questions? Contact the UQTk Developers at <uqtk-developers@software.sandia.gov>
#     Sandia National Laboratories, Livermore, CA, USA
#=====================================================================================
#=====================================================================================

# Script-example for running model inference
# User should play with this and not use as-is

UQBIN=${UQTK_INS}/bin

## Inference script
## Run ${UQBIN}/model_inf -h for the full set of arguments

PMTYPE=$1 #pcs
#DIMFIX=`awk 'END{print NR}' fixindnom.dat`
DIMFIX=0
DIMALL=`awk 'END{print NR}' prange.dat`
DIM=`expr $DIMALL - $DIMFIX`
PCORD=2 # This is not the model surrogate order, this is PC order for embedded NISP, should be 0 for classical inference
NMCMC=5000 #500000
GAMMA=0.05
EE=dataerr_calib.dat #-0.1 #dataerr_calib.dat # negative means infer the data noise

## Classical calibration
#${UQBIN}/model_inf -l classical -o 0 -i uniform -t xdata_all.txt -f $PMTYPE -c LU -s pci -d $DIM -z -x xdata_calib.txt -y ydata_calib.txt   -m $NMCMC -g $GAMMA -u 5 -e $EE #-a -1 -b 1 # -e $EE -j chainstart.dat

printf "0\n2" > de_params.dat # embedding in the first two parameters
## Embedded model error
${UQBIN}/model_inf -l gausmarg -r de_params.dat -o $PCORD -i uniform_LUpci -t xdata_all.txt -f $PMTYPE -c LU -s pci -d $DIM -z -x xdata_calib.txt -y ydata_calib.txt -e $EE  -m $NMCMC -g $GAMMA -u 5 -a -1 -b 1 #-j chainstart.dat
