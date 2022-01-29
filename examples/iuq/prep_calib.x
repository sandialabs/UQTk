#!/bin/bash -e
#=====================================================================================
#
#                      The UQ Toolkit (UQTk) version 3.1.2
#                          Copyright (2022) NTESS
#                        https://www.sandia.gov/UQToolkit/
#                        https://github.com/sandialabs/UQTk
#
#     Copyright 2022 National Technology & Engineering Solutions of Sandia, LLC (NTESS).
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

PMTYPE=$1 # Usually pcs
NOUT=$2 # Needed only if PMTYPE=pcs

UQPC=${UQTK_INS}/examples/uqpc


# Given the indices of calibrating data, extract the appropriate sub-matrixes of data and model surrogate
paste ind_calib.dat xdata_all.txt     | awk '$1==1{for (i=2; i<=NF; i++) printf("%lg ",$i); printf("\n")}' > xdata_calib.txt
paste ind_calib.dat ydata_all.txt     | awk '$1==1{for (i=2; i<=NF; i++) printf("%lg ",$i); printf("\n")}' > ydata_calib.txt
paste ind_calib.dat dataerr_all.dat   | awk '$1==1{for (i=2; i<=NF; i++) printf("%lg ",$i); printf("\n")}' > dataerr_calib.dat
${UQPC}/transpose_file.x ind_calib.dat > ind_select.tmp
if [[ $PMTYPE == "pc" ]]; then
    cp pccf_all_all.dat pccf_all_pred.dat
    cat ind_select.tmp pccf_all_all.dat | awk '{if (NR==1) for (i=1; i<=NF; i++) indsel[i]=$i}; {if (NR>1) for (i=1; i<=NF; i++) if (indsel[i]==1) printf("%lg ",$i); if (NR>1) printf("\n")}' > pccf_all.dat
elif [[ $PMTYPE == "pcs" ]]; then
    for (( i=0; i<$NOUT ; i++ )); do
        cp pccfp.${i}_all.dat pccfp.${i}_pred.dat
        cp mindexp.${i}_all.dat mindexp.${i}_pred.dat
    done
    j=0
    i=0
    while read  p; do
        if [[ $p == "1" ]]; then
            cp pccfp.${i}_all.dat pccfp.${j}.dat
            cp mindexp.${i}_all.dat mindexp.${j}.dat
            j=`expr $j + 1`
        fi
        i=`expr $i + 1`
    done < ind_calib.dat

fi
