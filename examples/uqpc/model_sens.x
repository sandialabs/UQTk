#!/bin/bash
#=====================================================================================
#
#                      The UQ Toolkit (UQTk) version 3.1.3
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


# Sensitivity evaluation black-box for given PC coefficients
# model_sens.x <pc_coef_file> <sensitivity_file>
# pc_coef_file     : input, size SxK, K is pc basis size, and S is number of samples
# sensitivity_file : output, size dx(2S), where d is dimensionality,
#                            and columns alternate main and total sensitivities for all S samples

rm -rf ms ts
${UQTK_INS}/examples/uqpc/transpose_file.x $1 > pccfs.dat


NSAM=`echo | awk 'END{print NR}' $1`

for (( i=1; i<=$NSAM ; i++ )); do

	awk '{print $i}' i=$i pccfs.dat > pcf_this.dat
	${UQTK_INS}/bin/pce_sens -x LU -m mi.dat -f pcf_this.dat > pcsens.log
	${UQTK_INS}/examples/uqpc/transpose_file.x mainsens.dat >> ms
	${UQTK_INS}/examples/uqpc/transpose_file.x totsens.dat >> ts
done

paste ms ts > $2
