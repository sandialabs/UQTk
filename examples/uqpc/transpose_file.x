#!/bin/sh -e
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

# Transpose a matrix: assumes all lines have same number
# of fields

# Example: transpose_file.x file_in > file_out

shopt -s expand_aliases
alias awk="awk -v OFMT='%.15e'"

exec awk '
NR == 1 {
	n = NF
	for (i = 1; i <= NF; i++)
		row[i] = $i
	next
}
{
	if (NF > n)
		n = NF
	for (i = 1; i <= NF; i++)
		row[i] = row[i] " " $i
}
END {
	for (i = 1; i <= n; i++)
		print row[i]
}' ${1+"$@"}
