#!/bin/bash
#=====================================================================================
#
#                      The UQ Toolkit (UQTk) version 3.1.5
#                          Copyright (2024) NTESS
#                        https://www.sandia.gov/UQToolkit/
#                        https://github.com/sandialabs/UQTk
#
#     Copyright 2024 National Technology & Engineering Solutions of Sandia, LLC (NTESS).
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
#     Questions? Contact the UQTk Developers at https://github.com/sandialabs/UQTk/discussions
#     Sandia National Laboratories, Livermore, CA, USA
#=====================================================================================

# python
/bin/rm -rf *.pyc

# pdf files
/bin/rm -rf forUQ*.pdf

# log files
/bin/rm -rf log*.dat

# intermediate files
/bin/rm -rf jointsens.dat mainsens.dat totsens.dat varfrac.dat
/bin/rm -rf PCEspls_100k.dat output_val.dat output_val_pc.dat
/bin/rm -rf solution.dat sp_mindex.1.dat PCcoeff_quad.dat mipc.dat
/bin/rm -rf surf_rxn.in.parsed.xml surf_rxn.in.xml
/bin/rm -rf input.dat output.dat input_fcn.dat output_fcn.dat input_val.dat
/bin/rm -rf xy*.dat
