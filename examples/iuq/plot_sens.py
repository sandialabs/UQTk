#!/usr/bin/env python
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


# Import necessary libraries
import os
import sys
import argparse
import numpy as np


# if (sys.version_info.major == 2):
#     import cPickle as pick
# elif (sys.version_info.major == 3):
#     import pickle as pick
# else:
#     print("Only Python 2 or 3 are supported. Exiting.")
#     sys.exit()


import PyUQTk.plotting.surrogate as ss


# Parse input arguments
usage_str = 'Script to plot sensitivity bars given indices of plotted outputs.'
parser = argparse.ArgumentParser(description=usage_str)
parser.add_argument("-i", "--ind_plot", dest="indplot_file", type=str,
                    default=None, help="File for indices of plotted outputs")
parser.add_argument("-s", "--sens", dest="sens_file",
                    type=str, default='sens.dat', help="Sensitivity file")
parser.add_argument("-p", "--pnames", dest="pnames_file",
                    type=str, default=None, help="Parameter names file")
parser.add_argument("-o", "--outnames", dest="outnames_file",
                    type=str, default=None, help="Output names file")
parser.add_argument("-c", "--colors", dest="colors_file",
                    type=str, default=None, help="Colors file")
args = parser.parse_args()

indplot_file = args.indplot_file
sens_file = args.sens_file
outnames_file = args.outnames_file
pnames_file = args.pnames_file
colors_file = args.colors_file

if sens_file is not None:
    sens = np.loadtxt(sens_file)
    nout = sens.shape[0]
    npar = sens.shape[1]
else:
    print("Sensitivity file needs to be given, with flag -s. Exiting now.")
    sys.exit()

if indplot_file is not None:
    ind_plot = np.loadtxt(indplot_file, dtype=int)
else:
    ind_plot = np.arange(nout)


if outnames_file is not None:
    with open(outnames_file) as f:
        outnames = f.read().splitlines()
else:
    outnames = [str(i) for i in range(1, nout + 1)]

if pnames_file is not None:
    with open(pnames_file) as f:
        pnames = f.read().splitlines()
else:
    pnames = ['Param # ' + str(i) for i in range(1, npar + 1)]

nout_plot = ind_plot.shape[0]

pars = np.arange(npar)
cases = np.arange(nout_plot)
outnames = [outnames[j] for j in ind_plot]

if colors_file is not None:
    colors = np.loadtxt(colors_file)
else:
    colors = []
    # colors = np.random.rand(npar, 3)
    # np.savetxt('colors.txt', colors)


ss.plot_sens(sens[ind_plot], pars, cases, vis="bar", reverse=False, par_labels=pnames,
             case_labels=outnames, ncol=4, legend_show=2, legend_size=26,
             grid_show=False, xlbl='', topsens=8, colors=colors,
             yoffset=-0.05, showplot=False)
