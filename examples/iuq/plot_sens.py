#!/usr/bin/env python
#=====================================================================================
#                     The UQ Toolkit (UQTk) version @UQTKVERSION@
#                     Copyright (@UQTKYEAR@) Sandia Corporation
#                     http://www.sandia.gov/UQToolkit/
#
#     Copyright (@UQTKYEAR@) Sandia Corporation. Under the terms of Contract DE-AC04-94AL85000
#     with Sandia Corporation, the U.S. Government retains certain rights in this software.
#
#     This file is part of The UQ Toolkit (UQTk)
#
#     UQTk is free software: you can redistribute it and/or modify
#     it under the terms of the GNU Lesser General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
#
#     UQTk is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU Lesser General Public License for more details.
#
#     You should have received a copy of the GNU Lesser General Public License
#     along with UQTk.  If not, see <http://www.gnu.org/licenses/>.
#
#     Questions? Contact Bert Debusschere <bjdebus@sandia.gov>
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
