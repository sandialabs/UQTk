#!/usr/bin/env python
#=====================================================================================
#
#                      The UQ Toolkit (UQTk) version @UQTKVERSION@
#                          Copyright (@UQTKYEAR@) NTESS
#                        https://www.sandia.gov/UQToolkit/
#                        https://github.com/sandialabs/UQTk
#
#     Copyright @UQTKYEAR@ National Technology & Engineering Solutions of Sandia, LLC (NTESS).
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
#=====================================================================================

import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt

import PyUQTk.plotting.fits as ft


plt.rc('legend', loc='best', fontsize=25)
plt.rc('lines', linewidth=4, color='r')
plt.rc('axes', linewidth=3, grid=True, labelsize=32)
plt.rc('xtick', labelsize=22)
plt.rc('ytick', labelsize=22)

##########################################################################
##########################################################################

varlabels = ['Surrogate error', 'Posterior uncertainty', 'Model error']
varcolors = ['grey', 'blue', 'lightblue']
# varlabels = ['Surrogate error', 'Posterior uncertainty']
# varcolors = ['grey', 'blue']
# custom_ylim = [-1, 10]


# Parse input arguments
usage_str = 'Script to plot 1d fit slices with variance decomposition.'
parser = argparse.ArgumentParser(description=usage_str)
parser.add_argument("-i", "--ind_plot", dest="indplot_file", type=str,
                    default=None, help="File for indices of plotted outputs")
parser.add_argument("-x", "--xdata", dest="xdata_file",
                    type=str, default=None, help="Xdata file")
parser.add_argument("-y", "--samples", dest="samples_file",
                    type=str, default=None, help="Samples file")
parser.add_argument("-m", "--mean_file", dest="mean_file",
                    type=str, default=None, help="Output mean file")
parser.add_argument("-v", "--var_file", dest="var_file",
                    type=str, default=None, help="Output variance file")
parser.add_argument("-d", "--ydata", dest="data_file",
                    type=str, default=None, help="Data file")
parser.add_argument("-s", "--ydata_std", dest="datastd_file",
                    type=str, default=None, help="Data st. deviation file")
parser.add_argument("-c", "--ncol", dest="ncol",
                    type=int, default=1, help="The relevant column of xdata file (count from 0)")
parser.add_argument("-r", "--interp", dest="interp",
                    type=str, default=None, help="Interpolation type",
                    choices=[None, 'linear', 'cubic'])
parser.add_argument("-t", "--title", dest="title",
                    type=str, default='', help="Title")
parser.add_argument("-u", "--xlabel", dest="xlabel",
                    type=str, default='', help="X label")
parser.add_argument("-w", "--ylabel", dest="ylabel",
                    type=str, default='', help="Y label")
parser.add_argument("-l", "--xticklabels", dest="xticklabels_file",
                    type=str, default=None, help="Xtick labels file")

args = parser.parse_args()

indplot_file = args.indplot_file
mean_file = args.mean_file
var_file = args.var_file
data_file = args.data_file
datastd_file = args.datastd_file
xdata_file = args.xdata_file
interp = args.interp
ncol = args.ncol
samples_file = args.samples_file

custom_title = args.title
custom_xlabel = args.xlabel
custom_ylabel = args.ylabel

if args.xticklabels_file is not None:
    with open(args.xticklabels_file) as f:
        custom_xticklabels = f.read().splitlines()

if mean_file is not None:
    mean = np.loadtxt(mean_file)
    nout = mean.shape[0]
else:
    print("Mean file needs to be given, with flag -m. Exiting now.")
    sys.exit()

if var_file is not None:
    var = np.loadtxt(var_file)
    if len(var.shape) == 1:
        var = var.reshape(-1, 1)
else:
    var = np.zeros((nout, 1))

if indplot_file is not None:
    ind_plot = np.loadtxt(indplot_file, dtype=int)
else:
    ind_plot = np.arange(nout)


if samples_file is not None:
    samples = np.loadtxt(samples_file)
    ysam = samples[:, ind_plot].T
else:
    ysam = None

if xdata_file is not None:
    xdata_all = np.loadtxt(xdata_file)
    if len(xdata_all.shape) == 1:
        assert(ncol == 0)
        xdata = xdata_all[ind_plot]
    else:
        xdata = xdata_all[ind_plot, ncol]
else:
    xdata = np.arange(nout)

nout_plot = ind_plot.shape[0]

fig = plt.figure(figsize=(12, 9))
thisax = plt.gca()

ft.plot_vars(thisax, xdata, mean[ind_plot], variances=var[ind_plot, :], ysam=ysam,
             stdfactor=1., varlabels=varlabels, varcolors=varcolors,
             interp=interp, connected=False)

thisax.set_xticks(xdata)

if data_file is not None:
    bcg_data = np.loadtxt(data_file, ndmin=2)
    for j in range(bcg_data.shape[1]):
        if datastd_file is not None:
            bcg_data_std = np.loadtxt(datastd_file, ndmin=2)
            thisax.errorbar(xdata, bcg_data[:, j], yerr=bcg_data_std[:, j],
                            ecolor='k', fmt='ko', label='Data', ms=12, zorder=50000)
        else:
            thisax.plot(xdata, bcg_data[:, j], 'ko', label='Data', ms=12, zorder=50000)


try:
    thisax.set_xticklabels(custom_xticklabels)
except NameError:
    pass


try:
    thisax.set_ylim(custom_ylim)
except NameError:
    pass

thisax.set_xlabel(custom_xlabel)
thisax.set_ylabel(custom_ylabel)
thisax.set_title(custom_title, size=32)


handles, labels = thisax.get_legend_handles_labels()
if ysam is not None:
    handles.append(plt.Line2D((0, 1), (0, 0), color='r'))
    labels.append(r'Prior ensemble')

# if pred_show:
#     if moderr:
#         handles.append(Rectangle((0, 0), 1, 1, fc=color_moderr))
#         labels.append(r'Model error')
#     if surrerr_show:
#         handles.append(Rectangle((0, 0), 1, 1, fc=color_surrerr))
#         labels.append(r'Surrogate error')
#     handles.append(Rectangle((0, 0), 1, 1, fc=color_dataerr))
#     labels.append(r'Posterior uncertainty')

thisax.legend(handles, labels, fontsize=14, ncol=1)
#thisax.grid(False)
plt.savefig('fit_vars.eps')
# show()
