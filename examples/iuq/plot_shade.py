#!/usr/bin/env python
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
#=====================================================================================

import argparse
import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt

import PyUQTk.plotting.fits as ft


plt.rc('legend', loc='best', fontsize=25)
plt.rc('lines', linewidth=4, color='r')
plt.rc('axes', linewidth=3, grid=True, labelsize=32)
plt.rc('xtick', labelsize=22)
plt.rc('ytick', labelsize=22)

##########################################################################
##########################################################################

# custom_xlabel = ''
# custom_ylabel = ''
# custom_ylim = [-1, 10]
# custom_title = args.title

# Parse input arguments
usage_str = 'Script to plot shaded PDF fits to a given slice of outputs.'
parser = argparse.ArgumentParser(description=usage_str)
parser.add_argument("-i", "--ind_plot", dest="indplot_file", type=str,
                    default=None, help="File for indices of plotted outputs")
parser.add_argument("-x", "--xdata", dest="xdata_file",
                    type=str, default=None, help="Xdata file")
parser.add_argument("-y", "--prior_out", dest="prior_output_file",
                    type=str, default=None, help="Prior output file")
parser.add_argument("-z", "--post_out", dest="post_output_file",
                    type=str, default=None, help="Posterior output file")
parser.add_argument("-d", "--ydata", dest="data_file",
                    type=str, default=None, help="Data file")
parser.add_argument("-s", "--ydata_std", dest="datastd_file",
                    type=str, default=None, help="Data st. deviation file")
parser.add_argument("-c", "--ncol", dest="ncol",
                    type=int, default=1, help="The relevant column of xdata file (count from 0)")
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
prior_output_file = args.prior_output_file
post_output_file = args.post_output_file
data_file = args.data_file
datastd_file = args.datastd_file
xdata_file = args.xdata_file
ncol = args.ncol
custom_title = args.title
custom_xlabel = args.xlabel
custom_ylabel = args.ylabel

if args.xticklabels_file is not None:
    with open(args.xticklabels_file) as f:
        custom_xticklabels = f.read().splitlines()

    #custom_xticklabels = np.loadtxt(args.xticklabels, dtype=str) # ideally read

prior_out_flag = False
if prior_output_file is not None:
    prior_output = np.loadtxt(prior_output_file)
    prior_out_flag = True
    nout = prior_output.shape[1]

post_out_flag = False
if post_output_file is not None:
    post_output = np.loadtxt(post_output_file)
    post_out_flag = True
    nout = post_output.shape[1]

assert(prior_out_flag or post_out_flag)
print(prior_output.shape)
print(post_output.shape)
assert(prior_output.shape[1] == post_output.shape[1])

if indplot_file is not None:
    ind_plot = np.loadtxt(indplot_file, dtype=int)
else:
    ind_plot = np.arange(nout)

if xdata_file is not None:
    xdata = np.loadtxt(xdata_file)[ind_plot, ncol]
else:
    xdata = np.arange(nout)

nout_plot = ind_plot.shape[0]

fig = plt.figure(figsize=(12, 9))
thisax = plt.gca()

if prior_out_flag:
    ft.plot_shade(thisax, xdata,
                  prior_output[:, ind_plot].T, cmap=cm.OrRd, grid_show=True)
if post_out_flag:
    ft.plot_shade(thisax, xdata,
                  post_output[:, ind_plot].T, grid_show=True)


thisax.set_xticks(xdata)

if data_file is not None:
    bcg_data = np.loadtxt(data_file, ndmin=2)
    for j in range(bcg_data.shape[1]):
        if datastd_file is not None:
            bcg_data_std = np.loadtxt(datastd_file, ndmin=2)
            thisax.errorbar(xdata, bcg_data[:, j], yerr=bcg_data_std[:, j],
                            ecolor='k', fmt='ko', label='Data', ms=12, zorder=100000)
        else:
            thisax.plot(xdata, bcg_data[:, j], 'ko', label='Data', ms=12, zorder=100000)

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

# thisax.text(1.1, 10, 'Prior', color='r', fontsize=33)
# thisax.text(1.1, 9, 'Posterior', color='g', fontsize=33)
handles, labels = thisax.get_legend_handles_labels()
handles.append(plt.Rectangle((0, 0), 1, 1, fc='r'))
labels.append('Prior')
handles.append(plt.Rectangle((0, 0), 1, 1, fc='g'))
labels.append('Posterior')
thisax.legend(handles, labels, fontsize=24, ncol=1)

# f.tight_layout()

plt.savefig('fit_shade.eps')
# show()
