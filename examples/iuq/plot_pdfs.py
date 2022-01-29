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


import os
import argparse
import numpy as np
import matplotlib.pyplot as plt

from PyUQTk.plotting.pdfs import plot_pdf1d, plot_pdf2d


plt.rc('legend', loc='best', fontsize=12)
plt.rc('lines', linewidth=4, color='r')
plt.rc('axes', linewidth=2, grid=True, labelsize=22)


usage_str = 'Script to plot PDFs given samples, either triangular or individual.'
parser = argparse.ArgumentParser(description=usage_str)
parser.add_argument('ind_show', type=int, nargs='*',
                    help="indices of requested parameters (count from 0)")
parser.add_argument("-p", "--samples_file", dest="samples_file", type=str, default='pchain.dat',
                    help="Samples file")
parser.add_argument("-n", "--names_file", dest="names_file", type=str, default=None,
                    help="Names file")
parser.add_argument("-l", "--nominal_file", dest="nominal_file", type=str, default=None,
                    help="Nominals file")
parser.add_argument("-g", "--prangeshow_file", dest="prangeshow_file", type=str, default=None,
                    help="Prior range file")
parser.add_argument("-t", "--plot_type", dest="plot_type", type=str, default='tri',
                    help="Plot type", choices=['tri', 'ind'])
parser.add_argument("-b", "--burnin", dest="burnin", type=int, default=0,
                    help="Samples burnin")
parser.add_argument("-e", "--every", dest="every", type=int, default=1,
                    help="Samples thinning")
parser.add_argument("-s", "--nsamxi", type=int, dest="nsam_xi", default=1,
                    help="Number of xi samples per posterior sample")
parser.add_argument("-x", "--lsize", type=int, dest="lsize", default=10,
                    help="Label size")
args = parser.parse_args()

# check if ind_show indeed counts from 1
# handle axes labels' size
# handle legends for ind plot_type
# add verbosity -v
# formatting ticklabels
# subplot or add_axes afterall?

ind_show = args.ind_show
samples_file = args.samples_file
plot_type = args.plot_type
names_file = args.names_file
burnin = args.burnin
every = args.every
nominal_file = args.nominal_file
prangeshow_file = args.prangeshow_file
nsam_xi = args.nsam_xi
lsize = args.lsize

plt.rc('axes', labelsize=lsize)
plt.rc('xtick', labelsize=lsize)
plt.rc('ytick', labelsize=lsize)

if (nsam_xi > 1):
    assert(burnin == 0)
    assert(every == 1)

samples_all = np.loadtxt(samples_file, ndmin=2)
if samples_all.shape[0] == 1:
    plot_type = 'ind'  # add message!!
    samples_all = samples_all.T
samples_all = samples_all[burnin::every, :]
npar_all = samples_all.shape[1]

if len(ind_show) == 0:
    ind_show = range(npar_all)
if len(ind_show) == 1:
    plot_type = 'ind'  # add message!!

samples = samples_all[:, ind_show]
npar = len(ind_show)

if names_file is not None and os.path.exists(names_file):
    with open(names_file) as f:
        names = f.read().splitlines()
        assert(len(names) == npar_all)
else:
    names = ['Param ' + str(i) for i in range(npar_all)]

if nominal_file is not None:
    nominals = np.loadtxt(nominal_file, ndmin=1)

show_range = False
if prangeshow_file is not None:
    prange = np.loadtxt(prangeshow_file, ndmin=2)
    show_range = True

if plot_type == 'tri':
    figs, axarr = plt.subplots(npar, npar, sharex='col', figsize=(15, 15))
else:
    axarr = []
    figs = []

for i in range(npar):
    print("Quantity # ", ind_show[i], "; Name : ", names[ind_show[i]])
    if plot_type == 'tri':
        thisax = axarr[i, i]
    else:
        fig = plt.figure(figsize=(8, 8))
        thisax = plt.gca()
        axarr.append(thisax)
        figs.append(fig)

    for isam_xi in range(nsam_xi):
        plot_pdf1d(samples[isam_xi::nsam_xi, i], pltype='kde', ax=thisax)
        thisax.set_ylim(bottom=0)
    # plot_pdf1d(samples[:, i], pltype='sam', ax=thisax)

    if nominal_file is not None:
        plot_pdf1d(np.array(nominals[ind_show[i]], ndmin=2), pltype='nom',
                   nom_height_factor=1., color='r', ax=thisax)
    if show_range:
        dr = prange[ind_show[i], 1] - prange[ind_show[i], 0]
        thisax.plot([prange[ind_show[i], 0], prange[ind_show[i], 0],
                     prange[ind_show[i], 1], prange[ind_show[i], 1]],
                    [0.0, 1.0 / dr, 1.0 / dr, 0.0], 'g--', linewidth=0.4, label='Prior', zorder=-100)

    x0, x1 = thisax.get_xlim()
    y0, y1 = thisax.get_ylim()
    # thisax.set_aspect((x1 - x0) / (y1 - y0))
    thisax.set_title('PDF of ' + names[ind_show[i]], fontsize=lsize)

    if plot_type == 'tri':
        if i == 0:
            thisax.set_ylabel(names[ind_show[i]])
        if i == npar - 1:
            thisax.set_xlabel(names[ind_show[i]])
        if i > 0:
            thisax.yaxis.set_ticks_position("right")
        # thisax.yaxis.set_label_coords(-0.12, 0.5)
    else:
        thisax.set_xlabel(names[ind_show[i]])
        #plt.savefig('pdf_' + str(ind_show[i]) + '.eps')
        plt.savefig('pdf_' + names[ind_show[i]] + '.eps')
        # plt.clf()

    for j in range(i):
        if plot_type == 'tri':
            thisax = axarr[i, j]
            axarr[j, i].axis('off')
        else:
            plt.figure(figsize=(8, 8))
            thisax = plt.gca()

        for isam_xi in range(nsam_xi):
            plot_pdf2d(samples[isam_xi::nsam_xi, j], samples[isam_xi::nsam_xi, i],
                       pltype='kde', ncont=10, ax=thisax)
        #plot_pdf2d(samples[:, j], samples[:, i], pltype='sam', mstyle='x', ax=thisax)

        if nominal_file is not None:
            plot_pdf2d(nominals[ind_show[j]], nominals[ind_show[i]],
                       pltype='sam', color='r', mstyle='x', ax=thisax)

        x0, x1 = thisax.get_xlim()
        y0, y1 = thisax.get_ylim()
        #thisax.set_aspect((x1 - x0) / (y1 - y0))

        if plot_type == 'tri':
            if j == 0:
                thisax.set_ylabel(names[ind_show[i]])
            if i == npar - 1:
                thisax.set_xlabel(names[ind_show[j]])
            if j > 0:
                thisax.yaxis.set_ticklabels([])

        else:
            thisax.set_ylabel(names[ind_show[i]])
            thisax.set_xlabel(names[ind_show[j]])
            #plt.savefig('pdf_' + str(ind_show[j]) + '_' + str(ind_show[i]) + '.eps')
            plt.savefig('pdf_' + names[ind_show[j]] + '_' + names[ind_show[i]] + '.eps')
            # plt.clf()


if plot_type == 'tri':
    plt.savefig('pdf_tri.eps')
    # plt.show()
