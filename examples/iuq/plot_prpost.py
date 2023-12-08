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

import os
import argparse
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

import PyUQTk.plotting.pdfs as uqp

plt.rc('legend', loc='best', fontsize=25)
plt.rc('lines', linewidth=4, color='r')
plt.rc('axes', linewidth=3, grid=True, labelsize=25)
plt.rc('xtick', labelsize=18)
plt.rc('ytick', labelsize=18)

##########################################################################
##########################################################################


# Parse input arguments
usage_str = 'Script to plot PDFs given samples, either triangular or individual.'
parser = argparse.ArgumentParser(description=usage_str)
parser.add_argument('ind_show', type=int, nargs='*',
                    help="indices of requested parameters")
parser.add_argument("-y", "--prior_out", dest="prior_output_file",
                    type=str, default=None, help="Prior output file")
parser.add_argument("-z", "--post_out", dest="post_output_file",
                    type=str, default=None, help="Posterior output file")
parser.add_argument("-d", "--ydata", dest="data_file",
                    type=str, default=None, help="Data file")
parser.add_argument("-n", "--names_file", dest="names_file", type=str, default=None,
                    help="Names file")
parser.add_argument("-t", "--title", dest="title",
                    type=str, default='', help="Title")

args = parser.parse_args()

ind_show = args.ind_show
prior_output_file = args.prior_output_file
post_output_file = args.post_output_file
data_file = args.data_file
names_file = args.names_file
custom_title = args.title

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
# assert(prior_output.shape[1] == post_output.shape[1])


if names_file is not None and os.path.exists(names_file):
    with open(names_file) as f:
        names = f.read().splitlines()
        assert(len(names) == nout)
else:
    names = ['Param ' + str(i) for i in range(nout)]

nout_show = len(ind_show)


fig = plt.figure(figsize=(12, 9))
thisax = plt.gca()
###############

for i in range(nout_show):

    fig = plt.figure(figsize=(8, 8))
    thisax = plt.gca()
    leglines = []
    leglabels = []
    if prior_out_flag:
        uqp.plot_pdf1d(prior_output[:, ind_show[i]], pltype='kde', color='r', ax=thisax)
        leglines.append(Line2D([0], [0], color='red', linewidth=2, linestyle='-'))
        leglabels.append('Prior')
    if post_out_flag:
        uqp.plot_pdf1d(post_output[:, ind_show[i]], pltype='kde', color='b', ax=thisax)
        leglines.append(Line2D([0], [0], color='blue', linewidth=2, linestyle='-'))
        leglabels.append('Posterior')

    if data_file is not None:
        bcg_data = np.loadtxt(data_file, ndmin=2)
        # uqp.plot_pdf1d(bcg_data, pltype='hist', color='g', ax=thisax)
        uqp.plot_pdf1d(bcg_data, pltype='nom', color='g', ax=thisax, nom_height_factor=20.0)
        leglines.append(Line2D([0], [0], color='green', linewidth=2, linestyle='--'))
        leglabels.append('Obs. data')
        # if len(bcg_data.shape) == 1:
        #     bcg_data = bcg_data.reshape(-1, 1)
        #     uqp.plot_pdf1d(bcg_data[ind_show[i],:], pltype='nom', color='g', nom_height_factor=1, ax=thisax)
        # else:
        #     uqp.plot_pdf1d(bcg_data[ind_show[i],:], pltype='hist', color='g', ax=thisax)

    thisax.set_ylim(bottom=0)
    thisax.set_xlabel(names[ind_show[i]])
    thisax.set_ylabel('PDF')
    # semi-HARDWIRED!!
    thisax.legend(leglines, leglabels, fontsize=22)

    #plt.savefig('pdf_prpost_' + str(ind_show[i]) + '.eps')
    plt.savefig('pdf_prpost_' + names[ind_show[i]] + '.eps')



try:
    thisax.set_title(custom_title, size=32)
except NameError:
    pass

# text(2, 8, 'Prior', color='r',fontsize=33)
# f.tight_layout()


