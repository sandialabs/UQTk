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
#     Questions? Contact the UQTk Developers at <uqtk-developers@software.sandia.gov>
#     Sandia National Laboratories, Livermore, CA, USA
#=====================================================================================
import os
import sys

sys.path.append(os.environ['UQTK_SRC'] + '/PyUQTk')
from PyUQTk.utils import pdf_kde

try:
    import numpy as np
except ImportError:
    print("Numpy was not found.")

import matplotlib.pyplot as plt

#############################################################
#############################################################
#############################################################


def plot_pdf1d(sams, pltype='hist', color='b', lw=1, nom_height_factor=10., ax=None):
    if ax is None:
        ax = plt.gca()

    if pltype == 'kde':
        ngrid = 111
        a = 3.
        pnom = np.mean(sams)
        pdel = np.std(sams)
        pgrid = np.linspace(pnom - a * pdel, pnom + a * pdel, ngrid)
        _, pdf = pdf_kde.get_pdf(sams, pgrid, verbose=0, method='Python')
        ax.plot(pgrid, pdf, color, linewidth=lw, label='')
    elif pltype == 'hist':
        n, bins, patches = ax.hist(sams, bins=20,
                                   density=True,
                                   facecolor=color, alpha=1.)
    elif pltype == 'sam':
        ax.plot(sams, np.zeros_like(sams), 'ro', ms=3, markerfacecolor='None', label='')
    elif pltype == 'nom':
        y1, y2 = ax.get_ylim()
        y2f = nom_height_factor
        for sam in sams:
            ax.plot([sam, sam], [y1, y2 / y2f], '--', color=color, linewidth=lw, label='')

    else:
        print("Plot type is not recognized. Exiting")
        sys.exit()

    #ax.set_ylim(bottom=0)
    ax.grid(False)

    return

#############################################################


def plot_pdf2d(samsx, samsy, pltype='kde', ncont=10, color=None, lwidth=2, mstyle='o', ax=None):
    if ax is None:
        ax = plt.gca()

    if pltype == 'kde':
        ngrid = 100
        a = 3.
        pnomx = np.mean(samsx)
        pdelx = np.std(samsx)
        pnomy = np.mean(samsy)
        pdely = np.std(samsy)

        x = np.linspace(pnomx - a * pdelx, pnomx + a * pdelx, ngrid)
        y = np.linspace(pnomy - a * pdely, pnomy + a * pdely, ngrid)
        X, Y = np.meshgrid(x, y)
        pgrid = np.vstack((X.flatten(), Y.flatten())).T  # pgrid.shape is (33^2,2)
        _, pdf = pdf_kde.get_pdf(np.vstack((samsx, samsy)).T, pgrid, verbose=0, method='Python')

        if color is None:
            ax.contour(X, Y, pdf.reshape(X.shape), ncont, linewidths=lwidth)
        else:
            ax.contour(X, Y, pdf.reshape(X.shape), ncont, colors=color, linewidths=lwidth)

    elif pltype == 'sam':
        ax.plot(samsx, samsy, color=color, marker=mstyle, linestyle='None')
    else:
        print("Plot type is not recognized. Exiting")
        sys.exit()

    ax.grid(False)

    return

#############################################################


def plot_pdfs(ind_show=[], samples_file='pchain.dat', plot_type='tri',
              names_file=None, burnin=0, every=1,
              nominal_file=None, nsam_xi=1, prangeshow_file=None):



    return figs, axarr
