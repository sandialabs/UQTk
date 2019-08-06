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
from __future__ import print_function

from scipy.interpolate import interp1d
from scipy.stats.mstats import mquantiles

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

try:
    import numpy as np
except ImportError:
    print("Numpy was not found. ")

try:
    import matplotlib as mpl
    import matplotlib.pyplot as plt
except ImportError:
    print("Matplotlib was not found. ")


def plot_vars(ax, xdata, ydata, variances=None, ysam=None, stdfactor=1.,
              varlabels=None, varcolors=None, grid_show=True,
              connected=True, interp=None, offset=(None, None)):

    if varcolors is None:
        cmap = mpl.cm.Greys

    if variances is None:
        variances = np.empty((0, 0))
    if ysam is not None:
        plot_ens(ax, xdata, ysam,
                 connected=connected, interp=interp,
                 color='r', lw=0.3, offset=offset)

    shift, scale = offset
    if shift is None:
        shift = np.zeros_like(ydata)
    if scale is None:
        scale = np.ones_like(ydata)

    ydata_ = (ydata - shift) / scale

    if variances.shape[1] > 0:
        variances_ = variances / (scale.reshape(-1, 1) * scale.reshape(-1, 1))
    else:
        variances_ = variances
    nvariances = variances_.shape[1]

    if interp is not None:
        connected = True
        xmin = xdata.min()
        xmax = xdata.max()
        xdel = xmax - xmin
        ngr = 100
        xdata_ = np.linspace(xmin - 0.0 * xdel, xmax + 0.0 * xdel, ngr)

        interp_fcn = interp1d(xdata, ydata_, kind=interp)
        ydata_ = interp_fcn(xdata_)

        tmp = np.empty((ngr, nvariances))
        for i in range(nvariances):
            interp_fcn = interp1d(xdata, variances_[:, i], kind=interp)
            tmp[:, i] = interp_fcn(xdata_)
        variances_ = tmp.copy()

    else:
        xdata_ = xdata.copy()

    cvars = np.cumsum(variances_, axis=1)
    nvars = variances_.shape[1]
    normalize = mpl.colors.Normalize(vmin=0.0, vmax=0.5)
    if varlabels is None:
        varlabels = ['Var ' + str(i) for i in range(1, nvars + 1)]

    if connected:
        plt.gca().plot(xdata_, ydata_, color='orange',
                   marker='None', linestyle='-',
                   label='Mean prediction', zorder=10000)
        for ii in range(nvars):
            if varcolors is None:
                varcolor = cmap(normalize(0.1 + ii * 0.3 / nvars))
            else:
                varcolor = varcolors[ii]
            plt.gca().fill_between(xdata_,
                               ydata_ - stdfactor * np.sqrt(cvars[:, ii]),
                               ydata_ + stdfactor * np.sqrt(cvars[:, ii]),
                               color=varcolor,
                               label=varlabels[ii], zorder=1000 - ii)
    else:
        plt.gca().plot(xdata_, ydata_, color='orange',
                   marker='o', linestyle='None',
                   label='Mean prediction', zorder=100000)
        for ii in range(nvars):
            if varcolors is None:
                varcolor = cmap(normalize(0.1 + ii * 0.3 / nvars))
            else:
                varcolor = varcolors[ii]
            plt.gca().errorbar(xdata_, ydata_,
                           yerr=stdfactor * np.sqrt(cvars[:, ii]),
                           color=varcolor,
                           fmt='o', elinewidth=13,
                           label=varlabels[ii], zorder=10000 - ii)

    # set the current axis to ax
    plt.sca(ax)

    ax.grid(grid_show)
    # plt.close(fig)
    # plt.show()

    return ax


def plot_shade(ax, xdata, ydata, nq=51, cmap=mpl.cm.BuGn,
               bounds_show=False, grid_show=True):
    nx = xdata.shape[0]
    assert(nx == ydata.shape[0])

    mq = mquantiles(ydata, prob=[float(i + 1) / float(nq)
                                 for i in range(nq - 1)], axis=1)
    print(mq.shape, ydata.shape)
    plt.sca(ax)

    normalize = mpl.colors.Normalize(vmin=0.01, vmax=0.5)

    for k in range(int(nq / 2)):
        plt.fill_between(xdata, mq[:, k], mq[:, k + 1],
                         color=cmap(normalize(0.01 + k * .02)))
    for k in range(int(nq / 2), nq - 2):
        plt.fill_between(xdata, mq[:, k], mq[:, k + 1],
                         color=cmap(normalize(0.5 - (k - nq / 2) * 0.02)))
    if bounds_show:
        plt.plot(xdata, mq[:, 0], linewidth=2, color="grey")
        plt.plot(xdata, mq[:, -1], linewidth=2, color="grey")

    ax.grid(grid_show)
    # plt.close(fig)
    # plt.show()

    return ax


def plot_ens(ax, xdata, ydata, color='b', lw=2, ms=1,
             grid_show=True, label='',
             connected=True, interp=True, offset=(None, None)):

    nx = ydata.shape[0]
    nsam = ydata.shape[1]

    shift, scale = offset
    if shift is None:
        shift = np.zeros_like(ydata)
    if scale is None:
        scale = np.ones_like(ydata)

    shift, scale = offset
    if shift is None:
        shift = np.zeros((nx, 1))
    else:
        shift = shift.reshape((nx, 1))
    if scale is None:
        scale = np.ones((nx, 1))
    else:
        scale = scale.reshape((nx, 1))
    ydata_ = (ydata - shift) / scale

    if interp is not None:
        connected = True
        xmin = xdata.min()
        xmax = xdata.max()
        xdel = xmax - xmin
        ngr = 100
        xdata_ = np.linspace(xmin - 0.0 * xdel, xmax + 0.0 * xdel, ngr)

        tmp = np.empty((ngr, nsam))
        for i in range(nsam):
            interp_fcn = interp1d(xdata, ydata_[:, i], kind=interp)
            tmp[:, i] = interp_fcn(xdata_)
        ydata_ = tmp.copy()

    else:
        xdata_ = xdata.copy()

    if connected:
        for i in range(nsam):
            plt.gca().plot(xdata_, ydata_[:, i],
                       color=color, linewidth=lw)
    else:
        for i in range(nsam):
            plt.gca().plot(xdata_, ydata_[:, i], color=color,
                       linestyle='None', marker='o', markersize=ms,
                       label=label, zorder=1000)

    # set the current axis to ax
    plt.sca(ax)

    ax.grid(grid_show)
    # plt.close(fig)
    # plt.show()

    return ax
