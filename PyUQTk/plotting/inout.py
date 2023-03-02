#!/usr/bin/env python
#=====================================================================================
#
#                      The UQ Toolkit (UQTk) version 3.1.3
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
import os
import sys

try:
    import numpy as np
except ImportError:
    print("Numpy was not found.")

try:
    import matplotlib
except ImportError:
    print("Matplotlib was not found.")

import matplotlib.pyplot as plt

from pylab import *

# sys.path.append(os.environ['UQTK_INS'])
import PyUQTk.utils.colors as ut

rc('legend', loc='upper left', fontsize=12)
rc('lines', linewidth=1, color='r')
rc('axes', linewidth=3, grid=True, labelsize=22)
rc('xtick', labelsize=20)
rc('ytick', labelsize=20)

#############################################################
def parallel_coordinates(parnames, values, labels, savefig=[]):
    """
    Plots parallel coordinates.
    Arguments:
        * parnames : list of d parameter names
        * values   : (d,N) array of N data points with d parameters
        * labels   : list of N labels/categories, one per point
        * savefig  : figure name to save. If [], then show the plot
    """

    # Start the figure
    fig=figure(figsize=(14,7))
    fig.add_axes([0.1,0.15,0.8,0.8])
    ax = gca()

    # Categorize
    ulabels = np.unique(labels)
    n_labels = len(ulabels)

    # Set colors
    cmap = plt.get_cmap('prism')
    colors = cmap(np.arange(n_labels)*cmap.N/(n_labels+1))

    # Plot
    class_id = np.searchsorted(ulabels, labels)
    lines = plt.plot(values[:,:], 'ko-',ms=6,linewidth=0.7)
    [ l.set_color(colors[c]) for c,l in zip(class_id, lines) ]

    # Gridification
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_position(('outward', 5))
    ax.spines['bottom'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.yaxis.set_ticks_position('both')
    ax.xaxis.set_ticks_position('none')

    plt.xticks(np.arange(len(parnames)), parnames)
    plt.grid(axis='x', ls='-')

    leg_handlers = [ lines[np.where(class_id==id)[0][0]]
                    for id in range(n_labels)]
    ax.legend(leg_handlers, ulabels, frameon=False, loc='upper left',
                    ncol=len(labels),
                    bbox_to_anchor=(0, -0.03, 1, 0))

    # Show or save
    if (savefig==[]):
        plt.show()
    else:
        plt.savefig(savefig)
        plt.clf()


#############################################################

def plot_xx(d1,d2,parnames, values, labels, savefig=[]): #(x1,x2,inputs,labels,pnames,outfigdir='.'):
    """
    Plots one-dimension versus another with various labels.
    Arguments:
        * d1       : first dimension to plot
        * d2       : second dimension to plot
        * parnames : list of d parameter names
        * values   : (d,N) array of N data points with d parameters
        * labels   : list of N labels/categories, one per point
        * savefig  : figure name to save. If [], then show the plot
    """

    # Start the figure
    fig=figure(figsize=(12,12))
    fig.add_axes([0.1,0.15,0.8,0.8])
    ax = gca()

    # Categorize
    ulabels = np.unique(labels)
    n_labels = len(ulabels)

    # Set colors
    cmap = plt.get_cmap('prism')
    colors = cmap(np.arange(n_labels)*cmap.N/(n_labels+1))

    # Plot
    class_id = np.searchsorted(ulabels, labels)
    for id in range(n_labels):
        plt.plot(values[class_id==id,d1],values[class_id==id,d2], 'o',color=colors[id],ms=7,label=ulabels[id])



    ax.legend(frameon=False, loc='upper left',
                  ncol=len(labels),
                  bbox_to_anchor=(0, -0.06, 1, 0))

    ax.set_xlabel(parnames[d1])
    ax.set_ylabel(parnames[d2])

    # Show or save
    if (savefig==[]):
        plt.show()
    else:
        plt.savefig(savefig)

    return fig

#############################################################

def plot_xy(x,y,pname, outname, label='', savefig=[]):
    """
    Plots one array versus another.
    Arguments:
        * x        : array for x-axis
        * y        : array for y-axis
        * pname    : xlabel
        * outname  : ylabel
        * label    : legend
        * savefig  : figure name to save. If [], then show the plot
    """

    # Start the figure
    fig=figure(figsize=(12,8))
    ax = gca()

    # Plot
    plt.plot(x,y,'o',ms=11,label=label)

    # Set labels
    ax.set_xlabel(pname)
    ax.set_ylabel(outname)

    # Show or save
    if (savefig==[]):
        plt.show()
    else:
        plt.savefig(savefig)
        #plt.clf()


    return fig

#############################################################

def plot_xxy(x1,x2,y,pnames, outname, label='', savefig=[]):
    """
    Plots one array versus another.
    Arguments:
        * x1       : array for x1-axis
        * x2       : array for x2-axis
        * y        : array for y-axis
        * pnames   : list of xlabels
        * outname  : ylabel (vertical axis)
        * label    : legend
        * savefig  : figure name to save. If [], then show the plot
    """

    # Start the figure
    fig=figure(figsize=(12,8))
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(x1,x2,y,c='k',label=label)

    # Set labels
    ax.set_xlabel(pnames[0])
    ax.set_ylabel(pnames[1])
    ax.set_zlabel(outname)

    # Show or save
    if (savefig==[]):
        plt.show()
    else:
        plt.savefig(savefig)
        #plt.clf()


    return fig
