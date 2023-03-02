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
#
# Does statistical analysis on samples from an MCMC chain.
from __future__ import print_function # To make print() in Python 2 behave like in Python 3

import os
import sys
import string
import numpy as np
import math
import matplotlib.pyplot as plt
import random as random

try:
    import tools as tools
except ImportError:
    print("Could not import the tools module")


################################################################################
def log_tick_formatter(val, pos=None):
    '''
    see https://github.com/matplotlib/matplotlib/issues/209
    '''
    return "{:.2e}".format(10**val)
################################################################################
def save_plot(name_base, jpeg = False):
    """
    saves the plot as either a pdf or jpegs
    to run after the plot is made
    inputs:
        name_base: name base of the file
        jpeg: set to True to save the plot as a jpeg [default False]
    output:
        saved plot
    """
    if jpeg:
        try:
            plt.savefig(name_base + ".jpg")
        except:
            print('Cannot save as jpeg')
            print('Might not have the required library')
            print('Please install pillow. Can be done by running:')
            print('     pip install pillow')
    else:
        plt.savefig(name_base + ".pdf")
################################################################################
def plot_parameter_graphs(name_base, all_samples, v_names, interactive = False, jpeg = False):
    '''
     Creates the plots of the chains of each parameter

     inputs:
        name_base: the name base with which to save the graph
        all_samples: 2D numpy array with the chain samples of all the parameter_sample
        v_names: a list of all the variable names
        interactive: Set to True to have the plots pop up when running the script [default False]
        jpeg: set to True to save the plot as a jpeg [default False]
     output:
        shows the chain plots of each parameter, also saves as a pdf or jpeg
    '''
    #get names for graphs
    names = v_names + ["Acceptance Probability", "Log Posterior"]
    plt.figure()
    plt.suptitle("Parameter graphs")
    for i in range(len(names)):
        plt.subplot(math.ceil(len(names) / 2.0) , 2 ,i + 1)
        plt.plot(all_samples[:,0],all_samples[:, i + 1])
        plt.title(names[i])
        plt.tight_layout()
        plt.tick_params(axis='both', labelsize=10)

    save_plot(name_base, jpeg)

    if interactive:
        plt.show()

    plt.close()
################################################################################
def model_agreement_graph_xy(name_base, data_obs, data_pred, my_model, interactive=False, jpeg = False):
    '''
     Create a scatter plot of observed vs. predicted x values

     inputs:
        name_base: the name base with which to save the graph
        data_obs: observational data x, y, one experiment per row
        data_pred: fit function prediction data, x, y, one experiment per row
        my_model: name of the model to graph
        interactive: Set to True to have the plots pop up when running the script [default False]
        jpeg: set to True to save the plot as a jpeg [default False]
     output:
        shows model agreement xy graph, also saves as a pdf or jpeg
    '''

    fig = plt.figure(figsize=(8,6))
    ax = fig.add_axes([0.10, 0.10, 0.8, 0.8])

    plt.scatter(data_obs[:,1], data_pred[:,1], c='r', marker='o')

    # Add a line x = y, courtesy of https://stackoverflow.com/questions/25497402/adding-y-x-to-a-matplotlib-scatter-plot-if-i-havent-kept-track-of-all-the-data
    lims = [
        np.min([ax.get_xlim(), ax.get_ylim()]),  # min of both axes
        np.max([ax.get_xlim(), ax.get_ylim()]),  # max of both axes
    ]

    # now plot both limits against eachother
    ax.plot(lims, lims, 'k-')
    ax.set_aspect('equal')
    ax.set_xlim(lims)
    ax.set_ylim(lims)

    # Label Axes
    plt.xlabel("$y$ Observed", size=16)
    plt.ylabel("$y$ Predicted", size=16)
    # Add title
    plt.suptitle("Predicted vs. Observed $y$\n " + my_model, size=20)
    # Change tick size
    plt.tick_params(axis='both', labelsize=14)
    # Pad tick labels
    plt.gca().tick_params(pad=6)
    # Save figure
    save_plot(name_base + "_model_data_agreement_xy", jpeg = jpeg)\

    # Show figure
    if interactive:
        plt.show()

    plt.close()

    return
###########################################################
def model_agreement_graph_real_data_pred_line(name_base, coeff, real_coeff, data_obs,my_model,
    interactive=False, jpeg = False):
    '''
     Create a plot with data, prediction line, and real line

     inputs:
        name_base: the name base with which to save the graph
        coeff: MAP parameters, of some kind of coeffceint predictions
        real_coeff: teh actual coefficents saved in get_data.py
        data_obs: 2D numpy array of the real data to also plot.
        my_model: name of the model to graph
        interactive: Set to True to have the plots pop up when running the script [default False]
        jpeg: set to True to save the plot as a jpeg [default False]
     output:
        Saves plot as pdf or jpg, can also show plot
    '''

    fig = plt.figure(figsize=(8,6))
    ax = fig.add_axes([0.10, 0.10, 0.8, 0.8])

    #get prdiction fit
    x = np.linspace(0,1)
    y = tools.evaluate_function(coeff, x)

    #get real fit
    real_y = tools.evaluate_function(real_coeff, x)

    #plot everything
    plt.plot(data_obs[:,0], data_obs[:,1], 'o', label="data")
    plt.plot(x, y, label = "Inferred")
    plt.plot(x, real_y, '--', label = "Truth")
    plt.title(my_model + " fit")
    plt.legend()

    #save plots
    save_plot(name_base + "_model_data_agreement_xy_with_real", jpeg = jpeg)

    # Show figure
    if interactive:
        plt.show()

    plt.close()

    return
################################################################################
def plot_with_chains(name_base, map_params, chains, real_coeff, data_obs, my_model,
    num_chains = 100, interactive=False, jpeg = False):
    '''
     Create a plot with data, prediction lines(with many samples of the chain
     and the MAP parameters), and real line

     inputs:
        name_base: the name base with which to save the graph
        map_params: MAP parameters, of some kind of coeffceint predictions
        chains: The MCMC chain, include only the parts that are coefficents (typically (1, -3)),
            not sample number, sigma, posterior probability and acceptance rate
        real_coeff: the actual coefficents saved in get_data.py
        data_obs: 2D numpy array of the real data to also plot.
        my_model: name of the model to graph
        num_chains: number of mcmc chains to plot. [default 100]
        interactive: Set to True to have the plots pop up when running the script [default False]
        jpeg: set to True to save the plot as a jpeg [default False]
     output:
        Saves plot as pdf or jpg, can also show plot
    '''

    fig = plt.figure(figsize=(8,6))
    ax = fig.add_axes([0.10, 0.10, 0.8, 0.8])
    x = np.linspace(0,1)

    #plot data
    plt.plot(data_obs[:,0], data_obs[:,1], 'yo', label="data")

    #plot 100 random samples from MCMC chain
    for i in range(num_chains):
        index = random.randint(0, len(chains) - 1)
        coeff = chains[index, :]
        y = tools.evaluate_function(coeff, x)
        plt.plot(x,y,'b', alpha = 0.1)

    #plot real equation
    real_y = tools.evaluate_function(real_coeff, x)
    plt.plot(x,real_y, 'g', label = "True Polynomial")

    #plot map value equation
    map_y = tools.evaluate_function(map_params, x)
    plt.plot(x,map_y, 'r', label = "MAP coefficeints")

    plt.legend()
    plt.title(my_model + " fit")

    #save plot
    save_plot(name_base + "_model_data_agreement_xy_many", jpeg = jpeg)

    # Show figure
    if interactive:
        plt.show()

    plt.close()

################################################################################
def plot_with_bars(name_base, x_vals, mean_list, stdv_list, title, y_label,
    model_name_list, true_solution, interactive=False, jpeg = False):
    '''
     Plots all models with error bars,
        can be any values passed in, (i.e. true values, deivatives, etc. )
     inputs:
        file_name: the name with which to save the graph
        x_vals: 1D numpy array with x values
        mean_list: 2D numpy array with mean of y_axis over all models, i.e. the derivative
        stdv_list: 2D numpy array with standard deviation of y_axis over all models, i.e. the derivative
        title: a string for the title of the graph
        y_label: the label for the y-axis (x-axis will be x)
        model_name_list: a list of the model names in order, will appear in the legend
        true_solution: 1D array of the true solution
        interactive: Set to True to have the plots pop up when running the script [default False]
        jpeg: set to True to save the plot as a jpeg [default False]
     output:
        Makes graph of the quanity with error bars, saves as pdf or jpeg
    '''

    fig = plt.figure(figsize=(8,6))
    ax = fig.add_axes([0.20, 0.10, 0.75, 0.8])

    plt.plot(x_vals, true_solution, color = "black", linewidth=2.0, label= "True solution")
    for i in range(len(mean_list)):
        plt.errorbar(x_vals,mean_list[i],yerr=stdv_list[i],linewidth=2.0,capsize=5,capthick=2, label= model_name_list[i])

    # Label Axes
    plt.xlabel("$x$", size=16)
    plt.ylabel(y_label, size=16)
    # Add title
    plt.title(title, size=20)
    # Change tick size
    plt.tick_params(axis='both', labelsize=14)
    # Pad tick labels
    plt.gca().tick_params(pad=6)
    # Create legend
    plt.legend(loc='upper right', prop={'size':12})
    # Save figure
    save_plot(name_base, jpeg = jpeg)

    # Show figure
    if interactive:
        plt.show()

    plt.close()

    return
################################################################################
def plot_all_models(name_base, x_vals, mean_list, title, y_label, model_name_list,
    real_coeff, data_obs, interactive=False, jpeg = False):
    '''
     Plots all models with real solution and data
     inputs:
        file_name: the name with which to save the graph
        x_vals: 1D numpy array with x values
        mean_list: 2D numpy array with mean of y_axis over all models, i.e. the derivative
        title: a string for the title of the graph
        y_label: the label for the y-axis (x-axis will be x)
        model_name_list: a list of the model names in order, will appear in the legend
        real_coeff: a list of the real coefficents to evaulate the function
        data_obs: 2D numpy array of the real data to also plot.
        interactive: Set to True to have the plots pop up when running the script [default False]
        jpeg: set to True to save the plot as a jpeg [default False]
     output:
        Makes graph for the model, true solution, and data, saves as pdf or jpeg
    '''

    fig = plt.figure(figsize=(8,6))
    # ax = fig.add_axes([0.20, 0.10, 0.75, 0.8])

    #make the x values you are plotting over be the entire range you evaulate the true solution at
    x_values = np.linspace(0, 1)
    plt.plot(x_values, tools.evaluate_function(real_coeff, x_values),
        linewidth=2.0, label= "True solution", color = "black")
    plt.scatter(data_obs[:,0], data_obs[:,1], color = "black")

    for i in range(len(mean_list)):
        plt.plot(x_vals,mean_list[i],linewidth=2.0,label= model_name_list[i])

    # Label Axes
    plt.xlabel("$x$", size=16)
    plt.ylabel(y_label, size=16)
    # Add title
    plt.title(title, size=20)
    # Change tick size
    plt.tick_params(axis='both', labelsize=14)
    # Pad tick labels
    plt.gca().tick_params(pad=6)
    # Create legend
    plt.legend(loc='upper right', prop={'size':12})
    # Save figure
    save_plot(name_base, jpeg = jpeg)

    # Show figure
    if interactive:
        plt.show()

    plt.close()

    return

################################################################################
def plot_evidence(name_base, evidence_list, my_models, title,
    interactive = False, jpeg = False):
    '''
     Plots the evidence values in a bar graph

     inputs:
        name_base: the name base for the file to save
        evidence_list: a list of evidence values to plot
        my_models: a list of the name of the models, to label the evidence values
        title: title for the plot
        interactive: Set to True to have the plots pop up when running the script [default False]
        jpeg: set to True to save the plot as a jpeg [default False]
     output:
        plots evidence value bar graph, also saves as a pdf
        '''
    fig = plt.figure(figsize=(8,6))

    plt.bar(np.arange(len(evidence_list)), evidence_list)

    #label the x-axis
    plt.xticks(np.arange(len(evidence_list)), my_models, rotation='vertical')
    # Tweak spacing to prevent clipping of tick-labels
    plt.subplots_adjust(bottom=0.15)

    #label the y-axis
    plt.ylabel('log evidence')

    # Add title
    plt.suptitle(title, size=20)

    # Save figure
    save_plot(name_base, jpeg = jpeg)

     # Show figure
    if interactive:
        plt.show()

    plt.close()

################################################################################
def plot_shaded_error(name_base, x_vals, mean_list, stdv_list, title, y_label,
    model_name_list, real_coeff, data_obs, interactive=False, jpeg = False):
    '''
     Create a plot of all models with error
     also plots true solution, and data points

     inputs:
        name_base: the name base for the file to save
        x_vals: 1D numpy array with x valuesthe x values
        mean_list: 2D numpy array with mean of the value to plot over all models
        std_list: 2D numpy array with standard deviation of the value over all models
        title: a string for the title of the graph
        y_label: the label of the y axis
        model_name_list: a list of the model names in order
        real_coeff: the real coefficents of the model, to calculate the true solution
        data_obs: observational data x, y, one experiment per row
        interactive: Set to True to have the plots pop up when running the script [default False]
        jpeg: set to True to save the plot as a jpeg [default False]
     output:
        shows a plot over all the models with the errors shaded in
         also saves as a pdf
        '''
    fig = plt.figure(figsize=(8,6))
    ax = fig.add_axes([0.20, 0.10, 0.75, 0.8])

    #make the x values you are plotting over be the entire range you evaulate the true solution at
    x_values = np.linspace(0, 1)
    plt.plot(x_values, tools.evaluate_function(real_coeff, x_values),linewidth=2.0, label= "True solution", color = 'black')
    plt.scatter(data_obs[:,0], data_obs[:,1])

    #set the order of colors, so they line up with the default.
    colors = (["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf"])

    for i in range(len(mean_list)):
        plt.plot(x_vals,mean_list[i], '-', color=colors[i], label= model_name_list[i])
        plt.fill_between(x_vals,mean_list[i] - stdv_list[i], mean_list[i] + stdv_list[i], color = colors[i], alpha = 0.4)


    # Label Axes
    plt.xlabel("$x$", size=16)
    plt.ylabel(y_label, size=16)
    # Add title
    plt.title(title, size=20)
    # Change tick size
    plt.tick_params(axis='both', labelsize=14)
    # Pad tick labels
    plt.gca().tick_params(pad=6)
    # Create legend
    plt.legend(loc='upper right', prop={'size':12})

    # Save figure
    save_plot(name_base, jpeg = jpeg)

    # Show figure
    if interactive:
        plt.show()

    plt.close()

    return
