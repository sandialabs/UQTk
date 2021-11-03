#!/usr/bin/env python
#=====================================================================================
#
#                      The UQ Toolkit (UQTk) version 3.1.1
#                          Copyright (2021) NTESS
#                        https://www.sandia.gov/UQToolkit/
#                        https://github.com/sandialabs/UQTk
#
#     Copyright 2021 National Technology & Engineering Solutions of Sandia, LLC (NTESS).
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
import getopt
import math
import matplotlib.pyplot as plt
import scipy
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import xml.etree.ElementTree as ET # for the xml tree and parsing
import unittest
import random as random

import matplotlib.ticker as mticker


import PyUQTk.inference.postproc as uqtkinfpp
import PyUQTk.utils.pdf_kde as uqtk_kde

import PyUQTk.uqtkarray as uqtkarray

try:
    import PyUQTk.inference.postproc as uqtkinfpp
    import PyUQTk.inference.evidence_solvers as evd
except ImportError:
    print("Could not import the UQTk postprocessing tools")
    sys.exit(1)


try:
    import tools as tools
except ImportError:
    print("Could not import the tools module")

try:
    import graph_tools as graph_tools
except ImportError:
    print("Could not import the tools module")

try:
    import evidence
except ImportError:
    print("Could not import the evidence module")


################################################################################
help_string = """
Usage:
  post.py [-h] [--ix <xml_input_file_name>] [-p] [-g] [-d] [-v <verbosity>] [--interactive] [--stats <stats_type>] [--jpeg] [--evidence] [--gaussian][-w <write file name>]
What:
  Inference postprocessing.
Where
  -h = print help info
  --ix = XML input file name [default "input.xml"]
  -p = make posterior plots
  -g = make parameter graphs
  -d = calculate derivatives and graph
  -v = verbosity level [default 1]
  --interactive = show plots interactively [defaults to not]
  --stats = desired MCMC statistics to generate [none by default]
            corr: autocorrelation functions
            conv: convergence tests
            all : all
  --jpeg = option to save plots and figures as jpeg [default is pdf]
  --evidence = Computes the evidence values [default False]
  --gaussian = Computes the gaussian evidence values (not very reliable calculation)[default False]
  -w = name of the file that gets written to [necessary if computing evidence]
Assumes the following:
    * The mcmc file name is taken from the xml input file
    * The file is in ASCII format
    * The first column contains the MCMC step number
    * The second to last column contains the acceptance probability for the jump proposed in this step
    * The last column contains the log posterior of the state in this step
    * The columns in between contain the sampled states
    * Unless the argument labels == False, the first line contains labels for each column
"""

if __name__ == "__main__":
    #
    # Process inputs
    #
    try:
        opts,v_names = getopt.getopt(sys.argv[1:],"hpgv:w:d",["ix=","interactive","stats=", "jpeg", "evidence"])
    except getopt.GetoptError as err:
        print(str(err))
        print(help_string)
        sys.exit(1)

    # Default values
    input_file_name = "input.xml"
    labels_present = True
    np_kde = 50 #points per dimension for plotting distributions (from KDE)
    post_plots = False
    parameter_graphs = False
    demo_verbose = 1 # set > 0 to get more verbose output (and higher numbers generate even more output)
    interactive = False
    desired_stats=[] # default list of desired stats
    jpeg = False #whether to save figures as jpegs
    to_compute_evidence = False #default, to not calculate evidence
    to_compute_gaussian = False #default, to not calculate gaussian evidence
    write_file_name = "" #write file, needed for evidence output
    calc_deriv = False

    #read in from command line
    for o,a in opts:
        if o == "-h":
            print(help_string)
            sys.exit(0)
        elif o == "--ix":
            input_file_name = a
        elif o == "-v":
            demo_verbose = float(a)
        elif o == "-p":
            post_plots = True
        elif o == "-g":
            parameter_graphs = True
        elif o == "--interactive":
            interactive = True
        elif o == "--stats":
            desired_stats.append(a)
        elif o == "-d":
            calc_deriv = True
        elif o == "--jpeg":
            jpeg = True
        elif o == "--evidence":
            to_compute_evidence = True
        elif o == "--gaussian":
            to_compute_gaussian = True
        elif o == "-w":
            write_file_name = a
        else:
            assert False, "Unhandled command line parsing option. Use -h flag to get usage info."

    # Error checking
    for stats in desired_stats:
        if stats not in ["corr", "conv", "all"]:
            print("Invalid statistics type",stats,"requested")
            sys.exit(1)

    #open xml file
    try:
        full_tree = ET.parse(input_file_name)
        full_tree_root = full_tree.getroot()
        print("\nRead in the file",input_file_name,"for run settings.")
    except:
        print("Error reading in xml tree from input file. Make sure the file",input_file_name,)
        print("is in the current directory and has proper xml syntax.")
        sys.exit(1)

    # Get the case information from the xml file
    run_info_tree = full_tree_root.find("case_info")
    my_case = run_info_tree.get("case")
    run_mode = run_info_tree.get("mode")
    if run_mode == "multi-model": # Run multiple models
        my_models = [model.get("name") for model in run_info_tree.find("models").findall("model")]
    else:  # Run for only first model listed
        my_models = [run_info_tree.find("models").find("model").get("name")]

    if my_models == []:
        print("\nERROR: no models found in xml input file\n")
        sys.exit(1)

    # Get xml tree section with case information
    case_info_tree = full_tree_root.find(".//case_types/" + my_case)
    if case_info_tree is None:
        print("\nERROR: case info not found in xml input file\n")
        sys.exit(1)

    if demo_verbose > 0:
        print("\nPostprocessing case",my_case,"with the models",my_models)

    # Read in data from file
    input_file_name = case_info_tree.get("input_file")
    data_obs = np.genfromtxt(input_file_name,comments='#',delimiter=',')
    coeff_input_file_name = case_info_tree.get("coeff_input_file")
    real_coeff = np.genfromtxt(coeff_input_file_name,comments='#',delimiter=',')
    #get error_level to use as sigma in fitting
    error_file_name = case_info_tree.get("error_file")
    error_level = np.genfromtxt(error_file_name,comments='#',delimiter=',')

    #get error_level to use as sigma in fitting
    get_data_info_tree = full_tree_root.find("get_data_info")
    error_level = float(get_data_info_tree.get("error_level"))

    #lists to store for multiple models
    model_name_list = []
    gaussian_evidence_list = []
    importance_evidence_list = []
    derivative_mean_list = []
    derivative_std_list = []
    y_values_mean_list = []
    y_values_std_list = []

    # Run for all models in list:
    for my_model in my_models:
        if demo_verbose > 0:
            print("\nRunning for model",my_model)

        #make the model object
        model = tools.make_model_object(my_model, case_info_tree, data_obs, error_level)

        # Create an easier to use list of parameter names that is sorted according to the same
        # index used in all other species arrays.
        param_names_dict = model.model_info["param_names_dict"]
        param_names_list = sorted(param_names_dict.keys(), key=param_names_dict.__getitem__)

        #
        # Process MCMC writeSamples
        #

        # # Post processing settings
        #get relavent info from model object
        n_skip = model.model_info["n_skip"]
        stride = model.model_info["stride"]

        if demo_verbose > 1:
            print("n_skip:",n_skip)
            print("stride:",stride)

        # Base name of file for outputting results
        samples_file_name = model.model_info["output_file"] #MCMC file name
        out_file_base = samples_file_name + ".sk" + str(n_skip) + ".st" + str(stride)

        # Import variables of interest from the MCMC data file
        print("UQTk extract all vars")
        all_samples, v_names = uqtkinfpp.extract_all_vars(samples_file_name,n_skip,demo_verbose,stride,labels=labels_present)

        # Get statistics
        print("UQTk get_mcmc_stats")
        map_params = uqtkinfpp.get_mcmc_stats(all_samples,v_names,out_file_base,demo_verbose,desired_stats)


        print("\nMAP parameter set:")
        for idx in range(len(map_params)):
            print(v_names[idx],":",map_params[idx])

        # Save MAP parameter info to a file
        MAP_info = np.stack((v_names,map_params),axis=-1)
        np.savetxt(out_file_base+"_map_params.dat",MAP_info,fmt='%-8s')

        # Plot posteriors
        # first and last three columns do not contain variable samples
        if post_plots:
            uqtkinfpp.plot_all_posteriors(all_samples[:,1:(1+len(v_names))],v_names,np_kde,out_file_base,demo_verbose,dense=True)

        # Make graphs for each parameter
        if parameter_graphs:
            graph_tools.plot_parameter_graphs(out_file_base + "_parameter_graphs",
                all_samples, v_names, interactive = interactive, jpeg = jpeg)

        # Make graphs for the agreement between data and observations:
        model.model_info["fit_params"] = map_params

        y_pred = model.prediction()

        data_pred = model.make_pred_array(y_pred)

        #use helper function to make various plots
        graph_tools.model_agreement_graph_xy(out_file_base, data_obs, data_pred,
            my_model, interactive=interactive, jpeg=jpeg)
        graph_tools.model_agreement_graph_real_data_pred_line(out_file_base,
            map_params, real_coeff, data_obs, my_model, interactive=interactive,jpeg=jpeg)
        graph_tools.plot_with_chains(out_file_base, map_params,
            all_samples[:, 1:-2], real_coeff, data_obs, my_model,interactive=interactive,jpeg=jpeg)

        #graph over many samples with std
        # Select samples from MCMC chain by randomly choosing row indices in MCMC sample file
        n_samp = 2000
        mcmc_ind_subset = np.random.choice(all_samples[:,1:-2].shape[0],n_samp,replace=False)
        x = np.linspace(0,1)

        #make array to save all derivatives
        y_values = np.zeros((n_samp, len(x)))

        for i_samp in range(n_samp):
            samp_idx = mcmc_ind_subset[i_samp]
            #calculate average
            y = tools.evaluate_function(all_samples[samp_idx,1:-2], x)
            y_values[i_samp,:] = y

        # Compute ensemble mean and standard deviation
        # y_values_mean_list.append(np.mean(y_values,axis=0))
        y_values_std_list.append(np.std(y_values,axis=0))

        #we are making y_values_mean be the map parameters instead
        y_values_mean_list.append(tools.evaluate_function(map_params, x))

        #Calculate and plot derivatives
        if calc_deriv:
            # Select samples from MCMC chain by randomly choosing row indices in MCMC sample file
            n_samp = 2000
            mcmc_ind_subset = np.random.choice(all_samples[:,1:-2].shape[0],n_samp,replace=False)
            x = np.linspace(0,1)

            #make array to save all derivatives
            derivatives = np.zeros((n_samp, len(x)))

            for i_samp in range(n_samp):
                samp_idx = mcmc_ind_subset[i_samp]
                #calculate average
                y = model.compute_derivative(x, all_samples[samp_idx,1:-2])
                derivatives[i_samp,:] = y

            # Compute ensemble mean and standard deviation
            derivative_mean_list.append(np.mean(derivatives,axis=0))
            derivative_std_list.append(np.std(derivatives,axis=0))

            model_name_list.append(my_model)

            #find true value of derivative
            true_derivative = tools.evaluate_derivative(real_coeff, x)

            #Plot for single-model
            if run_mode != "multi-model":
                graph_tools.plot_with_bars(out_file_base+"_derivatives", x,
                    derivative_mean_list, derivative_std_list,
                    "Uncertainty in the derivative of \n" + my_case + "  " + my_model,
                    "$\partial y / \partial x$", model_name_list, true_derivative,
                    interactive=False, jpeg = jpeg)

            #multi model will be handled at the end

        #calcaulate the evidence values for each model
        if to_compute_evidence:
            (gaussian_evidence, importance_evidence) = evidence.compute_evidence_single_model(model, all_samples, v_names)
            gaussian_evidence_list.append(gaussian_evidence)
            importance_evidence_list.append(importance_evidence)

    if run_mode == "multi-model":
        #find true value of derivative
        true_solution = tools.evaluate_function(real_coeff, x)
        #plot multiple fits with std
        graph_tools.plot_with_bars(my_case+"_all_fits_with_error", x,
            y_values_mean_list, y_values_std_list,
            "Uncertainty of the fit\n" + my_case,
            "$y$", my_models, true_solution, interactive=False, jpeg = jpeg)

        #plot with shaded error bars
        graph_tools.plot_shaded_error(my_case+"_all_fits_with_error_shaded", x,
            y_values_mean_list, y_values_std_list,
            "Uncertainty of the fit\n" + my_case,
            "$y$", my_models, real_coeff, data_obs, interactive=False, jpeg = jpeg)

        #plot multiple fits all on one graph, no error bars
        graph_tools.plot_all_models(my_case+"_all_fits", x,
            y_values_mean_list,
            "Fits of all models\n" + my_case,
            "$y$", my_models, real_coeff, data_obs, interactive=False, jpeg = jpeg)


        #graph derivative for multiple models, with error bars
        if calc_deriv:
            #find true value of derivative
            true_derivative = tools.evaluate_derivative(real_coeff, x)

            graph_tools.plot_with_bars(my_case+"_derivatives", x,
                derivative_mean_list, derivative_std_list,
                "Uncertainty in the derivative\n" + my_case,
                "$\partial y/ \partial x$", model_name_list, true_derivative,
                interactive=False, jpeg = jpeg)

        if to_compute_evidence: #plot the evidence values

            #plot importance evidence
            graph_tools.plot_evidence(my_case + "_importance_evidence", importance_evidence_list,
                my_models, "Importance Evidence Values for " + my_case,
                interactive = False, jpeg = jpeg)

        #plot gaussian evidence
        #not very reliable method to calculate evidence values, and therefore not default
        if to_compute_gaussian:
            graph_tools.plot_evidence(my_case + "_gaussian_evidence", gaussian_evidence_list,
                my_models, "Gaussian Evidence Values for " + my_case,
                interactive = False, jpeg = jpeg)
