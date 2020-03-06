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
# fitting the polynomial for example

from __future__ import print_function # To make print() in Python 2 behave like in Python 3


import math
import sys
import getopt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import xml.etree.ElementTree as ET # for the xml tree and parsing



try:
    import numpy as np
except ImportError:
    print("Numpy module could not be found")


try:
    import PyUQTk.inference as uqtkinf
except ImportError:
    print("PyUQTk inference module not found")

try:
    import tools as tools
except ImportError:
    print("Could not import the tools module")
################################################################################
help_string = """
Usage:
  fit.py [-h] [-v <verbosity>] [--ix <xml_input_file_name>] -w <output_file_name>
What:
  Demonstrate fit of thermodynamic properties of redox-active metal oxide
Where
  -h = print help info
  -v = verbosity level [default 1]
  --ix = XML input file name [default "input.xml"]
  -w = name of output file that summary results are appended to [default "output.txt"]
"""

if __name__ == "__main__":
    #
    # Process inputs
    #
    try:
        opts,extra_arguments = getopt.getopt(sys.argv[1:],"hv:w:",["ix="])

    except getopt.GetoptError as err:
        print(str(err))
        print(help_string)
        sys.exit(1)

    # Default values
    demo_verbose = 1 # set > 0 to get more verbose output (and higher numbers generate even more output)
    output_file_name = "output.txt"
    input_file_name = "input.xml"

    #read in command line inputs
    for o,a in opts:
        if o == "-h":
            print(help_string)
            sys.exit(0)
        elif o == "-v":
            demo_verbose = float(a)
        elif o == "--ix":
            input_file_name = a
        elif o == "-w":
            output_file_name = a
        else:
            assert False, "Unhandled command line parsing option. Use -h flag to get usage info."


    #error checking
    if output_file_name is "":
        print("\nError: Output file name must be specified.\n")
        print(help_string)
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

    # open output file to append to
    fd = open(output_file_name, 'a')

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
        print("\nRunning case ",my_case," with the models",my_models)

    # Read in data from file that was made in get_data.py
    input_file_name = case_info_tree.get("input_file")
    input_data = np.genfromtxt(input_file_name,comments='#',delimiter=',')

    #get error_level to use as sigma in fitting
    error_file_name = case_info_tree.get("error_file")
    error_level = np.genfromtxt(error_file_name,comments='#',delimiter=',')


    # Run for all models in list:
    for my_model in my_models:
        if demo_verbose > 0:
            print("\nRunning for model",my_model)

        #make model object
        model = tools.make_model_object(my_model, case_info_tree, input_data, error_level)


        #
        # Infer the model parameters
        #
        inference_out = tools.fit_dram(model)

        #print out information about the MCMC run
        if demo_verbose > 0:
            print("MCMC sample size:", inference_out['chain'].shape)
            print("Overall acceptance rate:", inference_out['accr'])
            print("Fraction of samples outside of prior bounds:", 1.0 - inference_out['accb'])
            print("MAP parameter set:")

            # Create an easier to use list of parameter names that is sorted according to the same
            # index used in all other species arrays.
            param_names_dict = model.model_info["param_names_dict"]
            param_names_list = sorted(param_names_dict.keys(), key=param_names_dict.__getitem__)

            for p_name in param_names_list:
                print(p_name,":",'{: .6e}'.format(inference_out['cmap'][param_names_dict[p_name]]))

        # write to output file
        fd.write("\nModel:" + my_model +  "\n")
        fd.write("MCMC sample size: " +  str(inference_out['chain'].shape) + "\n")
        fd.write("Overall acceptance rate: " + str(inference_out['accr']) + "\n")
        fd.write("Fraction of samples outside of prior bounds:" + str(1.0 - inference_out['accb']) + "\n")
        fd.write("MAP parameter set: \n")
        for p_name in param_names_list:
            fd.write(p_name + ":" + str('{: .6e}'.format(inference_out['cmap'][param_names_dict[p_name]])) + "\n")

        #get output file name to write samples to
        output_file = model.model_info["output_file"]

        # transform last two columns from [like, prior] to posteriors
        meta_info = tools.meta_info_correct_posterior(inference_out['minfo'])
        # Write out sample files in proper format for UQTk postprocessing
        tools.writeSamples(output_file,inference_out['chain'],[inference_out['cmap'],inference_out['pmap']],meta_info,param_names_list)

    # Close output file
    fd.close()
