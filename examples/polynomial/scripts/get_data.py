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
# Gets data form a random polynomial with random noise
from __future__ import print_function # To make print() in Python 2 behave like in Python 3


import math
import sys
import getopt
import random
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import xml.etree.ElementTree as ET # for the xml tree and parsing

try:
    import tools as tools
except ImportError:
    print("Could not import the tools module")


try:
    import numpy as np
except ImportError:
    print("Numpy module could not be found")


try:
    import PyUQTk.inference as uqtkinf
except ImportError:
    print("PyUQTk inference module not found")

################################################################################
def add_noise(data, level):
    """
    Add noise to a list of data

    inputs:
        data: list to add to
        level: size of the std of the Noise

    outputs:
        list of same size, with noise added
    """

    error_data = data + np.random.normal(scale = level, size = np.shape(data))

    return error_data

################################################################################
help_string = """

Usage:
    get_data.py [-h] [-g] [--ix <input.xml>] [-e]
What:
    generates data from a polynomial with noise to be used for the tutorial
Where
  -h = print help info
  --ix = XML input file name [default "input.xml"]
  -g = to show the graph of the function and chosen data points [default False]
  -e = to run with the example polynomial. Will still give random points along that polynomial.
Assumes the following:
    * assumes there's an xml file with the following components
        size_range
        power
        error_level
        number_of_points
        input_file
        coeff_input_file

    * See example in input.xml
"""

if __name__ == "__main__":

    #
    # Process inputs
    #
    try:
        opts,extra_arguments = getopt.getopt(sys.argv[1:],"hge",["ix="])

    except getopt.GetoptError as err:
        print(str(err))
        print(help_string)
        sys.exit(1)

    # Default values
    output_file_name = "output.txt"
    input_file_name = "input.xml"
    to_graph = False
    example = False

    for o,a in opts:
        if o == "-h":
            print(help_string)
            sys.exit(0)
        elif o == "--ix":
            input_file_name = a
        elif o == "-g":
            to_graph = True
        elif o == "-e":
            example = True
        else:
            assert False, "Unhandled command line parsing option. Use -h flag to get usage info."


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

    # Get xml tree section with case information
    case_info_tree = full_tree_root.find(".//case_types/" + my_case)
    if case_info_tree is None:
        print("\nERROR: case info not found in xml input file\n")
        sys.exit(1)

    # Read in file to save data and coefficents to
    input_file_name = case_info_tree.get("input_file")
    coeff_input_file_name = case_info_tree.get("coeff_input_file")
    error_file_name = case_info_tree.get("error_file")

    #get relatvent info to make data
    get_data_info_tree = full_tree_root.find("get_data_info")
    size_range = int(get_data_info_tree.get("size_range"))
    power = int(get_data_info_tree.get("power"))
    error_level = float(get_data_info_tree.get("error_level"))
    number_of_points = int(get_data_info_tree.get("number_of_points"))

    #get random coefficeints for function
    if (size_range == 1):
        coeff = -1 + 2 * np.random.random(size = power+1)
    else:
        coeff = np.random.randint(-size_range, size_range, size = power + 1)

    #This is an override of the randomly chosen ceofficents fir the example
    if example:
        coeff = [0, 0.15, -0.65, 0.5]
        error_level = 0.005

    #get random x-values
    x_values = np.random.rand(number_of_points, 1)

    #calculate y_values that correspond to random x_values
    y_values = tools.evaluate_function(coeff, x_values)

    #add some noise to the y_values
    y_values = add_noise(y_values, error_level)

    if to_graph:
        #plot the real function vs the noisy data points, for show
        x = np.linspace(0,1)

        real_y = tools.evaluate_function(coeff, x)

        plt.plot(x_values, y_values, 'o', label = 'Random data')
        plt.plot(x, real_y, '--', label = 'True Polynomial')
        plt.xlabel('x')
        plt.ylabel('y')
        # Create legend
        plt.legend(loc='upper right', prop={'size':12})
        plt.show()

        plt.close()

    #combine to one array
    x_y = np.zeros((len(x_values), 2))
    for i in range(len(x_values)):
        x_y[i,0] = x_values[i]
        x_y[i,1] = y_values[i]

    #Save file for use in fitting
    np.savetxt(input_file_name, x_y, delimiter=",")
    np.savetxt(coeff_input_file_name, coeff, delimiter=",")
    np.savetxt(error_file_name, np.array([error_level]), delimiter=",") #save error level because of hard coded example
    print("Data saved to ", input_file_name)
    print("Coefficents saved to ", coeff_input_file_name)
    print("coefficents = ", coeff)
    print("error level saved to ", error_file_name)
