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

import sys
import numpy as np

# File manipulation utilities:

def extract_vars(samples_file_name,n_burnin,v_names,debug,stride=1):
    """From a file with samples in ascii format, with
    the first line containing the label for each column, extract
    the columns with the labels in v_names and return them
    in a numpy array. Remove n_burnin samples from the top.
    Only read one in every stride number of lines after that."""

    # Open text file with all samples,
    samples_file = open(samples_file_name,"r")

    #sample_lines = samples_file.readlines()
    #samples_file.close()

    # Extract first line with the column labels and find the column
    # numbers corresponding to the variables of interest.
    #labels_line = sample_lines[0].rstrip('\n')
    labels_line = samples_file.readline().rstrip('\n')
    col_labels = [lbl for lbl in labels_line.split()]

    v_indices = []
    for s_v in v_names:
        try:
            i_v = col_labels.index(s_v)
            v_indices.append(i_v)
        except ValueError:
            print("Variable", s_v, "is not found in the list of labels", col_labels)
            sys.exit(1)

    if (debug > 0):
        print("Column labels in file",samples_file_name,"are:",col_labels)
        for i_v in range(len(v_names)):
            print("The column number of",v_names[i_v],"is:",v_indices[i_v])

    # Read subsequent lines, leaving out the first n_burnin, and only one
    # in every stride lines after that
    samples_list = []
    line_no = 0
    done = 0
    while not done:
        line = samples_file.readline()
        if (line == ""):
            done = 1
        else:
            line_no += 1
            if (line_no > n_burnin and (line_no - n_burnin) % stride == 0):
                records = line.split()
                num_records = [float(s) for s in records]
                samples_list.append(num_records)

    # Close the file
    samples_file.close()

    # Remove the last line if is has a value < 0 (i.e. -1) in the acceptance_prob column
    try:
        i_ap = col_labels.index("acceptance_prob")
    except ValueError:
        i_ap = -1

    # If this is a file with acceptance probabilities
    if (i_ap >= 0):
        # And the last line has a negative acceptance probability
        if(samples_list[-1][i_ap] < 0):
            # Remove the last line
            del samples_list[-1]
            if (debug > 0):
                print("The last sample line has been deleted as it contained the MAP values")

    # Convert list to array
    steady_samples = np.array(samples_list)

    # Remove burn-in samples from the top
    #if (n_burnin > 0):
    #    steady_samples = all_samples[n_burnin:,:]
    #else:
    #    steady_samples = all_samples

    #if (debug > 0):
    #    print "Removed", n_burnin, "burn-in samples"

    # Extract all columns of interest
    samples_cols = []
    for i_v in v_indices:
        samples_cols.append(steady_samples[:,i_v])

    samples = np.array(samples_cols).T
    if (debug > 0):
        print("Shape of samples array:",samples.shape)

    n_samples = len(samples[:,0])
    n_vars = len(samples[0,:])

    if (debug > 0):
        print("Read in", n_samples, "regular samples of", n_vars, "variables from file", samples_file_name)

    return samples

def extract_all_vars(samples_file_name,n_burnin,debug,stride=1):
    """From a file with samples in ascii format, with the first line containing the label
    for each column, extract all variables and return them in a numpy array. Remove
    n_burnin samples from the top. Only read one in every stride number of lines after that."""

    # Open text file with all samples
    samples_file = open(samples_file_name,"r")
    #sample_lines = samples_file.readlines()
    #samples_file.close()

    # Extract first line with the column labels and find the column
    # numbers corresponding to the variables of interest.
    #labels_line = sample_lines[0].rstrip('\n')
    labels_line = samples_file.readline().rstrip('\n')
    col_labels = [lbl for lbl in labels_line.split()]

    # Identify the MCMC vars, knowing that the first column is the step
    # number and the last two columns are acceptance and posterior prob
    n_cols = len(col_labels)
    n_vars = n_cols - 3

    v_names = col_labels[1:1+n_vars]

    if (debug > 0):
        print("Column labels in file", samples_file_name, "are:", col_labels)
        print("MCMC chain variables are", v_names)

    # Read subsequent lines, leaving out the first n_burnin, and only one
    # in every stride lines after that
    samples_list = []
    line_no = 0
    done = 0
    while not done:
        line = samples_file.readline()
        if (line == ""):
            done = 1
        else:
            line_no += 1
            if (line_no > n_burnin and (line_no - n_burnin) % stride == 0):
                records = line.split()
                num_records = [float(s) for s in records]
                samples_list.append(num_records)

    # Close the file
    samples_file.close()

    # Remove the last line if is has a value < 0 (i.e. -1) in the acceptance_prob column
    try:
        i_ap = col_labels.index("acceptance_prob")
    except ValueError:
        i_ap = -1

    # If this is a file with acceptance probabilities
    if (i_ap >= 0):
        # And the last line has a negative acceptance probability
        if(samples_list[-1][i_ap] < 0):
            # Remove the last line
            del samples_list[-1]
            if (debug > 0):
                print("The last sample line has been deleted as it contained the MAP values")

    # Convert list to array
    samples = np.array(samples_list)

    # Remove burn-in samples from the top
    #if (n_burnin > 0):
    #    samples = all_samples[n_burnin:,:]
    #else:
    #    samples = all_samples

    #if (debug > 0):
    #    print "Removed", n_burnin, "burn-in samples"

    n_samples = len(samples[:,0])

    if (debug > 0):
        print("Read in", n_samples, "regular samples of", n_vars, "variables from file", samples_file_name)

    return samples, v_names


