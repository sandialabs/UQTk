#!/usr/bin/env python
#=====================================================================================
#
#                      The UQ Toolkit (UQTk) version 3.1.0
#                          Copyright (2020) NTESS
#                        https://www.sandia.gov/UQToolkit/
#                        https://github.com/sandialabs/UQTk
#
#     Copyright 2020 National Technology & Engineering Solutions of Sandia, LLC (NTESS).
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

# Testing evidence solver with data from awsm_fit
# finds evidence 3 different ways: harmonic, Gaussian prosterior, and Importance sampling
# prints out evidence that is found in a varity of ways
# functions that calculate evidence were made by Xun


from __future__ import print_function # To make print() in Python 2 behave like in Python 3

import numpy as np
import sys
import os
import getopt
import xml.etree.ElementTree as ET # for the xml tree and parsing
import unittest

try:

    import PyUQTk.inference.postproc as uqtkinfpp
    import PyUQTk.inference.evidence_solvers as evd
except ImportError:
    print("Could not import the UQTk postprocessing tools")
    sys.exit(1)

#import awsm_tools, needed to calculate log_like in Importance sampleing
try:
    import tools
except ImportError:
    print("Could not import the tools module")
    sys.exit(1)



################################################################################
def Importance(posterior_samples, n_importance_samples, ln_prior, param_names_dict, model):
    '''
    Performs all stages of the Importance Sampling method to find evidence

    inputs:
     posterior_samples: an array with samples in different rows and parameters on the columns
     n_importance_samples: an integer for how many importance samples you want to take
     ln_prior: a number (not array) and is used as a uniform prior
     model_name: name of the model to use
     param_names_dict: the result of model_info.get("param_names_dict"), where model_info is the output of awsm_tools.get_model_info
     *also need an input.xml file to get some values for the postAWSM, this file is the same as used for awsm_fit.py

    output:
     return value is the evidence
     '''

    #importance likelihood stage 1
    #gets the importance samples and the log of the pdf
    importance_samples, importance_samples_ln_PDF = evd.ImportanceLikelihoodMC_PosteriorSamples(posterior_samples, n_importance_samples, [], [], [], 1)


    #importance likelihood stage 1.5
    # taken from e-mail from Xun
    # Computes the ln_likelihood of the importance samples.
    importance_samples_ln_likelihood = np.zeros(n_importance_samples)


    for i in range(n_importance_samples):
        # Computes the ln-likelihood.
        parameter_sample = importance_samples[i]
        #handle when sigma is less than zero and causes error when ln is taken
        # sig_index = param_names_dict['sigma']
        # if(parameter_sample[sig_index] < 0):
        #     #### not exactly sure if this is the right thing to do
        #     #importance_samples_ln_likelihood[i] = 0.0
        #     # We should really be doing this via the prior, since the prior on \sigma should exclude these values
        #     importance_samples_ln_likelihood[i] = - float("inf")
        # else:
            # Compute ln-likelihood from postAWSM (since this routine currently does not
            # factors in the prior, as it is uniform anyway. So the posterior probabilities returned
            # are actually likelihood values.)
        [importance_samples_ln_likelihood[i], ln_prior] = model.postAWSM(parameter_sample, {})
        #importance_samples_ln_likelihood[i] = awsm_tools.postAWSM(parameter_sample,postinfo)

    # Construct the prior array of correct size. Will all be same value because of uniform prior
    # (Note, this is only valid if the sampled value actually lies inside the valid prior range!
    # This can not be assumed for the values sampled from the Gaussian approximation. NEED TO FIX THIS!)
    importance_samples_ln_prior = np.array([ln_prior] * n_importance_samples)

    #stage 2, calculates evidence
    Importance_evidence = evd.ImportanceLikelihoodMC_PosteriorSamples([], [], importance_samples_ln_prior, importance_samples_ln_likelihood, importance_samples_ln_PDF, 2)

    return Importance_evidence
################################################################################
#making class for unit testing
class testing_evidence(unittest.TestCase):

    def test_True(self):
        #This just tests the testing method
        self.assertTrue(True)

#functions that will be made into methods
def test_almost_equal(a, b):
    def test(self):
        self.assertAlmostEqual(a, b)
    return test

def test_one_percent(actual, calculated):
    def test(self):
        one_percent = actual * 0.01
        self.assertGreaterEqual(actual + one_percent, calculated)
        self.assertLessEqual(actual - one_percent, calculated)
    return test

################################################################################

def compute_evidence_single_model(model, all_samples, v_names): #(my_model, opts, all_samples, v_names):
    '''
      Calculated the evidence of an MCMC chain in two different ways:
        Gaussian, and Importance Sampling
        Also writes it to an output file
        Calculates for a single model

    Inputs:
        my_model: which model to calculate the evidence for
        opts: the opts variable from awsm_post.py
        all_samples, v_names: the output of uqtkinfpp.extract_all_vars, determined in awsm_post.py

    Output:
        Returns the evidence values to the command line and output file.
        Also preforms some simple testing of evidence values
    '''

    # Defaults
    n_importance_samples = 100000
    model_info  = model.model_info


    # Create an easier to use list of parameter names that is sorted according to the same
    # index used in all other species arrays.
    param_names_dict = model.model_info["param_names_dict"]
    param_names_list = sorted(param_names_dict.keys(), key=param_names_dict.__getitem__)

    # # open file to append to
    # fd = open(write_file_name, 'a')

    ln_likelihood = all_samples[:,-1]
    posterior_samples = all_samples[:, 1:len(v_names) + 1]

    #get ln_prior vector
    # Get ranges of inferred variables from model_info
    spllo = model_info.get("spllo") # Low limits on sampling ranges for parameters
    splhi = model_info.get("splhi") # High limits on sampling ranges for parameters

    # Calculate log prior
    ln_prior = 0.0
    for i in range(len(spllo)):
        ln_prior -= np.log(splhi[i] - spllo[i])

    # Put log prior in array of same size as the log likelihood
    ln_prior_array = np.array([ln_prior] * len(ln_likelihood))


    #compute the evidence using posterior Gaussian
    gaussian_evidence = evd.PosteriorGaussian_PosteriorSamples(posterior_samples, ln_prior_array, ln_likelihood)


    #compute evidence using Importance sampling
    importance_evidence = Importance(posterior_samples, n_importance_samples, ln_prior, param_names_dict, model)

    print("Importance Sampling: ", importance_evidence)

    return (gaussian_evidence, importance_evidence)
