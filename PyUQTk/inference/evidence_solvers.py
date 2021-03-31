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

#=====================================================================================
#     For more info, see also the UQTk manual section on "Bayesian Evidence Estimation"
#=====================================================================================

# Import python libraries.
import numpy as npy
import os
import scipy.stats as stats
import subprocess



################################################################################
# Function: LikelihoodMC_PriorSamples
#
# Inputs:
# ln_likelihood --- vector of ln-likelihood values corresponding to prior samples
#
# Outputs:
# ln-evidence estimate
################################################################################

def LikelihoodMC_PriorSamples(ln_likelihood):

    # Performs Monte Carlo marginalization of the likelihood using
    # prior samples.
    running_sum = 0.0
    n_samples = ln_likelihood.shape[0]
    for i in range(n_samples):
        running_sum = running_sum + npy.exp(ln_likelihood[i])

    return npy.log(running_sum / n_samples)



################################################################################
# Function: ImportanceLikelihoodMC_PosteriorSamples
#
# The function works in two stages. The first stage involves
# constructing the biasing distribution and generating samples from that
# distribution.
#
# Stage 1 inputs:
# posterior_samples --- array of posterior samples (each row is a sample)
# n_importance_samples --- number samples requested from the biasing
#                          distribution
# stage --- set to 1 for stage 1
#
# Stage 1 outputs:
# importance_samples --- array of samples from the biasing distribution
#                        (each row is a sample)
# importance_samples_ln_PDF --- vector of ln-PDF values of the biasing
#                               distribution  corresponding to these samples
#
# At this point, the user needs to externally compute and provide the
# ln-prior and ln-likelihood values for these samples and pass them back
# into the function. The second stage can then estimate the ln-evidence.
#
# Stage 2 inputs:
# ln_prior --- vector of ln-prior values corresponding to the biasing
#              samples generated in stage 1
#
# ln_likelihood --- vector of ln-likelihood values corresponding to the
#                   biasing samples generated in stage 1
# ln_importance_input --- pass back in the output importance_samples_ln_PDF
#                         generated from stage 1 without modifications
# stage --- set to 2 for stage 2
#
# Stage 2 outputs:
# ln-evidence estimate
################################################################################

def ImportanceLikelihoodMC_PosteriorSamples(posterior_samples, n_importance_samples, ln_prior, ln_likelihood, ln_importance_input, stage):

    # Performs Monte Carlo marginalization of the likelihood using
    # importance sampling. Here the importance distribution is the
    # Gaussian approximation to the posterior constructed from
    # posterior sample moments. Note this is in two stages. Stage 1:
    # construct importance distribution and generate importance
    # samples. Stage 2: estimate the ln-evidence once the
    # ln-likelihood and ln-prior are computed (and thus the importance
    # weights).
    if (stage == 1):
        posterior_sample_mean = npy.mean(posterior_samples, axis = 0)
        posterior_sample_cov = npy.cov(posterior_samples, rowvar = False)
        if (posterior_samples.shape[1] == 1):
            posterior_sample_cov = posterior_sample_cov.reshape(1, 1)

        # Generates samples from importance distribution.
        importance_samples = npy.zeros((n_importance_samples, posterior_samples.shape[1]))
        importance_samples_ln_PDF = npy.zeros(n_importance_samples)
        for i in range(n_importance_samples):
            importance_samples[i, :] = npy.random.multivariate_normal(posterior_sample_mean, posterior_sample_cov.reshape(posterior_samples.shape[1], posterior_samples.shape[1]))
        # Computes the ln-PDF values of importance sample in bulk.
        importance_samples_ln_PDF = stats.multivariate_normal.logpdf(importance_samples, posterior_sample_mean, posterior_sample_cov)

        return importance_samples, importance_samples_ln_PDF

    elif (stage == 2):
        # Computes the importance MC estimate.
        running_sum = 0.0
        n_samples = ln_likelihood.shape[0]
        # Need to make this running sum robust against overflow by updating the running average instead
        for i in range(n_samples):
            running_sum = running_sum + npy.exp(ln_likelihood[i] + ln_prior[i] - ln_importance_input[i])

        return npy.log(running_sum / n_samples)
    else:
        print('Error: ImportanceLikelihoodMC_PosteriorSamples misuse.')
        exit


################################################################################
# Function: PosteriorGaussian_PosteriorSamples
#
# Inputs:
# posterior_samples --- array of posterior samples (each row is a sample)
# ln_prior --- vector of ln-prior values corresponding to the posterior samples
# ln_likelihood --- vector of ln-likelihood values corresponding to the posterior
#                   samples
#
# Outputs:
# ln-evidence estimate
################################################################################

def PosteriorGaussian_PosteriorSamples(posterior_samples, ln_prior, ln_likelihood):

    # Computes the ln-posterior from a Gaussian approximation using
    # posterior sample moments.
    posterior_sample_mean = npy.mean(posterior_samples, axis = 0)
    posterior_sample_cov = npy.cov(posterior_samples, rowvar = False)

    ln_posterior = stats.multivariate_normal.logpdf(posterior_samples, posterior_sample_mean, posterior_sample_cov)

    # Computes the resulting ln-evidence.
    ln_evidence_values = ln_prior + ln_likelihood - ln_posterior

    # Extracts the mean of the ln-evidence estimates from all
    # posterior samples (in a perfect world, they should all be the
    # same value).
    return npy.mean(ln_evidence_values)



################################################################################
# Function: Harmonic_PosteriorSamples
#
# Inputs:
# ln_likelihood --- vector of ln-likelihood values corresponding to the
#                   posterior samples
#
# Outputs:
# ln-evidence estimate
################################################################################

def Harmonic_PosteriorSamples(ln_likelihood):

    # Performs Harmonic estimate using posterior samples.
    running_sum = 0.0
    n_samples = ln_likelihood.shape[0]
    for i in range(n_samples):
        running_sum = running_sum + 1.0 / npy.exp(ln_likelihood[i])

    return -npy.log(running_sum / n_samples)
