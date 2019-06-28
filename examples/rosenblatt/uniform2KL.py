#!/usr/bin/env python
# Script to read in a set of multi-D uniform samples (on [0,1]^D) and use the Rosenblatt transformation
# to convert them to an arbitrary distribution, as defined by a set of random variable samples.
# Notes, this script is for a particular case where the samples first need to be divided by the square root
# of the number of samples to get unit variance.

#
# Reads samples from file arb_rvs.dat, with one n-dimensional RV sample per line
# Numerical parameters to set
#	n_test_samples = 500 # number of test samples to generate for conversion with inverse Rosenblatt
#   ndim = min(3,nvars)  # number of PCA modes to convert

import math

try:
	import numpy as np
except ImportError:
	print "Numpy module could not be found"

try:
	from scipy import stats
except ImportError:
	print "Scipy stats module could not be found"

try:
	import PyUQTk.uqtkarray as uqtkarray
	import PyUQTk.quad as uqtkquad
	import PyUQTk.pce as uqtkpce
	import PyUQTk.tools as uqtktools
except ImportError:
	print "PyUQTk array, quad, PCE, or tools module not found"

try:
	import data_utils
except ImportError:
	print "data_utils library module not found"

#
# Main Script
#

# Script parameters
verbose = 1 # set > 0 to get more verbose output (and higher numbers generate even more output)
n_test_samples = 1000 # number of test samples to generate from PCEs
n_modes = 7 # only consider first n_modes dimensions


# Read in data file with one multi-D sample per row
rvs_infile = np.genfromtxt("arb_rvs.dat")
nvars = rvs_infile.shape[1] # Number of KL variables sampled in file
nsamp = rvs_infile.shape[0] # Number of samples in the input file

# Use only subset of dimensions
ndim = min(n_modes,nvars)
print "Considering",ndim,"modes"
# Store the samples we want to keep, and apply normalizing factor N^{1/2}
rvs_in = rvs_infile[:,0:ndim]*np.sqrt(nsamp)
# Write rescaled input samples to file
np.savetxt("arb_rvs_unitvar.dat",rvs_in)


#
# Re-create samples from the arbitrary RV distribution using the inverse Rosenblatt transformation,
# starting from uniformly distributed samples.
#
unif_samples = np.random.uniform(0.0,1.0,size=(n_test_samples, ndim))
np.savetxt("UnifSamples.dat",unif_samples)

if verbose > 1:
	print "Shape of numpy array with uniform samples:",unif_samples.shape

# Use the inverse Rosenblatt map to convert them to the distribution of the Arbitrary samples
rvs_mapped = data_utils.Unif2ArbSamples(rvs_in,unif_samples,verbose)


# Write out mapped samples to file
sample_file_name = "InvRosRVSamplesFromUnif.dat"
np.savetxt(sample_file_name,rvs_mapped)
