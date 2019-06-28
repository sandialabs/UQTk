#!/usr/bin/env python
# =====================================================================================
#                      The UQ Toolkit (UQTk) version 3.0.4
#                      Copyright (2017) Sandia Corporation
#                      http://www.sandia.gov/UQToolkit/
#
#      Copyright (2017) Sandia Corporation. Under the terms of Contract DE-AC04-94AL85000
#      with Sandia Corporation, the U.S. Government retains certain rights in this software.
#
#      This file is part of The UQ Toolkit (UQTk)
#
#      UQTk is free software: you can redistribute it and/or modify
#      it under the terms of the GNU Lesser General Public License as published by
#      the Free Software Foundation, either version 3 of the License, or
#      (at your option) any later version.
#
#      UQTk is distributed in the hope that it will be useful,
#      but WITHOUT ANY WARRANTY; without even the implied warranty of
#      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#      GNU Lesser General Public License for more details.
#
#      You should have received a copy of the GNU Lesser General Public License
#      along with UQTk.  If not, see <http://www.gnu.org/licenses/>.
#
#      Need help with UQTk? Check out the resources on http://www.sandia.gov/UQToolkit/
#      or e-mail uqtk-users@software.sandia.gov
#      (subscription details listed at http://www.sandia.gov/UQToolkit/)
#      Other questions? Contact Bert Debusschere <bjdebus@sandia.gov>
#      Sandia National Laboratories, Livermore, CA, USA
# =====================================================================================
# Rosenblatt demonstration and testing script

# Map a Gaussian distribution to a target distribution based on samples from this distribution
# and then evaluate how well the agreement is.

# Make the below described functionality an option.
# Script to read in a set of multi-D Gaussian samples and use the inverse Rosenblatt transformation
# to convert them to an arbitrary target distribution, as defined by a set of samples from this arbitrary
# target distribution.


#
# Reads samples from file arb_rvs.dat, with one n-dimensional RV sample per line
# Numerical parameters to set
#	n_test_samples = 500 # number of test samples to generate for conversion with inverse Rosenblatt
#   ndim = min(3,nvars)  # number of PCA modes to convert

import math
import sys
import getopt
import matplotlib.pyplot as plt

try:
	import numpy as np
except ImportError:
	print "Numpy module could not be found"

try:
	import scipy.stats as sts
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
	import RosenblattUtils
except ImportError:
	print "RosenblattUtils library module not found"


################################################################################
def map_normal_to_beta(n_dim=3,target_alpha=1,target_beta=4,n_target_samples=1000,n_std_samples=1000,n_pdf_points=200,pdf_range=4,verbose=0):

	# Draw samples from a target Beta distribution
	target_samples = np.random.beta(target_alpha,target_beta,size=(n_target_samples,n_dim))
	np.savetxt("TargetSamples.dat",target_samples)

	if verbose > 0:
		print "Generated",target_samples.shape[0],"samples of the target distribution in",target_samples.shape[1],"dimensions."

	# Draw samples from a standard normal distribution to map to the target distribution
	std_samples = np.random.normal(0.0,1.0,size=(n_std_samples,n_dim))
	np.savetxt("StdSamples.dat",std_samples)

	if verbose > 0:
		print "Generated",std_samples.shape[0],"samples of the standard distribution in",std_samples.shape[1],"dimensions."


	# Use the inverse Rosenblatt map to convert them to the distribution of the KL samples
	# TODO: replace this with conversion to uniform and then mapping via uniform to arb function
	mapped_samples = RosenblattUtils.Gauss2KLRVs(target_samples,std_samples,verbose)

	# Write out mapped samples to file
	mapped_sample_file_name = "StdMappedToTarget.dat"
	np.savetxt(mapped_sample_file_name,mapped_samples)

	# Compute the mean and standard deviations along each dimension
	pdf_target_means = np.mean(target_samples,0)
	pdf_target_stdvs = np.std(target_samples,0)

	# Compute the PDF of the target distribution along each dimension
	# Set up the x-locations to evaluate the pdf at
	pdf_target_x_locs = np.zeros((n_pdf_points,n_dim)) # Locations to evaluate the PDF in
	for i_dim in range(n_dim):
		pdf_target_x_locs[:,i_dim] = np.linspace(pdf_target_means[i_dim] - pdf_range*pdf_target_stdvs[i_dim],
		                                        pdf_target_means[i_dim] + pdf_range*pdf_target_stdvs[i_dim], n_pdf_points)

	# Evaluate PDFs on those locations
	pdf_target = RosenblattUtils.compute_pdf_kde(target_samples,pdf_target_x_locs)
	pdf_mapped = RosenblattUtils.compute_pdf_kde(mapped_samples,pdf_target_x_locs)

	# Plot the marginal PDFs in each dimension
	if verbose > 0:
		print "Plotting marginal pdfs in each dimension"
		RosenblattUtils.plot_pdf_comparison_per_dim(pdf_target_x_locs,pdf_target,"Target",pdf_mapped,"Mapped")

	# Compute Kullback-Leibler divergence between the distributions corresponding to the
	# input Target samples and the samples that have been mapped to the target distribution.
	# Do this for each dimension individually

	kl_divs = np.zeros(n_dim)
	for i_dim in range(n_dim):
		kl_divs[i_dim] = sts.entropy(pdf_target[:,i_dim],pdf_mapped[:,i_dim])

	kl_divs_ave = np.mean(kl_divs)

	if verbose > 0:
		print '\nKullback-Leibler Divergences between target and mapped distributions per dimension:'
		for i_dim in range(n_dim):
			print "Dimension",i_dim,": ","{0:10.1e}".format(kl_divs[i_dim])
		print "Average     : "          ,"{0:10.1e}".format(kl_divs_ave)

	# Compute and compare statistical moments
	stats_error_norm = RosenblattUtils.compare_sample_stats(target_samples,mapped_samples,verbose)

	if verbose > 0:
		print '\nNorm of difference between moments in target and mapped distribution:'
		print '(Based on mean, variance, skewness, kurtosis, '
		print 'and normalized by number of dimensions and the number of moments)'
		print 'Moments Error:','{0:10.3e}'.format(stats_error_norm)

	return kl_divs_ave, stats_error_norm
################################################################################
help_string = """
Usage:
  RosenblattDemo.py [-h] [-d <dimension>][--nt <# target samples>] [--ns <# std samples>] [--nr <# replicas>] [--verbose]
What:
  Demonstrate / Test Rosenblatt transformation by mapping known distributions into each other and
  computing metrics of the differences between the target and the mapped distributions
Where
  -h = print help info
  -d = dimension to use [defaults to 3]
  --nt = number of target distribution samples to use [defaults to 1000]
  --ns = number of standard distribution samples (to be mapped to target) to use [defaults to 1000]
  --nr = number of replica tests to run [defaults to 1]
  --verbose = to switch on more verbose output [defaults to not verbose]
"""

if __name__ == "__main__":
    #
    # Process inputs
    #
    try:
        opts,extra_arguments = getopt.getopt(sys.argv[1:],"hd:",["nt=","ns=","nr=","verbose"])
    except getopt.GetoptError, err:
        print str(err)
        print help_string
        sys.exit(1)

    # Default values
    # Script parameters
    demo_verbose = 0 # set > 0 to get more verbose output (and higher numbers generate even more output)
    n_pdf_points = 200 # number of points to evalute RV PDFs at in order to plot them or to compute Kullback-Leibler divergences
    pdf_range = 4.0 # Sample from  (mean - 4*stdv to mean + 4*stdv) to cover the support of the full PDF.
                    # (If PDF is higlhy skewed, it might be better to sample from min to max)
    n_replicas = 1  # Number of replicas to run of the tests

    # Target distribution parameters
    n_dim = 3 # Dimensionality of the target distribution
    n_target_samples = 1000 # number of samples to take from the target distribution
    target_alpha = 1 # Beta distribution shape parameter alpha
    target_beta = 4  # Beta distribution shape parameter beta

    # Standard distribution parameters
    n_std_samples = 1000 # number of samples from standard distribution to map to the target distribution

    for o,a in opts:
      if o == "-h":
        print help_string
        sys.exit(0)
      elif o == "-d":
        n_dim = int(a)
      elif o == "--nt":
        n_target_samples = int(a)
      elif o == "--ns":
        n_std_samples = int(a)
      elif o == "--nr":
        n_replicas = int(a)
      elif o == "--verbose":
        demo_verbose = 1
      else:
        assert False, "Unhandled command line parsing option. Use -h flag to get usage info."

    # error checking
    if (n_target_samples < 0):
      print "The number of target distribution samples needs to be >= 0"
      print help_string
      sys.exit(1)

    if (n_std_samples < 0):
      print "The number of standard distribution samples needs to be >= 0"
      print help_string
      sys.exit(1)

    if (n_dim < 0):
      print "The dimensionality needs to be >= 0"
      print help_string
      sys.exit(1)

    if (n_replicas < 0):
      print "The number of replicas needs to be >= 0"
      print help_string
      sys.exit(1)

    # Run Rosenblatt demo
    for i_test in range(n_replicas):
        if demo_verbose > 0:
            print "Testing replica #: ",i_test
        # Run the tests
        kl_divs_ave, stats_error_norm = map_normal_to_beta(n_dim,target_alpha,target_beta,n_target_samples,n_std_samples,n_pdf_points,pdf_range,demo_verbose)
        # Output the error metrics
        print kl_divs_ave, stats_error_norm
