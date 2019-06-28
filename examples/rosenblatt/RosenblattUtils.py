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
#
# Import libraries
#
import math
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


################################################################################
# Function definitions
################################################################################
def Unif2ArbSamples(rv_arb_in,rv_unif_in,verbose=0):
	"""Given a set of original RV samples (rv_samples_in) from an arbitrary distribution,
	build a Rosenblatt map and use it to convert a set of uniformly distributed
	input samples (rv_samples_in) into a set of samples distributed according
	to the original set of RV samples.
	Input:
		rv_arb_in  : numpy array with input RV samples from arbitrary distribution.
					 Each line is a sample. The columns
		        	 represent the dimensions of the input RV.
	    rv_unif_in : numpy array with uniformly distributed samples (on [0,1]) to be mapped to the
		             distribution of rv_arb_in. Each line is a sample. The columns
					 represent the dimensions of the input RV and
					 should match the ones in rv_arb_in.
		verbose    : verbosity level (more output for higher values).

	Returns        : numpy array with mapped samples. Same dimensionality as rv_unif_in.

	NOTE:
		* Maybe rename to InverseRosenblatt
		* Do a test that the input std samples are between 0 and 1 (is their std dev also 1 ??)

	"""

	# Algorithm parameters
	bw = -1 # KDE bandwidth for Rosenblatt (on interval 0.1)
	iiout = 500 # interval for output to screen

	# Dimensionality and number of samples of input RVs
	ndim_arb_in = rv_arb_in.shape[1]
	nsamp_arb_in = rv_arb_in.shape[0]

	# Set up transpose of input data for the inverse Rosenblatt transformation in a UQTk array
	# Rosenblatt functions want input with dimensions of [ndim,nsamp]
	rv_arb_in_uqtk = uqtkarray.numpy2uqtk(rv_arb_in.T)

	# Dimensionality and number of uniform samples to be converted to the KL RV
	# distribution using the inverse Rosenblatt map
	ndim_in = rv_unif_in.shape[1]
	nsamp_in = rv_unif_in.shape[0]

	# check that dimensionality matches the KL RV dimensionality
	if ndim_in != ndim_arb_in:
		print "Unif2ArbSamples: dimensionality of input uniform RVs does not match dimensionality of target RVs!"
		quit(1)

	# Set up numpy array for mapped sample points
	mappedRVs_np = np.zeros((nsamp_in,ndim_in))

	# Map all samples in the uniform input to the distribution given by the data
	# using the inverse Rosenblatt transformation
	for ipt in range(nsamp_in):

		# print "Converting sample point #",ipt

		# Set up working arrays
		unif_uqtk_1s = uqtkarray.dblArray1D(ndim_in,0.0)
		mappedRV_uqtk_1s = uqtkarray.dblArray1D(ndim_in,0.0)

		# Store current sample on [0,1] in work array
		for idim in range(ndim_in):
			unif_uqtk_1s[idim] = rv_unif_in[ipt,idim]

		# Map each point from uniform[0,1] to the distribution given by the original samples via inverse Rosenblatt
		if bw > 0:
			uqtktools.invRos(unif_uqtk_1s,rv_arb_in_uqtk,mappedRV_uqtk_1s,bw)
		else:
			uqtktools.invRos(unif_uqtk_1s,rv_arb_in_uqtk,mappedRV_uqtk_1s)

		# Store results
		for idim in range(ndim_in):
			mappedRVs_np[ipt,idim] = mappedRV_uqtk_1s[idim]

		# Screen diagnostic output
		if verbose > 1:
			if ((ipt+1)%iiout == 0) or ipt==0 or (ipt+1)==nsamp_in:
				print "Inverse Rosenblatt of uniform input samples:",(ipt+1),"/",nsamp_in,"=",(ipt+1)*100/nsamp_in,"% completed"

	# Return samples mapped to target distribution
	return mappedRVs_np

################################################################################
def Gauss2KLRVs(kl_rvs_in,gauss_rvs_in,verbose=0):
	"""Given a set of original KL RV samples (kl_rvs_in), build a Rosenblatt map and use it
	to convert a set of normally distributed input samples (gauss_rvs_in) into
	a set of samples distributed according to the original set of KL RV samples.
	Input:
		kl_rvs_in	: numpy array with input RV samples. Each line is a sample. The columns
		              represent the dimensions of the input RV
	    gauss_rvs_in: numpy array with normally distributed samples to be mapped to the
		              distribution of kl_rvs_in. Each line is a sample. The columns
					  represent the dimensions of the input RV. Number of columns
					  should match the ones in kl_rvs_in
		verbose     : verbosity level (more output for higher values)

	Returns         : numpy array with mapped samples. Same dimensionality as gauss_rvs_in
	"""

	# Algorithm parameters
	bw = -1 # KDE bandwidth for Rosenblatt (on interval 0.1)
	iiout = 500 # interval for output to screen

	# Dimensionality and number of samples of input RVs
	ndim_kl = kl_rvs_in.shape[1]
	nsamp_kl = kl_rvs_in.shape[0]

	# Set up transpose of input data for the inverse Rosenblatt transformation in a UQTk array
	# Rosenblatt functions want input with dimensions of [ndim,nsamp]
	kl_samp_uqtk = uqtkarray.numpy2uqtk(kl_rvs_in.T)

	# Dimensionality and number of Gaussian samples to be converted to the KL RV
	# distribution using the inverse Rosenblatt map
	ndim_in = gauss_rvs_in.shape[1]
	nsamp_in = gauss_rvs_in.shape[0]

	# check that dimensionality matches the KL RV dimensionality
	if ndim_in != ndim_kl:
		print "Gauss2KLRVs: dimensionality of input Gaussian RVs does not match KL RV target dimensionality"
		quit(1)

	# Set up numpy array for mapped sample points
	mappedRVs_np = np.zeros((nsamp_in,ndim_in))

	# Map all samples in Gaussian input to the distribution given by the data
	# using the inverse Rosenblatt transformation
	for ipt in range(nsamp_in):

		# print "Converting sample point #",ipt

		# Set up working arrays
		unif_uqtk_1s = uqtkarray.dblArray1D(ndim_in,0.0)
		mappedRV_uqtk_1s = uqtkarray.dblArray1D(ndim_in,0.0)

		# First map each point to uniform[0,1]
		# PCtoPC maps to [-1,1], which then gets remapped to [0,1]
		for idim in range(ndim_in):
			unif_uqtk_1s[idim] = (uqtktools.PCtoPC(gauss_rvs_in[ipt,idim],"HG",0.0,0.0,"LU",0.0,0.0)+1.0)/2.0

		# Map each point from uniform[0,1] to the distribution given by the original samples via inverse Rosenblatt
		if bw > 0:
			uqtktools.invRos(unif_uqtk_1s,kl_samp_uqtk,mappedRV_uqtk_1s,bw)
		else:
			uqtktools.invRos(unif_uqtk_1s,kl_samp_uqtk,mappedRV_uqtk_1s)

		# Store results
		for idim in range(ndim_in):
			mappedRVs_np[ipt,idim] = mappedRV_uqtk_1s[idim]

		# Screen diagnostic output
		if verbose > 1:
			if ((ipt+1)%iiout == 0) or ipt==0 or (ipt+1)==nsamp_in:
				print "Inverse Rosenblatt of Gaussian input samples:",(ipt+1),"/",nsamp_in,"=",(ipt+1)*100/nsamp_in,"% completed"

	# Return samples mapped to target distribution
	return mappedRVs_np
################################################################################
def compute_pdf_kde(dist_samples,x_locs):
	"""Compute the PDF of a distribution given by samples by relying
	on Kernel Density Estimation along each dimension of the distribution.
	Input:
		dist_samples: samples of the distribution to compute the PDF of. Each row is one multi-D sample
		x_locs      : locations to evaluate PDF at in each dimension. Each column has the locations for
		              the corresponding dimension
	Return:
		pdf_eval    : 1D marginal PDF evaluated at the x-locs in each dimension. Each column has the
					  evaluated marginal PDF for the corresponding dimension.
	"""
	# Get the number of points and the dimension.
	n_pdf_points = x_locs.shape[0]
	n_dim = x_locs.shape[1]

	# Set up data structure for evaluated PDFs
	pdf_eval = np.zeros((n_pdf_points,n_dim))

	# Evaluated PDF with KDE
	for i_dim in range(n_dim):
		kern_i = sts.kde.gaussian_kde(dist_samples[:,i_dim])
		pdf_eval[:,i_dim] = kern_i(x_locs[:,i_dim])

	# Return result
	return pdf_eval
################################################################################
def compare_sample_stats(target_samples,mapped_samples,verbose=0):
	"""Compare the statistics between two sets of samples and compute the norm of
	the difference based on the mean, variance, skewness and kurtosis.
	Input:
		target_samples: set of samples of the first distribution. Each row is one multi-D sample.
		mapped_samples: set of samples of the second distribution. Each row is one multi-D sample.
		verbose       : generates more screen output as verbose > 0
	Output:
		Returns the root mean square of the differences in mean, variance,
		skewness and kurtosis between the two data sets.
	"""
	# Get the dimensionality of each sample:
	n_dim = target_samples.shape[1]

	# Compute and compare statistical moments
	target_dist_stats = sts.describe(target_samples, axis=0)
	mapped_dist_stats = sts.describe(mapped_samples, axis=0)

	if verbose > 0:
		print '\nComparison of statistics between target and mapped distributions per dimension:'
		print '                            Target      Mapped'
		print 'Number of data points:','{:10d}'.format(target_dist_stats[0]),'{:10d}'.format(mapped_dist_stats[0])
		print 'Minima:'
		for i_dim in range(n_dim):
		    print '        Dimension:',i_dim,': ','{0:10.1e}'.format(target_dist_stats[1][0][i_dim]), \
			                                      '{0:10.1e}'.format(mapped_dist_stats[1][0][i_dim])
		print 'Maxima:'
		for i_dim in range(n_dim):
		    print '        Dimension:',i_dim,': ','{0:10.1e}'.format(target_dist_stats[1][1][i_dim]), \
			                                      '{0:10.1e}'.format(mapped_dist_stats[1][1][i_dim])
		print 'Mean:'
		for i_dim in range(n_dim):
		    print '        Dimension:',i_dim,': ','{0:10.1e}'.format(target_dist_stats[2][i_dim]), \
			                                      '{0:10.1e}'.format(mapped_dist_stats[2][i_dim])
		print 'Variance:'
		for i_dim in range(n_dim):
		    print '        Dimension:',i_dim,': ','{0:10.1e}'.format(target_dist_stats[3][i_dim]), \
			                                      '{0:10.1e}'.format(mapped_dist_stats[3][i_dim])
		print 'Skewness:'
		for i_dim in range(n_dim):
		    print '        Dimension:',i_dim,': ','{0:10.1e}'.format(target_dist_stats[4][i_dim]), \
			                                      '{0:10.1e}'.format(mapped_dist_stats[4][i_dim])
		print 'Kurtosis:'
		for i_dim in range(n_dim):
		    print '        Dimension:',i_dim,': ','{0:10.1e}'.format(target_dist_stats[5][i_dim]), \
			                                      '{0:10.1e}'.format(mapped_dist_stats[5][i_dim])

	# Compute error norm between target and mapped distributions based on Mean, Variance, Skewness
	# and Kurtosis moments
	stats_error_norm = 0.0
	n_moments = 4
	i_start = 2
	i_end = 2+n_moments
	for i_stat in range(2,2+n_moments):
		stats_error_norm += np.sum(np.square(target_dist_stats[i_stat]-mapped_dist_stats[i_stat]))

	# normalize by number of moments and take square root
	norm_fac = float(n_dim * n_moments)
	stats_error_norm = math.sqrt(stats_error_norm/norm_fac)

	# Return the norm
	return stats_error_norm
################################################################################
def plot_pdf_comparison_per_dim(x_locs,pdf_1,pdf_name_1,pdf_2,pdf_name_2):
	"""Plot a comparison between PDFs, one plot per dimension comparing two
	distributions.
	Input
		x_locs		: Locations at which to plot PDFs [n_locations x n_dimensions]
		pdf_1	    : First set PDF values (same dimension as x_locs)
		pdf_name_1	: Name of the first PDF
		pdf_2       : Second set of PDF values (same dimension as x_locs)
		pdf_name_2  : Name of the second PDF

	Output
		Set of plots stored in files "pdf_comparison_d*.pdf"
	"""
	# Plot parameters
	lw1=2 # Line width
	fs1=18 # Font size

	# Get dimension of data sets
	n_dim = x_locs.shape[1]

	# Make plots for each dimension
	for i_dim in range(n_dim):
		fig = plt.figure(figsize=(8,6))
		ax=fig.add_axes([0.10, 0.10, 0.85, 0.85]) ;
		l1=plt.plot(x_locs[:,i_dim],pdf_1[:,i_dim],linewidth=lw1,label=pdf_name_1+" Dim "+str(i_dim))
		l2=plt.plot(x_locs[:,i_dim],pdf_2[:,i_dim],linewidth=lw1,label=pdf_name_2+" Dim "+str(i_dim))
		ax.set_xlabel("$x$",fontsize=fs1)
		ax.set_ylabel("$PDF(x)$",fontsize=fs1)
		plt.legend(loc="upper left")
		fig_file_name="pdf_comparison_d"+str(i_dim)+".pdf"
		plt.savefig(fig_file_name)

	return
################################################################################
def map2pce(pc_model,rvs_in,verbose=0):
	"""Obtain PC representation for the random variables that are described by samples.
	Employ a Rosenblatt transformation to build a map between the input RVs and the space
	of the PC germ.
	Input:
		pc_model: object with properties of the PCE to be constructed
		rvs_in	: numpy array with input RV samples. Each line is a sample. The columns
		          represent the dimensions of the input RV.
		verbose : verbosity level (more output for higher values)

	Output:
		Numpy array with PC coefficients for each RV in the original rvs_in input
	"""

	# Dimensionality and number of samples of input RVs
	ndim = rvs_in.shape[1]
	nsamp = rvs_in.shape[0]

	# Algorithm parameters
	bw = -1 # KDE bandwidth for Rosenblatt (on interval 0.1)
	iiout = 500 # interval for output to screen

	# Number of PCE terms
	npce = pc_model.GetNumberPCTerms()

	# Get the default quadrature points
	qdpts = uqtkarray.dblArray2D()
	pc_model.GetQuadPoints(qdpts)

	totquat = pc_model.GetNQuadPoints()
	print "Total number of quadrature points =",totquat

	# Set up transpose of input data for the inverse Rosenblatt transformation in a UQTk array
	ydata_t = uqtkarray.dblArray2D(ndim,nsamp)
	ydata_t.setnpdblArray(np.asfortranarray(rvs_in.T))

	# Set up numpy array for mapped quadrature points
	invRosData = np.zeros((totquat,ndim))

	# Map all quadrature points in chosen PC set to the distribution given by the data
	# using the inverse Rosenblatt transformation
	for ipt in range(totquat):

		# print "Converting quadrature point #",ipt

		# Set up working arrays
		quadunif = uqtkarray.dblArray1D(ndim,0.0)
		invRosData_1s = uqtkarray.dblArray1D(ndim,0.0)

		# First map each point to uniform[0,1]
		# PCtoPC maps to [-1,1], which then gets remapped to [0,1]
		for idim in range(ndim):
			quadunif[idim] = (uqtktools.PCtoPC(qdpts[ipt,idim],pc_model.GetPCType(),pc_model.GetAlpha(),pc_model.GetBeta(),"LU",0.0,0.0)+1.0)/2.0

		# Map each point from uniform[0,1] to the distribution given by the original samples via inverse Rosenblatt
		if bw > 0:
			uqtktools.invRos(quadunif,ydata_t,invRosData_1s,bw)
		else:
			uqtktools.invRos(quadunif,ydata_t,invRosData_1s)

		# Store results
		for idim in range(ndim):
			invRosData[ipt,idim] = invRosData_1s[idim]

		# Screen diagnostic output
		if ((ipt+1)%iiout == 0) or ipt==0 or (ipt+1)==totquat:
			print "Inverse Rosenblatt for Galerkin projection:",(ipt+1),"/",totquat,"=",(ipt+1)*100/totquat,"% completed"

	# Get PC coefficients by Galerkin projection
	# Set up numpy array for PC coefficients (one column for each transformed random variable)
	c_k = np.zeros((npce,ndim))

	# Project each random variable one by one
	for idim in range(ndim):

		# UQTk array for PC coefficients for one variable
		c_k_1d = uqtkarray.dblArray1D(npce,0.0)

		# UQTk array for map evaluations at quadrature points for that variable
		invRosData_1d = uqtkarray.dblArray1D(totquat,0.0)
		# invRosData_1d.setnpdblArray(np.asfortranarray(invRosData[:,idim])
		for ipt in range(totquat):
			invRosData_1d[ipt]=invRosData[ipt,idim]

		# Galerkin Projection
		pc_model.GalerkProjection(invRosData_1d,c_k_1d)

		# Put coefficients in full array
		for ip in range(npce):
			c_k[ip,idim] = c_k_1d[ip]

	# Return numpy array of PC coefficients
	return c_k
################################################################################
def evaluate_pce(pc_model,pc_coeffs,germ_samples):
	"""
	Evaluate PCE at a set of samples of the germ of this PCE
	Input:
		pc_model: PC object with into about PCE
		pc_coeffs: numpy array with PC coefficients of the RVs to be evaluated.
		           Each column corresponds to one RV.
		germ_samples: numpy array with samples of the PCE grem at which the RVs
		              are to be evaluated. Each line is one sample. The number
					  of colums is the number of RVs.

	Output:
		Numpy array with PCE evaluations
	"""

	# Get data set dimensions etc.
	n_test_samples = germ_samples.shape[0]
	ndim = germ_samples.shape[1]
	npce = pc_model.GetNumberPCTerms()

	# Put PC germ samples in a UQTk array
	std_samples_uqtk = uqtkarray.dblArray2D(n_test_samples, ndim)
	std_samples_uqtk.setnpdblArray(np.asfortranarray(germ_samples))

	# Numpy array to store all RVs evaluated from sampled PCEs
	rvs_sampled = np.zeros((n_test_samples,ndim))

	# Evaluate PCE for RVs in each dimension
	for idim in range(ndim):

		# Create and fill UQTk array for PC coefficients
		c_k_1d_uqtk = uqtkarray.dblArray1D(npce,0.0)
		for ip in range(npce):
			c_k_1d_uqtk[ip] = pc_coeffs[ip,idim]

		# Create UQTk array to store outputs in
		rv_from_pce_uqtk = uqtkarray.dblArray1D(n_test_samples,0.0)

		# Evaluate the PCEs for reach input RV at those random samples
		pc_model.EvalPCAtCustPoints(rv_from_pce_uqtk,std_samples_uqtk,c_k_1d_uqtk)

		# Put evaluated samples in full 2D numpy array
		for isamp in range(n_test_samples):
			rvs_sampled[isamp,idim] = rv_from_pce_uqtk[isamp]

	# return numpy array of PCE evaluations
	return rvs_sampled
