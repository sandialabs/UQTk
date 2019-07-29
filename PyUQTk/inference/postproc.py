#!/usr/bin/env python
#=====================================================================================
#                     The UQ Toolkit (UQTk) version @UQTKVERSION@
#                     Copyright (@UQTKYEAR@) Sandia Corporation
#                     http://www.sandia.gov/UQToolkit/
#
#     Copyright (@UQTKYEAR@) Sandia Corporation. Under the terms of Contract DE-AC04-94AL85000
#     with Sandia Corporation, the U.S. Government retains certain rights in this software.
#
#     This file is part of The UQ Toolkit (UQTk)
#
#     UQTk is free software: you can redistribute it and/or modify
#     it under the terms of the GNU Lesser General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
#
#     UQTk is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU Lesser General Public License for more details.
#
#     You should have received a copy of the GNU Lesser General Public License
#     along with UQTk.  If not, see <http://www.gnu.org/licenses/>.
#
#     Questions? Contact Bert Debusschere <bjdebus@sandia.gov>
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
from scipy import stats, mgrid, c_, reshape, random, rot90
from matplotlib.ticker import MultipleLocator, FormatStrFormatter



try:
    import pymc
    have_pymc = True
# except ImportError:
#     import pymc3 as pymc
#     have_pymc = True
except ImportError:
    print("PyMC is required for some of the MCMC postprocessing codes.")
    print ("Will proceed without, but some convergence tests will not be available.")
    have_pymc = False

##################################################################################
# Compute autocorrelation a one-dimensional set of samples
# Main function is acor(X), where X is a numpy array of samples
##################################################################################
from numpy import *
from matplotlib.pyplot import *
def acor_in(X, MAXLAG, WINMULT):

    # compute mean of X
    L = len(X)
    mu = mean(X)
    Xhat = X - mu
    std = sqrt(var(X))

    iMax = L - MAXLAG

    # compute autocorrelation
    # sind = arange(MAXLAG+1)
    iind = arange(iMax)
    C = zeros(MAXLAG + 1)
    for s in range(MAXLAG+1):
        C[s] += sum(Xhat[iind]*Xhat[iind+s])
    C *= 1./iMax

    D = C[0] # diffusion coeff
    D += 2*sum(C[1:])

    sigma = sqrt(abs(D/L))
    tau = D/C[0]

    # print D, L, sigma, tau, tau*WINMULT, MAXLAG
    return C[0], D, L, sigma, tau, tau*WINMULT, Xhat
# take in 1d numpy array of samples X
def acor(X,MAXLAG = 10, WINMULT = 5):
    C0, D, L, sigma, tau, tauWINMULT, X = acor_in(X, MAXLAG, WINMULT)
    # print tau, sigma
    Ls = []
    S = []
    while tau*WINMULT >= MAXLAG:
        Lh = L/2
        Ls.append(Lh)
        j1,j2 = 0,1
        for i in range(Lh):
            X[i] = X[j1] + X[j2]
            j1 += 2
            j2 += 2
        _, D, L, sigma, tau, tauWINMULT, X = acor_in(X[:Lh], MAXLAG, WINMULT)
        S.append(sigma)
    if len(S) == 0:
        S.append(sigma)
        Ls.append(L)

    sigma = S[-1]
    Ls = 2*array(Ls[::-1])
    for i in range(len(S)):
        D = .25 * sigma**2 * Ls[i]
        tau = D/C0
        sigma = sqrt(D/Ls[i])

    return tau
###################################################################################################
def compute_group_auto_corr(v,maxlag):
  """Compute autocorrelation of v, an array where each column is a set of samples,
  for a lag ranging from 0 to maxlag-1. Ouputs numpy array with autocorrelation."""

  # Get dimensions of input array with samples
  n_pts = np.shape(v)[0]
  n_var = np.shape(v)[1]

  # Initialize array
  auto_corr = np.zeros((maxlag,n_var))

  # Get mean and variance of v for each variable over the samples provided
  v_m = v.mean(0)
  v_var = v.var(0)


  # Subtract the mean of v
  v_nm = v - v_m

  # Compute autocovariance of v over all variables
  for lag in range(maxlag):
    n_sum = n_pts - lag     # total number of terms in sum
    for i in range(n_sum):
      auto_corr[lag,:] += v_nm[i,:]*v_nm[i+lag,:]
    auto_corr[lag,:] /= float(n_sum)

  # Normalize by variance
  auto_corr /= v_var

  return auto_corr
###################################################################################################
def compute_auto_corr(v,maxlag):
  """Compute autocorrelation of v (1D vector of samples) for a lag ranging from 0 to maxlag-1.
  Ouputs numpy array with autocorrelation."""

  # Initialize array
  auto_corr = np.zeros(maxlag)

  # Get mean and variance of v
  v_m = v.mean()
  v_var = v.var()
  n_pts = len(v)

  # Subtract the mean of v
  v_nm = v - v_m

  # Compute autocovariance of v
  for lag in range(maxlag):
    n_sum = n_pts - lag     # total number of terms in sum
    for i in range(n_sum):
      auto_corr[lag] += v_nm[i]*v_nm[i+lag]
    auto_corr[lag] /= float(n_sum)

  # Normalize by variance
  auto_corr /= v_var

  return auto_corr
###################################################################################################
def plot_auto_corr(v,vname):
  """Plot autocorrelation (in v), for variable with name vname"""

  # Set up the figure
  fig = plt.figure(figsize=(8,6))
  ax = fig.add_axes([0.10, 0.10, 0.85, 0.85])
  l1 = plt.plot(v,label=vname)

  ax.set_xlabel("lag")
  ax.set_ylabel("autocorrelation")
  plt.ylim([-0.1,1])

  plt.legend()
  plt.savefig("corr_"+vname+".pdf")
  plt.close()

  return
###################################################################################################
def compute_effective_sample_size(n_sam,auto_corr):
  """Computes the effective sample size for a vector of samples
  by dividing the number of samples (n_sam) by the integral
  of the autocorrelation (auto_corr) between the samples. (i.e. the more
  correlated successive samples are, the less independent samples there are
  in the chain.)
  The algorithm is based on:
      Markov Chain Monte Carlo in Practice: A Roundtable Discussion
      Robert E. Kass, Bradley P. Carlin, Andrew Gelman and Radford M. Neal
      The American Statistician, Vol. 52, No. 2 (May, 1998), pp. 93-100
      Published by: American Statistical Association
      Article DOI: 10.2307/2685466
      Article Stable URL: http://www.jstor.org/stable/2685466
  """

  # Length of autocorrelation array
  n_ac = len(auto_corr)

  # Find the lag where the autocorrelation goes to zero (or below)
  i_zero = 1   # start at lag 1 since the autocorrelation has value 1.0 at lag 0 by definition
  done = False
  while (i_zero < n_ac and not done):
    if auto_corr[i_zero] > 0.0:
      i_zero += 1
    else:
      done = True

  if i_zero == n_ac:
    print("WARNING: Autocorrelation did not go to zero within range provided")

  # Integral relies on symmetry and the fact that autocorrelation is 1 at zero lag
  ESS = int(n_sam // (1.0+2.0*np.sum(auto_corr[1:i_zero])))

  return ESS
###################################################################################################
###################################################################################################
def read_remaining_lines(samples_file,n_burnin,stride):
    """Read all remaining lines in MCMC sample filename, leaving out the first n_burnin,
    and only one in every stride lines after that
    """
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

    return samples_list

###################################################################################################
def remove_MAP_line(samples_list,debug):
    """Remove the last line if is has a value < 0 (i.e. -1) in the acceptance_prob column (next to last)"""
    if(samples_list[-1][-2] < 0):
        # Remove the last line
        del samples_list[-1]
        if (debug > 0):
            print("The last sample line has been deleted as it contained the MAP values")

###################################################################################################
def extract_vars(samples_file_name,n_burnin,v_names,debug,stride=1):
    """From a file with samples in ascii format, with
    the first line containing the label for each column, extract
    the columns with the labels in v_names and return them
    in a numpy array. Remove n_burnin samples from the top.
    Only read one in every stride number of lines after that.
    Assumes that the first column is the MCMC step number, the next to last column is the acceptance
    probability for each step, and the last column is the posterior probability for each step. The
    last line is removed if it contains -1 for the acceptance probability (which means this
    line contains the MAP values)"""

    # Open text file with all samples,
    samples_file = open(samples_file_name,"r")

    # Extract first line with the column labels and find the column
    # numbers corresponding to the variables of interest.
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

    # Read subsequent lines
    samples_list = read_remaining_lines(samples_file,n_burnin,stride)

    # Close the file
    samples_file.close()

    # Remove MAP values, if present
    remove_MAP_line(samples_list,debug)

    # Convert list to array
    steady_samples = np.array(samples_list)


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

###################################################################################################
def extract_all_vars(samples_file_name,n_burnin,debug,stride=1,labels=True):
    """Extract samples and labels from an MCMC chain file.
    Assumes the following:
        * The file is in ASCII format
        * The first column contains the MCMC step number
        * The next to last column contains the acceptance probability for the jump proposed in this step
        * The last column contains the posterior probability of the state in this step
        * The columns in between contain the sampled states
        * Unless the argument labels == False, the first line contains labels for each column
        * If the last line has -1 in the acceptance probability column, then this line contains
          the MAP values. This line is removed before returning the samples to the calling program.
    Arguments:
        * samples_file_name: name of file to parse
        * n_burnin: number of lines to skip from the top
        * debug: higher values are more verbose in output
        * stride: stride to take in parsing sample lines. [default = 1]
        * labels: True if the file contains column labels in first line. False if not. [default = True]
                  If not column labels are present, they are manufactured as aa, ab, ac, ..., az, ba, bb, ...
    Returns:
        * A numpy array with samples (one sample of all parameters per row)
        * A list of variable names
    """

    # Open text file with all samples
    samples_file = open(samples_file_name,"r")


    if (labels): # Column labels are present in first line
        # Extract first line with the column labels
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
    else:
        # Extract first line to see how many columns we have
        first_line = samples_file.readline().rstrip('\n')
        first_line_items = [item for item in first_line.split()]

        # Identify the MCMC vars, knowing that the first column is the step
        # number and the last two columns are acceptance and posterior prob
        n_cols = len(first_line_items)
        n_vars = n_cols - 3

        # Generate variable names as aa, ab, ..., az, ba, bb, ...
        if (n_vars > 52*26):    # only 52 entries in string.letters. If need be, could go higher by allowing aA, aB, ... , Aa, ...
            print("In routine extract_all_vars: too many columns for automatic name generation")
            sys.exit(1)

        v_names = []
        for i_v in range(n_vars):
            name = ""
            name += string.ascii_letters[i_v // 26]
            name += string.ascii_letters[i_v%26]
            v_names.append(name)

        if (debug > 0):
            print("There are",n_cols," columns in file", samples_file_name)
            print("MCMC chain variables have been labeled", v_names)

        # Rewind file so the first line will be read just like the other sample lines
        samples_file.seek(0)

    # Read subsequent lines
    samples_list = read_remaining_lines(samples_file,n_burnin,stride)

    # Close the file
    samples_file.close()

    # Test to make sure the number of samples is not 0 (which can happen in n_burnin > number of samples in file)
    if len(samples_list) < 1:
        print("\nERROR: No valid samples were extracted. Check the sample file, or the specified n_burnin")
        sys.exit(1)

    # Remove MAP values, if present
    remove_MAP_line(samples_list,debug)

    # Convert list to array
    samples = np.array(samples_list)

    n_samples = samples.shape[0]

    if (debug > 0):
        print("Read in", n_samples, "regular samples of", n_vars, "variables from file", samples_file_name)

    return samples, v_names

###################################################################################################
# def effective_sample_sizes(var_samples,par_mean,par_cov):
#     """Computes the effective sample size for each column
#     by dividing the number of samples by the integral
#     of the autocorrelation between the samples. (i.e. the more
#     correlated successive samples are, the less independent samples there are
#     in the chain.)
#     The algorithm is based on:
#         Markov Chain Monte Carlo in Practice: A Roundtable Discussion
#         Robert E. Kass, Bradley P. Carlin, Andrew Gelman and Radford M. Neal
#         The American Statistician, Vol. 52, No. 2 (May, 1998), pp. 93-100
#         Published by: American Statistical Association
#         Article DOI: 10.2307/2685466
#         Article Stable URL: http://www.jstor.org/stable/2685466
#     """
#
#     # Number of variable samples in set
#     n_sam = var_samples.shape[0]
#     # Number of variables in this sample set
#     n_vars = var_samples.shape[1]
#
#     # Array to store effective sample sizes in
#     ess = []
#
#     # Cut-off point for autocorrelation
#     # Ideally, n_a should be chosen such that the autocorrelation goes to 0 at this lag.
#     # Chosing n_a too low will give inaccurate results (overpredicting ESS), but going
#     # to much higher lag will create a lot of noise in ESS estimate.
#     n_a = min(100,n_sam)
#     for i_v in range(n_vars):
#         # Subtract mean from current variable samples
#         v_nm = var_samples[:,i_v] - par_mean[i_v]
#         # Compute autocorrelation for this variable. np.autocorrelate returns vector with
#         # lag from -n_sam to n_sam, with the 0 shift in the middle. Only retain from lag 0 to n_a.
#         r_v = np.correlate(v_nm, v_nm, mode = 'full')[-n_sam:-n_sam+n_a]
#         # Devide by number of samples in each sum, and normalize by variance
#         # (note: 0 lag has n_sam samples in sum, lag i has (n_sam - i) samples in sum
#         r_a = r_v / (par_cov[i_v,i_v]*(np.arange(n_sam, n_sam-n_a, -1)))
#         # Plot autocorrelation to see if n_a is large enough
#         #pl1,=plt.plot(r_a)
#         #plt.show()
#         # Effective Sample Size (Number of samples devided by integral of autocorrelation)
#         # Integral relies on symmetry and the fact that r_a is 1 at zero lag
#         ess.append(n_sam / (1.0+2.0*np.sum(r_a[1:])))
#
#     return ess

###################################################################################################
def plot_all_posteriors(d0,vnames,np_kde,out_file_base,debug=0,dense=False,truth=[]):
    """
    Given a set of samples of random variables, this script plots a lower triangular
    matrix of marginalized densities. The diagonal contains the density of individual
    random variables, marginalized over all other variables. Plots below the diagonal
    contain the 2D density of the associated pair of random variables, marginalized over
    all other variables.
    For chains with many variables, the "dense" option can be selected, which plots the
    triangular set of densities for the full chain with minimum spacing and labels, so that
    it is less cluttered. In this mode, this function also writes out a set of plots
    with the same posterior information, but just for two variables at the time,
    which is easier to read.

    Arguments:
        d0      : Set of samples, one column per variable (no extra columns)
        vnames  : Variable names
        np_kde  : Number of grid points at which to compute posterior densities with KDE at
        out_file_base: Base name for output files with plots
        debug   : >0 writes more output to screen (and even more if >1) [Defaults to 0]
        dense   : Set to True if dense output desired [Defaults to False]. The "dense" output
                  format puts all plots in the triangular format up against each other, without
                  any axis labels or space in between them. It is useful when plotting the
                  posteriors of a chain with many variables.
    """
    # Some settings to connect with code Cosmin gave me
    nthin = 1                # take only every nthin state (for faster kde)
    nskip = 0                # entries to skip
    istart = 0               # number of column with first MCMC variable
    cend   = 0               # extra columns at end to be removed

    nvars=d0.shape[1]-istart-cend  # number of variables we will actually process
    print('Number of sample lines in file',d0.shape[0])
    print('Number of vars  we will process in file',nvars)

    # Section 2
    # set up 2D kde objects
    print("Setting up 2D KDE objects")
    kern_i_j=[]
    for j in range(istart+1,istart+nvars):
        for i in range(istart,j):
            if (debug > 2):
                print(i,j)
            kern_i_j.append(stats.kde.gaussian_kde(c_[d0[nskip::nthin,i],d0[nskip::nthin,j]].T))

    # Section 3
    # set up 2D meshes and evaluate kde objects on those meshes
    # no. of grid points is controlled with kde_idx, defaults to 100
    print("Evaluating 2D KDE objects on meshes. This may take a while ...")
    kde_idx = np_kde*1j # complex number to include end points
    xmesh=[]; ymesh=[]; zmesh=[];
    icount=0
    cov_idx = np.zeros((nvars,nvars),dtype=np.int) # 2D array to keep track of which index in xmesh etc. the
                                                       # the plots corresponding to vars i,j belong to
    for j in range(istart+1,istart+nvars):
        for i in range(istart,j):
            if (debug > 0):
                print("Computing 2D marginal distribution between variables:",i,",",j,":",vnames[i]," & ",vnames[j])
            x,y = mgrid[d0[nskip:,i].min():d0[nskip:,i].max():kde_idx, d0[nskip:,j].min():d0[nskip:,j].max():kde_idx]
            z   = reshape(kern_i_j[icount](c_[x.ravel(), y.ravel()].T).T, x.T.shape)
            xmesh.append(x);
            ymesh.append(y);
            zmesh.append(z);
            cov_idx[i,j] = icount
            icount = icount+1

    # Section 4
    # evaluate 1D pdfs
    print("Evaluating 1D marginal pdfs with KDE")
    xlin=[]; pdflin=[];
    for i in range(istart,istart+nvars):
        xlin.append(np.linspace(d0[nskip:,i].min(),d0[nskip:,i].max(),np_kde)) ;
        kernlin=stats.kde.gaussian_kde(d0[nskip::nthin,i]);
        pdflin.append(kernlin(xlin[i-istart]));

    # Formatting for plots
    fs1=20 # Font size
    lw1=2 # Line width

    if (not dense):
        # Section 5
        print("Assembling lower-triangular plots in non-dense format")

        # ds is the distance between subplots
        # xs,ys are the coordinates (normalized) of the subplot in the lower left corner
        # xe,ye are the distances left in the uppper right corner
        # fsizex, fsizey are figure sizes
        # ncont are no of contours for 2D pdfs
        xs=0.12; ys=0.1; ds=0.04
        xe=0.08; ye=0.05
        fsizex=12; fsizey=12;
        ncont=20;
        sx=(1-(nvars-1)*ds-xs-xe)/nvars;
        sy=(1-(nvars-1)*ds-ys-ye)/nvars;
        majorFormatter = FormatStrFormatter('%6.0e')

        figname=out_file_base+".lowertriangle.pdf"   # figure name

        fig = plt.figure(figsize=(fsizex,fsizey))

        # Section 5.1
        subs=[]
        # add diagonal plots
        for i in range(nvars):
            subs.append(fig.add_axes([xs+i*(sx+ds),ys+(nvars-1-i)*(sy+ds),sx,sy]))

        # add lower triangular plots
        for i in range(nvars-1):
            for j in range(i+1):
                if (debug > 2):
                    print(j,(nvars-2-i))
                subs.append(fig.add_axes([xs+j*(sx+ds),ys+(nvars-2-i)*(sy+ds),sx,sy]))

        subsnp=np.array(subs)

        # Plot 1D pdfs
        for i in range(nvars):
            subsnp[i].plot(xlin[i],pdflin[i])

        # Plot 2D pdfs
        for i in range(nvars*(nvars-1) // 2):
            subsnp[nvars+i].contour(xmesh[i],ymesh[i],zmesh[i],ncont)

        # Section 5.2
        # just a few ticks and ticklabels
        for subpl in subsnp:
        #     subpl.set_xticks([])
        #     subpl.set_yticks([])
            subpl.locator_params(tight=True, nbins=5)

        # for diagonal plots, put no ticks and lables on y-axis
        # and no grid on the plots
        for i in range(istart,istart+nvars):
        #     subsnp[i-istart].set_xticks([d0[nskip:,i].min(),d0[nskip:,i].max()]);
            subsnp[i-istart].set_yticks([])
            subsnp[i-istart].grid(False)

        # Set y labels on the right for diagonal plots
        for i in range(nvars):
            subsnp[i].yaxis.tick_right()
            subsnp[i].yaxis.set_label_position("right")
            subsnp[i].set_ylabel(vnames[i], fontsize=fs1)
            #subsnp[i].set_ylabel(r'$'+vnames[i]+'$', fontsize=fs1)

        plt.savefig(figname)
        plt.close()

    else:
        # Section 5
        # Dense plot format: print full lower-triangular matrix but w/o any white space, tick marks or labels.
        print("Assembling lower-triangular plots in dense format")

        # ds is the distance between subplots
        # xs,ys are the coordinates (normalized) of the subplot in the lower left corner
        # xe,ye are the distances left in the uppper right corner
        # fsizex, fsizey are figure sizes
        # ncont are no of contours for 2D pdfs
        xs=0.12; ys=0.1; ds=0.0
        xe=0.08; ye=0.05
        fsizex=12; fsizey=12;
        ncont=10;
        sx=(1-(nvars-1)*ds-xs-xe)/nvars;
        sy=(1-(nvars-1)*ds-ys-ye)/nvars;
        majorFormatter = FormatStrFormatter('%6.0e')

        figname=out_file_base+".lowertriangle-dense.pdf"   # figure name

        fig_d = plt.figure(figsize=(fsizex,fsizey))

        # Section 5.1
        subs=[]
        # add diagonal plots
        for i in range(nvars):
            subs.append(fig_d.add_axes([xs+i*(sx+ds),ys+(nvars-1-i)*(sy+ds),sx,sy]))

        # add lower triangular plots
        for i in range(nvars-1):
            for j in range(i+1):
                if (debug > 2):
                    print(j,(nvars-2-i))
                subs.append(fig_d.add_axes([xs+j*(sx+ds),ys+(nvars-2-i)*(sy+ds),sx,sy]))

        subsnp=np.array(subs)

        # Plot 1D pdfs along diagonals
        for i in range(nvars):
            subsnp[i].plot(xlin[i],pdflin[i])
            if len(truth) > 0:
              subsnp[i].axvline(truth[i],ls='--',color='tab:gray',alpha=.5)

        # Plot 2D pdfs
        for i in range(nvars*(nvars-1) // 2):
            subsnp[nvars+i].contour(xmesh[i],ymesh[i],zmesh[i],ncont)

        # Section 5.2
        # no ticks and ticklabels on most plots
        for subpl in subsnp:
            subpl.set_xticks([]);
            subpl.set_yticks([]);

        # Set variable names as title for diagonal marginal plots
        for i in range(nvars):
            subsnp[i].set_title(vnames[i], fontsize=fs1)
            # subsnp[i].yaxis.set_label_position("right")
            # # Plot variable names along diagonal on the right. Plot them at angle
            # # to make them easier to read. Hack to add some white space in front of label
            # # so that it does not overlap with plot frame.
            # subsnp[i].set_ylabel("     "+vnames[i], fontsize=fs1, rotation=45)
        #     #subsnp[i].set_ylabel(r'$'+vnames[i]+'$', fontsize=fs1)


        plt.savefig(figname)
        plt.close()

        print("Assembling marginal density plots for all pairs of MCMC variables")

        # ds is the distance between subplots
        # xs,ys are the coordinates (normalized) of the subplot in the lower left corner
        # xe,ye are the distances left in the uppper right corner
        # fsizex, fsizey are figure sizes
        # ncont are no of contours for 2D pdfs
        xs=0.12; ys=0.1; ds=0.04
        xe=0.08; ye=0.05
        fsizex=12; fsizey=12;
        ncont=20;
        nvars_sm=2
        sx=(1-(nvars_sm-1)*ds-xs-xe)/nvars_sm;
        sy=(1-(nvars_sm-1)*ds-ys-ye)/nvars_sm;
        fs1=20
        majorFormatter = FormatStrFormatter('%6.0e')


        # loop over all pairs of MCMC variables.
        for j in range(istart+1,istart+nvars):
            for i in range(istart,j):

                print("Plotting densities for variables",vnames[i],"and",vnames[j])
                figname=out_file_base + "." + vnames[i] + "-" + vnames[j] + ".pdf"

                fig_sm = plt.figure(figsize=(fsizex,fsizey))

                subs=[]
                # add diagonal plots
                subs.append(fig_sm.add_axes([xs        ,ys+(sy+ds),sx,sy]))  # marginal for var i
                subs.append(fig_sm.add_axes([xs+(sx+ds),ys        ,sx,sy]))  # marginal for var j

                # add lower triangular plot
                subs.append(fig_sm.add_axes([xs        ,ys        ,sx,sy]))  # marginal for vars i,j

                subsnp=np.array(subs)

                # Plot 1D pdfs
                subsnp[0].plot(xlin[i],pdflin[i])
                subsnp[1].plot(xlin[j],pdflin[j])

                # Plot 2D pdfs
                i_2D = cov_idx[i,j]
                subsnp[2].contour(xmesh[i_2D],ymesh[i_2D],zmesh[i_2D],ncont)

                # set just a few ticks and ticklabels
                for subpl in subsnp:
                    subpl.locator_params(tight=True, nbins=5)

                # no ticks and ticklabels on y axes on diagonals (first two plots in subsnp array)
                # no grid on diagonal plots
                for subpl in subsnp[0:2]:
                    subpl.set_yticks([])
                    subpl.grid(False)

                # for diagonal plots only put xmin and xmax
                #subsnp[0].set_xticks([d0[nskip:,i].min(),d0[nskip:,i].max()]);
                #subsnp[1].set_xticks([d0[nskip:,j].min(),d0[nskip:,j].max()]);


                # Set y labels on the right for diagonal plots
                #subsnp[0].yaxis.tick_right()
                subsnp[0].yaxis.set_label_position("right")
                subsnp[0].set_ylabel(vnames[i], fontsize=fs1)

                #subsnp[1].yaxis.tick_right()
                subsnp[1].yaxis.set_label_position("right")
                subsnp[1].set_ylabel(vnames[j], fontsize=fs1)

                # Write out figure
                plt.savefig(figname)
                plt.close()

###################################################################################################
def get_mcmc_stats(all_samples,v_names,out_file_base,debug,stats_types=['all']):
    """
    Generate statistics of the passed in MCMC samples.
    Assumes that the first column of all_samples contains the step number, and the last two
    columns contain the acceptance probability and the posterior probability for each sampled state.

    Inputs:
        all_samples   : Array with all samples (one sample set per row). Has step number in first
                        column and acceptance probability and posterior in last two columns
        v_names       : Actual variable names
        out_file_base : Base for output file names
        debug         : Writes out more info if number is larger
        stats_types   : List with types of statistics to compute. Default is ['all']
                        Valid types in the list are  'corr' (autocorrelations)
                                                     'conv' (chain convergence tests)
                                                     'all'  (computes all of these.)
                        Note: basic statistics (mean, map, covariance) are always computed.

    Outputs:
        Various statistics written to the screen
        Correlation functions written to pdf files
        Returns array of map values
    """

    # Number of variables, columns, samples in the file
    n_vars = len(v_names)
    n_cols = all_samples.shape[1]
    n_sam  = all_samples.shape[0]

    # Extract all MCMC chain variables in separate array
    var_samples = all_samples[:,1:1+n_vars]
    if (debug > 0):
        print(var_samples.shape)

    #
    # Compute basic stats
    #
    print("\nComputing basic statistics:\n")

    # Compute mean parameter values
    par_mean = np.mean(var_samples,axis=0,dtype=np.float64)

    #print "\nParameter mean values:\n"
    #for i_v in range(n_vars):
    #    print "  ", v_names[i_v], ":", par_mean[i_v]

    # Compute the covariance
    par_cov = np.cov(var_samples,rowvar=0)

    print("\nParameter covariances:\n")
    print(par_cov)

    # write out covariance matrix to file
    cov_file_name = out_file_base + ".covariance.dat"
    np.savetxt(cov_file_name,par_cov)

    # print the square root of the diagonal entries of the covariance
    #print "\nParameter standard deviations (proposal width estimates):\n"
    #for i_v in range(n_vars):
    #    print "  ", v_names[i_v], ":", math.sqrt(par_cov[i_v,i_v])

    #
    # Compute the MAP values
    # (could also get this from the last line of the MCMC output file
    # but this line is not always there; and it is more fun
    # to do it with Python)
    #

    # Sample index with max posterior prop (last column in MCMC file):
    i_map = all_samples[:,-1].argmax()

    print()
    print('%27s' % "Parameter :", '%15s' % "Mean Value", '%15s' % "MAP values", '%15s' % "Std. Dev.")
    for i_v in range(n_vars):
        print('%25s' % v_names[i_v], ":", '%15.8e' % par_mean[i_v], '%15.8e' % var_samples[i_map,i_v],end=" ")
        print('%15.8e' % math.sqrt(par_cov[i_v,i_v]))

    # Write mean and MAP to file
    mean_file_name = out_file_base + ".mean.dat"
    np.savetxt(mean_file_name,par_mean)

    map_file_name = out_file_base + ".map.dat"
    np.savetxt(map_file_name,var_samples[i_map,:])

    # Compute mean and standard deviation of acceptance probability
    print("\nAcceptance Probability:\n")

    # In some cases, the next to last column contains the ratio of posterior
    # values rather than the acceptance probability. First convert this number
    # to acceptance probabilities: acc_prob = min(alpha,1)
    # (This does no harm if the next to last column already contains the actual acceptance probability)
    acc_prob = np.minimum(all_samples[:,-2],np.ones_like(all_samples[:,-2]))
    # In some cases, a very large negative number is shown in the column for acceptance
    # probability to indicate a proposed value was out of bounds. In that case, replace
    # the value with 0. Again, this does no harm if the next to last column already contains
    # the actual acceptance probability.
    acc_prob = np.maximum(acc_prob,np.zeros_like(acc_prob))

    print("Mean     :",'%8.4f' % acc_prob.mean())
    print("Std. Dev.:",'%8.4f' % acc_prob.std(),"\n")

    #
    # Compute correlations
    #
    if 'corr' in stats_types or 'all' in stats_types:
        # #
        # # Compute effective sample size  (ESS)
        # #
        # print "\nEffective Sample Sizes:\n"
        #
        # ess = effective_sample_sizes(var_samples,par_mean,par_cov)
        #
        # for i_v in range(n_vars):
        #     print "  ",v_names[i_v],":",int(ess[i_v]),"out of",n_sam

        #
        # Compute autocorrelations and effective sample size  (ESS)
        #
        print("\nAutocorrelations and Effective Sample Sizes:\n")

        # Number of variable samples in this file
        n_sam = var_samples.shape[0]

        # Cut-off point for autocorrelation
        # Ideally, n_a should be chosen such that the autocorrelation goes to 0 at this lag.
        # Chosing n_a too low will give inaccurate results (overpredicting ESS), but going
        # to much higher lag will create a lot of noise in ESS estimate.
        n_a = min(1000,n_sam)

        # Autocorrelation computation
        auto_corr_vars = compute_group_auto_corr(var_samples,n_a)

        # Plotting and computation of effective sample size
        for i_v in range(n_vars):
          # Plot autocorrelation to see if n_a is large enough
          plot_auto_corr(auto_corr_vars[:,i_v],v_names[i_v])
          # Effective Sample Size (Number of samples divided by integral of autocorrelation)
          ESS = compute_effective_sample_size(n_sam,auto_corr_vars[:,i_v])
          print("  ",v_names[i_v],":",ESS,"out of",n_sam," ; skip factor:",n_sam // ESS)

        print("\n  See plots corr-*.pdf for autocorrelations of chain samples for all variables.")

    #
    # Compute convergenve statistics (which rely on PyMC)
    #
    if ('conv' in stats_types or 'all' in stats_types) and have_pymc:

        #
        # Compute Raftery-Lewis convergence test
        #
        print("\nComputing Raftery-Lewis criteria for all variables\n")
        quant = 0.025 # Quantile level to be estimated
        relacc = 0.01  # Error in quantile level (relative to the mean of the parameter)
        conf = 0.95  # Confidence in the achieved accuray
        print ("  Computing # of samples needed to compute quantile",quant)
        print ("  for an accuracy",relacc*100,"% relative to parameter mean, with confidence",conf*100,"%:\n")
        print ("  Variable name: # initial samples to skip, # additional samples to take, thinning factor")
        for i_v in range(n_vars):
            output = pymc.raftery_lewis(var_samples[:,i_v], q=quant, r=relacc*par_mean[i_v], s=conf, verbose=0)
            print("  ",'%25s' % v_names[i_v], ":", '%8d' % output[2],",",'%8d' % output[3],",",'%8d' % output[4])

        print("\n")

        quant = 0.5 # Quantile level to be estimated
        print ("  Computing # of samples needed to compute quantile",quant)
        print ("  for an accuracy",relacc*100,"% relative to parameter mean, with confidence",conf*100,"%:\n")
        print ("  Variable name: # initial samples to skip, # additional samples to take, thinning factor")
        for i_v in range(n_vars):
            output = pymc.raftery_lewis(var_samples[:,i_v], q=quant, r=relacc*par_mean[i_v], s=conf, verbose=0)
            print ("  ",'%25s' % v_names[i_v], ":", '%8d' % output[2],",",'%8d' % output[3],",",'%8d' % output[4])

        #
        # Geweke test
        #
        print ("\nComputing Geweke test for all variables\n")
        print ("Geweke Test temporarily disabled. Needs to be debugged.")
        # for i_v in range(n_vars):
        #     var_scores = pymc.geweke(var_samples[:,i_v], intervals=20)
        #     pymc.Matplot.geweke_plot(var_scores, v_names[i_v])
        # print "  See plots *-diagnostic.png"

        #
        # Autocorrelations (done above already)
        #
        # print "\nComputing autocorrelations for all variables\n"
        # for i_v in range(n_vars):
        #     pymc.Matplot.autocorrelation(var_samples[:,i_v], v_names[i_v])
        # print "  See plots *-acf.png"

    return var_samples[i_map,:]

###################################################################################################

help_string = """
Usage:
  postproc.py [-h] -i <mcmc filename> [--nb <burn-in samples>]  [-s <stride>] \
              [--post_kde_pts n_kde_pts] [--nolabels] \
              [--stats [<stats_types>]]
What:
  Compute elementary statistics of MCMC chain
Where
  -h = print help info
  -i = name of file containing MCMC data
  -s = stride with which to read the file [defaults to 1]
  --post_kde_pts = if present: plot posteriors with KDE using the number of points specified in each dimension [defaults to not]
  --nb = number of burn-in samples to be removed from the chain [defaults to 0]
  --nolabels Indicates that the MCMC data file does not contain column labels (in which case they are generated)
  --stats = string that indicates type of stats requested. [Defaults to 'all']
            Valid types are  'corr' (autocorrelations)
                             'conv' (chain convergence tests)
                             'all'  (computes all of these.)
                             'none' (no advanced stats)
            Notes: * Basic statistics (mean, map, covariance) are always computed.
                   * This option can be specified more than once to add more stats to be computed.

Assumes the following:
    * The file is in ASCII format
    * The first column contains the MCMC step number
    * The next to last column contains the acceptance probability for the jump proposed in this step
    * The last column contains the posterior probability of the state in this step
    * The columns in between contain the sampled states
    * Unless the argument labels == False, the first line contains labels for each column
"""

if __name__ == "__main__":
    #
    # Process inputs
    #
    try:
        opts,v_names = getopt.getopt(sys.argv[1:],"hi:s:",["nb=","post_kde_pts=","stats=","nolabels"])
    except getopt.GetoptError as err:
        print(str(err))
        print(help_string)
        sys.exit(1)

    # Default values
    samples_file_name=""
    labels_present = True
    plot_posteriors = False
    n_kde_pts = 0
    stats_requested = ["all"]
    stats_user_specified = False
    stride = 1

    for o,a in opts:
        if o == "-h":
            print(help_string)
            sys.exit(0)
        elif o == "-i":
            samples_file_name = a
        elif o == "-s":
            stride = int(a)
        elif o == "--nb":
            n_burnin = int(a)
        elif o == "--nolabels":
            labels_present = False
        elif o == "--post_kde_pts":
            n_kde_pts = int(a)
        elif o == "--stats":
            if not stats_user_specified: # This is the first occurrence of the "--stats" option
                                         # So we reset the default list so we can add all the user-requested options
                stats_user_specified = True
                stats_requested = []
            stats_requested.append(a)
        else:
            assert False, "Unhandled command line parsing option. Use -h flag to get usage info."

    # error checking
    if(samples_file_name==""):
      print("Sample file name must be specified")
      print(help_string)
      sys.exit(1)

    if (n_burnin < 0):
      print("The number of burn-in samples needs to be >= 0")
      print(help_string)
      sys.exit(1)

    if (stride < 0):
      print("The file read stride needs to be >= 0")
      print(help_string)
      sys.exit(1)

    if (n_kde_pts > 0):
        plot_posteriors = True
    else:
        print("Posterior graphs will not be computed. If posterior graphs are desired, specify the --post_kde_pts option.")

    # Base name of file for outputting results
    out_file_base = samples_file_name + ".nb" + str(n_burnin) + ".s" + str(stride)

    # Set to 1 to get more output to screen
    # Set to > 1 to get a lot of output to screen
    debug = 1

    #
    # Import variables of interest from the MCMC data file
    #
    all_samples, v_names = extract_all_vars(samples_file_name,n_burnin,debug,stride,labels=labels_present)
    n_vars = len(v_names)
    var_samples = all_samples[:,1:1+n_vars]

    # Get statistics
    get_mcmc_stats(all_samples,v_names,out_file_base,debug,stats_requested)


    # Plot posteriors if desired
    if plot_posteriors:
        print("\nPlotting Posteriors:\n")
        # Skip first and last two columns in samples. The columns in between contain variable samples
        plot_all_posteriors(var_samples,v_names,n_kde_pts,out_file_base,debug)
