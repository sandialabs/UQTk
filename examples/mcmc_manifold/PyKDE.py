import numpy as np
from scipy import stats
import sys
import os
import math

KDEfile = sys.argv[1]
samplefile = sys.argv[2]

prior = np.genfromtxt(KDEfile, unpack=True)
samples = np.genfromtxt(samplefile)

kernel = stats.gaussian_kde(prior)

# Remove existing file
if (os.path.isfile('tmcmc_lp.dat')):
    os.remove('tmcmc_lp.dat')

# Append to empty binary file
lpfile = open('tmcmc_lp.dat', 'ba')

# Don't want to overload memory, kernel call in sections
sections = samples.shape[0] / 5000

for i in range(math.ceil(sections)):
    startIndex = 5000*i
    endIndex = 5000*(i + 1)

    if (endIndex > samples.shape[0]):
        logsamples = samples[startIndex::, :]
    else:
        logsamples = samples[startIndex:endIndex, :]

    logprior = kernel.logpdf(logsamples.T)

    np.savetxt(lpfile, logprior, fmt="%.18e", delimiter = " ", newline = "\n")

lpfile.close()
