
import numpy as np
import scipy.linalg as la

import pymuqModeling as mm
import pymuqSamplingAlgorithms as ms


class MyGauss(mm.PyGaussianBase):

    def __init__(self, mu, gamma):
        """
        Constructs the Gaussian from a mean vector mu and covariance matrix gamma.
        """
        mm.PyGaussianBase.__init__(self,mu)

        self.gamma = gamma
        self.L = la.cho_factor(gamma,lower=True)[0]

    def SampleImpl(self, inputs):
        """ Overloading this function is optional. """
        z = np.random.randn(self.gamma.shape[0])
        return self.GetMean() + self.ApplyCovSqrt(z)

    def ApplyCovariance(self, x):
        return self.gamma @ x

    def ApplyPrecision(self, x):
        return la.cho_solve((self.L,True),x)

    def ApplyCovSqrt(self,x):
        return self.L @ x

    def ApplyPrecSqrt(self,x):
        return la.solve_triangular(self.L.T , x, lower=False)





###############################
# Construct the target density

tgtMean = np.ones((2,))

rho = 0.8
std1 = 1.5
std2 = 0.5
tgtCov = np.array([[std1*std1, rho*std1*std2],
                   [rho*std1*std2, std2*std2]])

tgtGauss = mm.Gaussian(tgtMean,tgtCov)
tgtDens = mm.Gaussian(tgtMean, tgtCov).AsDensity()

###############################
# Construct a custom Gaussian distribution to use in Crank-Nicholson proposal
gauss = MyGauss(tgtMean, tgtCov)

##############################
# Build the MCMC Sampler
opts = dict()
opts['NumSamples'] = 2000 # Number of MCMC steps to take
opts['BurnIn'] = 10 # Number of steps to throw away as burn in
opts['PrintLevel'] = 3 # in {0,1,2,3} Verbosity of the output
opts['Beta'] = 0.75 # Crank Nicholson parameter

# Construct the sampling problem from the target density
problem = ms.SamplingProblem(tgtDens)

# Construct the CrankNicholson proposal
pcnProp = ms.CrankNicolsonProposal(opts, problem, gauss)

# Use the proposal to construct a Metropolis-Hastings kernel
kern = ms.MHKernel(opts,problem,pcnProp)

# Construct the MCMC sampler using this transition kernel
sampler = ms.SingleChainMCMC(opts, [kern])

#################################
# Run the MCMC sampler
x0 = [tgtMean]
samps = sampler.Run(x0)

#################################
# Look at the results
print('\nEffective Sample Size =\n', samps.ESS(), '\n')
print('Sample mean=\n', samps.Mean(), '\n')
print('Sample Covariance=\n', samps.Covariance(), '\n')
