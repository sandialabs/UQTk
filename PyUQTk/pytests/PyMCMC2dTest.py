#=====================================================================================
#                     The UQ Toolkit (UQTk) version 3.0.4
#                     Copyright (2017) Sandia Corporation
#                     http://www.sandia.gov/UQToolkit/
#
#     Copyright (2017) Sandia Corporation. Under the terms of Contract DE-AC04-94AL85000
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
# include path to include PyUQTk
import sys
sys.path.append('../../')

try:
	from numpy import *
	from matplotlib.pyplot import *
except ImportError:
	"Need numpy and matplotlib to test PyUQTk"

try:
	import PyUQTk.array as uqtkarray
	import PyUQTk.mcmc as uqtkmcmc
except ImportError:
	print "PyUQTk array module not found"
	print "If installing in a directory other than the build directory, make sure PYTHONPATH includes the install directory"

import time

# temp = random.randn(1000.)
# a = uqtkarray.dblArray1D(1000,101.0)

class pyLikelihood(uqtkmcmc.LikelihoodBase):
	def eval(self,x):
		# a.getnpdblarray(temp)
		# temp = array(a.pint4py())
		x0 = x[0]
		x1 = x[1]
		return -(1-x0)*(1-x0) - 100*(x1 - x0*x0)*(x1 - x0*x0)

start = time.time()
# testing MCMC library
print '\n*****************\nTesting MCMC\n*****************\n'
print 'Setting LogPosterior function, L'
print 'L is defined in uqtk.cpp (Rosenbrock function)'
L = pyLikelihood()
print 'Testing logposterior function'
xstart = uqtkarray.dblArray1D(2,0)
print xstart
print 'L.eval(x) =  ', L.eval(xstart)

print 'Setting up the sampler'
mchain = uqtkmcmc.MCMC(L)
print 'Setting chain dim, type (ss), initial proposal covariance'
dim = 2
mchain.setChainDim(dim)
mchain.initMethod("ss")
g = uqtkarray.dblArray1D(dim,.1)
mchain.initChainPropCovDiag(g)
print 'Chain Setup:'
mchain.printChainSetup();
print 'Running chain to chain.dat ...'
nCalls = 100000
mchain.setOutputInfo("txt","chain.dat",nCalls,nCalls);
mchain.runChain(nCalls,xstart);
print 'loading samples and plotting'
thin = 25
samples = loadtxt('chain.dat')[3001:-1:thin,1:3]
figure()
plot(samples[:,0],samples[:,1],'.')

# get the likelihood information
print 'Getting samples into numpy array...'
# logprobs = zeros(nCalls);
mchain.getSamples()
samples = mchain.samples
npsamples = zeros((dim,nCalls))
samples.getnpdblArray(npsamples)
figure()
plot(npsamples[0][::thin],npsamples[1][::thin],'.g')

end = time.time()
print end - start

# show()
