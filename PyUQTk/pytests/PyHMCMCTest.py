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
#     Questions? Contact the UQTk Developers at https://github.com/sandialabs/UQTk/discussions
#     Sandia National Laboratories, Livermore, CA, USA
#=====================================================================================
from __future__ import print_function # To make print() in Python 2 behave like in Python 3

try:
	from numpy import *
except ImportError:
	print("Need numpy to test PyUQTk")
try:
	from matplotlib.pyplot import *
except ImportError:
	print("Need matplotlib to test PyUQTk")
try:
	from acor import *
except ImportError:
	print("Need acor to test PyUQTk")

try:
	import PyUQTk.array as uqtkarray
except ImportError:
	print("PyUQTk array module not found")
	print("If installing in a directory other than the build directory, make sure PYTHONPATH includes the install directory")
try:
	import PyUQTk.mcmc as uqtkmcmc
except ImportError:
	print("PyUQTk mcmc module not found")
	print("If installing in a directory other than the build directory, make sure PYTHONPATH includes the install directory")
try:
	from PyUQTk.inference.mcmc import *
except ImportError:
	print("PyUQTk inference.mcmc module not found")
	print("If installing in a directory other than the build directory, make sure PYTHONPATH includes the install directory")
try:
	from PyUQTk.inference.postproc import *
except ImportError:
	print("PyUQTk inference.postproc module not found")
	print("If installing in a directory other than the build directory, make sure PYTHONPATH includes the install directory")

try:
	import time
except ImportError:
	print("Time module needed to test PyUQTk")

'''
Use HMCMC to get samples from banana shaped function
'''

def U(q,a=1.0,b=100.0):
	'''
	U(q) = -log(prior(q)*Likelihood(q|data))
	'''
	q = copy(atleast_2d(q))
	return b*(q[:,1] - q[:,0]**2)**2 + (a - q[:,0])**2

def grad_U(q):
	'''
	grad_U(q) = gradient vector of U at q
	'''
	q = copy(atleast_2d(q))
	dUdx = (-400*q[:,0]*(q[:,1] - q[:,0]**2) - 2*(1 - q[:,0]))[0]
	dUdy = (200*(q[:,1] - q[:,0]**2))[0]
	return array([dUdx,dUdy])

# def HMCMC(U,grad_U,eps,L,q):
# 	current_q = copy(q) # save current

# 	# generate current p
# 	# propcov = 4*array([[ 0.01175383,  0.02065261],[ 0.02065261,  0.04296117]])
# 	p = random.randn(len(current_q))
# 	# p = random.multivariate_normal([0,0],propcov)
# 	current_p = copy(p) # save current p

# 	# make half step for momentum used for leap frog step
# 	p = p - eps * grad_U(q)/2.0

# 	for i in range(L):
# 		# p = p - eps * grad_U(q)/2.0
# 		q = q + eps*p
# 		# p = p - eps * grad_U(q)/2.0
# 		if (i != L-1): p = p - eps*grad_U(q)

# 	# make a half step for momentum at the end
# 	p = p - eps * grad_U(q)/2.0

# 	# negate the momentum to make a symmetric proposal
# 	p = -p

# 	# Evaluate potential and kinetic energy
# 	current_U = U(current_q)[0]
# 	current_K = sum(current_p**2)/2.0
# 	proposed_U = U(q)[0]
# 	proposed_K = sum(p**2)/2.0

# 	# Accept or reject the state at end of trajectory, returning either
# 	# the position at the end of the trajectory or the initial position

# 	if (log(random.rand()) < current_U-proposed_U+current_K-proposed_K):
# 		return q
# 	else:
# 		alpha = 0
# 		return current_q

fig = figure()
ax1 = fig.add_subplot(2,2,1)
ax2 = fig.add_subplot(2,2,2)
ax3 = fig.add_subplot(2,2,3)
ax4 = fig.add_subplot(2,2,4)

# Test U(q)
N = 80
qx = linspace(-1.5,2.5,N)
qy = linspace(-.5,5,N)
qx,qy = meshgrid(qx,qy)
qz =  exp(-U(array(zip(qx.flatten(),qy.flatten()))))
qz.shape = (N,N)

# Test grad_U
dUdx = zeros((N,N))
dUdy = zeros((N,N))
for i in range(N):
	for j in range(N):
		dU = grad_U([qx[i,j],qy[i,j]])
		dUdx[i,j] = dU[0]
		dUdy[i,j] = dU[1]

# Test HMCMC
print('\n*****************\nTesting HMCMC\n*****************\n')
samples1 = []
qstart = array([1.0,1.0])
q = copy(qstart)
samples1.append(copy(q))
eps = .01
L = 150
M = 500
nburn = 300
thin1 = 1
for i in range(M):
	q = HMCMC(U,grad_U,eps,L,q)
	if i > 0:
		samples1.append(copy(q))
samples1 = array(samples1)

# ax1.plot(qstart[0],qstart[0],'ok')
# ax1.quiver(qx,qy,-dUdx,-dUdy,alpha=.1)
# ax1.contour(qx,qy,qz,20,alpha=.4)
# ax1.plot(samples1[nburn::thin1,0],samples1[nburn::thin1,1],'*k',alpha=.1)

'''
Use AMCMC for banana shaped function
'''
class pyLikelihood(uqtkmcmc.LikelihoodBase):
	def eval(self,x):
		x0 = x[0]
		x1 = x[1]
		return -(1-x0)*(1-x0) - 100*(x1 - x0*x0)*(x1 - x0*x0)

# testing MCMC library
print('\n*****************\nTesting AMCMC\n*****************\n')
Like = pyLikelihood()
xstart = uqtkarray.dblArray1D(2,1.0)
mchain = uqtkmcmc.MCMC(Like)
dim = 2
mchain.setChainDim(dim)
mchain.initMethod("am")
g = uqtkarray.dblArray1D(dim,.1)
mchain.initChainPropCovDiag(g)
nCalls = L*M
thin2 = thin1*L
mchain.setWriteFlag(0)
mchain.setOutputInfo("txt","chain.dat",M,nCalls);
mchain.runChain(nCalls,xstart);
mchain.getSamples()
samples = mchain.samples
samples2 = zeros((dim,nCalls))
samples.getnpdblArray(samples2)
samples2 = samples2.T
propcov = uqtkarray.dblArray2D(2,2,0)
mchain.getChainPropCov(propcov)
m1 = propcov[0,0]; m2 = propcov[1,1]

# ax2.contour(qx,qy,qz,250,alpha=.4)
# ax2.plot(samples2[L*nburn::thin2,0],samples2[L*nburn::thin2,1],'*g',alpha=.1)

# plot mixing of samples
ax3.plot(samples1[nburn::thin1,0],'k',alpha=.4)
ax3.plot(samples2[L*nburn::thin2,0],'g',alpha=.4)
ax4.plot(samples1[nburn::thin1,1],'k',alpha=.4)
ax4.plot(samples2[L*nburn::thin2,1],'g',alpha=.4)

print('acor using HMCMC', acor(samples1[nburn::thin1,0]))
print('acor using MCMC', acor(samples2[L*nburn::thin2,0]))
