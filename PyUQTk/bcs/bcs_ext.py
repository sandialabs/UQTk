%pythoncode %{

import numpy as np
import matplotlib.pyplot as mpl
import uqtkarray as uqtkarray
import pce as uqtkpce
import tools as uqtktools
from uqtkarray import uqtk2numpy, numpy2uqtk
# BCS already added to path in compilation and install


class bcsreg:
	'''
	Class to compute the bcs regression coefficients for a scalar function of ndim dimensions. 
	'''
	def __init__(self,ndim,pcorder,pctype):
		'''
		Construction has the following inputs:
		ndim	: (int) number of input dimensions (features)
		pcorder	: (int) the initial order of the polynomial (changes in the algorithm)
		pctype	: ('LU','HG') type of polynomial basis functions, e.g., Legendre, Hermite
		
		'''
		self.ndim = ndim # int
		self.pcorder = pcorder # int
		self.pctype = pctype # 'LU', 'HG'

		# generate multi index
		self.__mindex_uqtk = uqtkarray.intArray2D()
		uqtktools.computeMultiIndex(self.ndim,self.pcorder,self.__mindex_uqtk);
		self.mindex = uqtk2numpy(self.__mindex_uqtk)
		self.__mindex0_uqtk = self.__mindex_uqtk # keep original 

		# get projection/ Vandermonde matrix
		self.__Phi_uqtk = uqtkarray.dblArray2D()

		# check if compiled
		self.__compiled = False
		self.compile()

	def compile(self,l_init=0.0,adaptive=0,optimal=1,scale=.1,verbose=0):
		'''
		Setting up variables for the BCS algorithm. Most of the variables do not need to be set. Default settings are sufficient for more cases. See the C++ code for more information about variables. 
		'''
		# now we begin BCS routine
		# set work variables
		self.__newmindex_uqtk = uqtkarray.intArray2D() # for uporder iteration
		self.__sigma2_p = uqtktools.new_doublep() # initial noise variance
		self.__lambda_init = uqtkarray.dblArray1D() # hierarchical prior parameter
		self.__adaptive, self.__optimal, self.__scale, self.__verbose = adaptive,optimal,scale,verbose
		self.__weights_uqtk = uqtkarray.dblArray1D() # weights/ coefficients for basis
		self.__used_uqtk = uqtkarray.intArray1D() # index of weights retained (nonzero)
		self.__errbars_uqtk = uqtkarray.dblArray1D() # error bars for each weight
		self.__nextbasis_uqtk = uqtkarray.dblArray1D() # if adaptive
		self.__alpha_uqtk = uqtkarray.dblArray1D() # prior hyperparameter (1/gamma)
		self.__lambda_p = uqtktools.new_doublep() 
		
		uqtktools.doublep_assign(self.__lambda_p,l_init)

		self.__compiled = True

	def leastsq(self,X,y):
		'''
		perform simple least squares based on the original 
		pc order. 
		'''
		# convert input to uqtk arrays
		self.__X_uqtk = numpy2uqtk(X)
		self.__y_uqtk = numpy2uqtk(y)

		# get vandermonde matrix w.r.t. original pc basis
		self.__V_uqtk = uqtkarray.dblArray2D()
		self.__pcmodel0 = uqtkpce.PCSet("NISPnoq",self.__mindex0_uqtk,self.pctype,0.0,1.0) # initiate
		self.__pcmodel0.EvalBasisAtCustPts(self.__X_uqtk,self.__V_uqtk)
		self.Vandermonde = uqtk2numpy(self.__V_uqtk)
		self.__sol = np.linalg.lstsq(self.Vandermonde,y)
		return self.__sol[0], self.__sol[1]

	def fit(self,X,y,tol=1e-8,sigsq=None,upit=0):
		'''
		Train bcs model coefficients with X and y data
		X 		:	2d numpy of inputs/ feature data 
		y 		: 	1d numpy array of labels/ outputs
		tol 	:	tolerance (smaller means we keep more coefficients)
		sigsq	: 	initial noise set automatically based on y data
		upit 	: 	(int) number of iterations to add higher order terms 	  

		returns the polynomial coefficient (weights), and the mulitindex. One can also return the sensitivity indices by calling self.sens
		'''
		if self.__compiled == False:
			print "Need to compile first!"

		# convert numpy test data into uqtk data types
		self.__X_uqtk = numpy2uqtk(X)
		self.__y_uqtk = numpy2uqtk(y)
		self.Xtrain = X
		self.ytrain = y

		if sigsq == None: 
			self.__sigma2 = np.var(y)/1e2
		else: self.__sigma2 = sigsq
		uqtktools.doublep_assign(self.__sigma2_p,self.__sigma2)

		self.__tol = tol
		self.__upit = upit

		# begin uporder iterations
		for iter in range(self.__upit+1):

			# get projection/ Vandermonde matrix
			self.__pcmodel = uqtkpce.PCSet("NISPnoq",self.__mindex_uqtk,self.pctype,0.0,1.0) # initiate with new mindex
			self.__pcmodel.EvalBasisAtCustPts(self.__X_uqtk,self.__Phi_uqtk)
			self.__Phi = uqtk2numpy(self.__Phi_uqtk)

			# resest sigma parameter (if not, may get seg fault)
			uqtktools.doublep_assign(self.__sigma2_p,self.__sigma2)

			# change to uqtkbcs.BCS if testing outside source
			BCS(self.__Phi_uqtk,self.__y_uqtk,self.__sigma2_p,self.__tol,self.__lambda_init,self.__adaptive,self.__optimal,self.__scale,self.__verbose,self.__weights_uqtk,self.__used_uqtk,self.__errbars_uqtk,self.__nextbasis_uqtk,self.__alpha_uqtk,self.__lambda_p)

			# add new mulitindex to newmindex
			uqtkarray.subMatrix_row_int(self.__mindex_uqtk,self.__used_uqtk,self.__newmindex_uqtk)

			if iter < self.__upit :
				# redefine mindex = newmindex if still iterating
				self.__newmindex_added_uqtk = uqtkarray.intArray2D()
				uqtktools.upOrder(self.__newmindex_uqtk,self.__newmindex_added_uqtk)
				self.__mindex_uqtk = self.__newmindex_added_uqtk
				print "New mindex basis: ", uqtk2numpy(self.__mindex_uqtk)[len(self.__newmindex_uqtk):]

		# return new multiindex to create new pce model
		self.__pcmodel_new = uqtkpce.PCSet("NISPnoq",self.__newmindex_uqtk,self.pctype,0.0,1.0)
		self.mindex = uqtk2numpy(self.__newmindex_uqtk)
		eff_dim = self.ndim - sum(sum(self.mindex,0) == 0)

		self.weights = uqtk2numpy(self.__weights_uqtk)
		self.weight_index = uqtk2numpy(self.__used_uqtk)
		self.error_bars = uqtk2numpy(self.__errbars_uqtk)

		# get main effect sensitivity indices
		self.__main_eff_uqtk = uqtkarray.dblArray1D()
		self.__tot_eff_uqtk = uqtkarray.dblArray1D()
		self.__joint_eff_uqtk = uqtkarray.dblArray2D()
		self.__pcmodel_new.ComputeMainSens(self.__weights_uqtk,self.__main_eff_uqtk)
		self.__pcmodel_new.ComputeTotSens(self.__weights_uqtk,self.__tot_eff_uqtk)
		self.__pcmodel_new.ComputeJointSens(self.__weights_uqtk,self.__joint_eff_uqtk)
		self.main_eff = uqtk2numpy(self.__main_eff_uqtk)
		self.tot_eff = uqtk2numpy(self.__tot_eff_uqtk)
		self.joint_eff = uqtk2numpy(self.__joint_eff_uqtk)
		self.senssum = {"main effect": self.main_eff, "total effect": self.tot_eff, "joint effect": self.joint_eff}
		return self.weights, self.mindex

	def predict(self,Xtest):
		'''
		Predict values after training the data
		Xtest 	: 	2d numpy array 

		returns 1d numpy scalar array of predictions
		'''
		self.__Xtest_uqtk = numpy2uqtk(Xtest)
		self.__ytest_uqtk = uqtkarray.dblArray1D()
		self.__pcmodel_new.EvalPCAtCustPoints(self.__ytest_uqtk,self.__Xtest_uqtk,self.__weights_uqtk)
		self.__ytest = uqtk2numpy(self.__ytest_uqtk)
		return self.__ytest
	def getsens(self):
		'''
		return sensitivities as dictionary
		'''
		return self.senssum

%}
