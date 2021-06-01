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
#     Questions? Contact the UQTk Developers at <uqtk-developers@software.sandia.gov>
#     Sandia National Laboratories, Livermore, CA, USA
#=====================================================================================
import sys
sys.path.append('../pyuqtkarray')
sys.path.append('../pce')
sys.path.append('../pyuqtkarray_tools')
sys.path.append('../tools')
sys.path.append('../bcs')


import numpy as np
import matplotlib.pyplot as mpl
import pyuqtkarray as uqtkarray
import pyuqtkarray_tools as uqtkarray_tools
import _pce as uqtkpce
import _tools as uqtktools
import _bcs as bcs

# BCS already added to path in compilation and install

# cross validation splitting
def kfold_split(nsamples,nfolds,seed=13):
    '''
    return dictionary of training and testing pairs using k-fold cross-validation
    '''
    # returns split data where each data is one fold left out
    KK=nfolds
    rn = np.random.RandomState(seed)

    indp=rn.permutation(nsamples)
    split_index=np.array_split(indp,KK)


    cvindices = {}
    # create testing and training folds
    for j in range(KK):
        fold = j
        newindex = [split_index[i] for i in range(len(split_index)) if i != (fold)]
        train_ind = np.array([],dtype='int64')
        for i in range(len(newindex)): train_ind = np.concatenate((train_ind,newindex[i]))
        test_ind = split_index[fold]
        cvindices[j] = {'train index': train_ind, 'val index': test_ind}

    return cvindices

def kfoldCV(x,y,nfolds=3,seed=13):
    '''
     Splits data into training/testing pairs for kfold cross-val
    x is a data matrix of size n x d1, d1 is dim of input
    y is a data matrix of size n x d2, d2 is dim of output

    Test:
    xtest = array([i*ones(6) for i in range(1,501)])
    xtest = np.array([i*np.ones(6) for i in range(1,501)])
    ytest = xtest[:,:2]
    K,ci = kfoldCV(xtest,ytest)

    for k in K.keys():
        a = K[k]['xtrain'][:,0]
        b = K[k]['xval'][:,0]
        test = sum(np.in1d(a,b)) + sum(np.in1d(b,a))
        print test
     '''
    n,d1 = x.shape
    ynew = np.atleast_2d(y)
    if len(ynew) == 1: ynew = ynew.T # change to shape (n,1)
    _,d2 = ynew.shape
    cv_idx = kfold_split(n,nfolds,seed)

    kfold_data = {}
    for k in cv_idx.keys():
        kfold_data[k] = {
        'xtrain': x[cv_idx[k]['train index']],
        'xval': x[cv_idx[k]['val index']],
        'ytrain': np.squeeze(ynew[cv_idx[k]['train index']]),
        'yval': np.squeeze(ynew[cv_idx[k]['val index']])
        } # use squeeze to return 1d array

        # set train and test to the same if 1 fold
        if nfolds == 1:
            kfold_data[k]['xtrain'] = kfold_data[k]['xval']
            kfold_data[k]['ytrain'] = kfold_data[k]['yval']

    return kfold_data

class bcsreg:
    '''
    Class to compute the bcs regression coefficients for a scalar function of ndim dimensions.
    '''
    def __init__(self,ndim,pcorder,pctype):
        '''
        Construction has the following inputs:
        ndim    : (int) number of input dimensions (features)
        pcorder    : (int) the initial order of the polynomial (changes in the algorithm)
        pctype    : ('LU','HG') type of polynomial basis functions, e.g., Legendre, Hermite

        '''
        self.ndim = ndim # int
        self.pcorder = pcorder # int
        self.pctype = pctype # 'LU', 'HG'

        # generate multi index
        self.__mindex_uqtk = uqtkarray.intArray2D()
        uqtktools.computeMultiIndex(self.ndim,self.pcorder,self.__mindex_uqtk);
        self.mindex = uqtkarray_tools.uqtk2numpy(self.__mindex_uqtk)
        self.__mindex0_uqtk = self.__mindex_uqtk # keep original

        # get projection/ Vandermonde matrix
        self.__Phi_uqtk = uqtkarray.dblArray2D()

        # check if compiled
        self.__compiled = False
        self.compile()

        self.__cv_flag = False

    def compile(self,l_init=0.0,adaptive=0,optimal=1,scale=.1,verbose=0):
        '''
        Setting up variables for the BCS algorithm. Most of the variables do not need to be set. Default settings are sufficient for more cases. See the C++ code for more information about variables.
        '''
        # now we begin BCS routine
        # set work variables
        self.__newmindex_uqtk = uqtkarray.intArray2D() # for uporder iteration
        self.__sigma2_p = uqtkarray.dblArray1D(1) # initial noise variance
        self.__lambda_init = uqtkarray.dblArray1D() # hierarchical prior parameter
        self.__adaptive, self.__optimal, self.__scale, self.__verbose = adaptive,optimal,scale,verbose
        self.__weights_uqtk = uqtkarray.dblArray1D() # weights/ coefficients for basis
        self.__used_uqtk = uqtkarray.intArray1D() # index of weights retained (nonzero)
        self.__errbars_uqtk = uqtkarray.dblArray1D() # error bars for each weight
        self.__nextbasis_uqtk = uqtkarray.dblArray1D() # if adaptive
        self.__alpha_uqtk = uqtkarray.dblArray1D() # prior hyperparameter (1/gamma)
        self.__lambda_p = uqtkarray.dblArray1D(1,l_init)
        print(self.__lambda_p[0])

        #uqtktools.doublep_assign(self.__lambda_p,l_init)
        #self.__lambda_p.assign(0,l_init)

        self.__compiled = True

    def leastsq(self,X,y):
        '''
        perform simple least squares based on the original
        pc order.
        '''
        # convert input to uqtk arrays
        self.__X_uqtk = uqtkarray_tools.numpy2uqtk(X)
        self.__y_uqtk = uqtkarray_tools.numpy2uqtk(y)

        # get vandermonde matrix w.r.t. original pc basis
        self.__V_uqtk = uqtkarray.dblArray2D()
        self.__pcmodel0 = uqtkpce.PCSet("NISPnoq",self.__mindex0_uqtk,self.pctype,0.0,1.0) # initiate
        self.__pcmodel0.EvalBasisAtCustPts(self.__X_uqtk,self.__V_uqtk)
        self.Vandermonde = uqtkarray_tools.uqtk2numpy(self.__V_uqtk)
        self.__sol = np.linalg.lstsq(self.Vandermonde,y)
        return self.__sol[0], self.__sol[1]

    def fit0(self,X,y,tol=1e-8,sigsq=None,upit=0):
        '''
        Train bcs model coefficients with X and y data
        X         :    2d numpy of inputs/ feature data
        y         :     1d numpy array of labels/ outputs
        tol     :    tolerance (smaller means we keep more coefficients)
        sigsq    :     initial noise set automatically based on y data
        upit     :     (int) number of iterations to add higher order terms

        returns the polynomial coefficient (weights), and the mulitindex. One can also return the sensitivity indices by calling self.sens
        '''
        if self.__compiled == False:
            print("Need to compile first!")

        # convert numpy test data into uqtk data types
        self.__X_uqtk = uqtkarray_tools.numpy2uqtk(X)
        self.__y_uqtk = uqtkarray_tools.numpy2uqtk(y)
        self.Xtrain = X
        self.ytrain = y

        if sigsq == None:
            self.__sigma2 = np.var(y)/1e2
        else: self.__sigma2 = sigsq
        #uqtktools.doublep_assign(self.__sigma2_p,self.__sigma2)
        self.__sigma2_p.assign(0,self.__sigma2)
        print(self.__sigma2_p[0])

        self.__tol = tol
        self.__upit = upit

        # begin uporder iterations
        for iter in range(self.__upit+1):

            # get projection/ Vandermonde matrix
            self.__pcmodel = uqtkpce.PCSet("NISPnoq",self.__mindex_uqtk,self.pctype,0.0,1.0) # initiate with new mindex
            self.__pcmodel.EvalBasisAtCustPts(self.__X_uqtk,self.__Phi_uqtk)
            self.__Phi = uqtkarray_tools.uqtk2numpy(self.__Phi_uqtk)

            # resest sigma parameter (if not, may get seg fault)
            #uqtktools.doublep_assign(self.__sigma2_p,self.__sigma2)
            self.__sigma2_p.assign(0,self.__sigma2)

            # change to uqtkbcs.BCS if testing outside source
            bcs.BCS(self.__Phi_uqtk,self.__y_uqtk,self.__sigma2_p,self.__tol,self.__lambda_init,self.__adaptive,self.__optimal,self.__scale,self.__verbose,self.__weights_uqtk,self.__used_uqtk,self.__errbars_uqtk,self.__nextbasis_uqtk,self.__alpha_uqtk,self.__lambda_p)
            print(self.__lambda_p[0])
            # add new mulitindex to newmindex
            uqtkarray.subMatrix_row_int(self.__mindex_uqtk,self.__used_uqtk,self.__newmindex_uqtk)

            if iter < self.__upit :
                # redefine mindex = newmindex if still iterating
                self.__newmindex_added_uqtk = uqtkarray.intArray2D()
                uqtktools.upOrder(self.__newmindex_uqtk,self.__newmindex_added_uqtk)
                self.__mindex_uqtk = self.__newmindex_added_uqtk
                # print("New mindex basis: ", uqtk2numpy(self.__mindex_uqtk)[len(self.__newmindex_uqtk):])

        # return new multiindex to create new pce model
        self.__pcmodel_new = uqtkpce.PCSet("NISPnoq",self.__newmindex_uqtk,self.pctype,0.0,1.0)
        self.mindex = uqtkarray_tools.uqtk2numpy(self.__newmindex_uqtk)
        # eff_dim = self.ndim - sum(sum(self.mindex,0) == 0)

        self.weights = uqtkarray_tools.uqtk2numpy(self.__weights_uqtk)
        self.weight_index = uqtkarray_tools.uqtk2numpy(self.__used_uqtk)
        self.error_bars = uqtkarray_tools.uqtk2numpy(self.__errbars_uqtk)

        # get main effect sensitivity indices
        self.__main_eff_uqtk = uqtkarray.dblArray1D()
        self.__tot_eff_uqtk = uqtkarray.dblArray1D()
        self.__joint_eff_uqtk = uqtkarray.dblArray2D()
        self.__pcmodel_new.ComputeMainSens(self.__weights_uqtk,self.__main_eff_uqtk)
        self.__pcmodel_new.ComputeTotSens(self.__weights_uqtk,self.__tot_eff_uqtk)
        self.__pcmodel_new.ComputeJointSens(self.__weights_uqtk,self.__joint_eff_uqtk)
        self.main_eff = uqtkarray_tools.uqtk2numpy(self.__main_eff_uqtk)
        self.tot_eff = uqtkarray_tools.uqtk2numpy(self.__tot_eff_uqtk)
        self.joint_eff = uqtkarray_tools.uqtk2numpy(self.__joint_eff_uqtk)
        self.senssum = {"main effect": self.main_eff, "total effect": self.tot_eff, "joint effect": self.joint_eff}
        return self.weights, self.mindex

    def predict(self,Xtest,verbose=1):
        '''
        Predict values after training the data
        Xtest     :     2d numpy array

        returns 1d numpy scalar array of predictions
        '''
        if self.__cv_flag == True:
            if verbose == 1:
                print("CV should only be used for hyperparameter tuning. Once tuned, run with nk = 1. ")

        self.__Xtest_uqtk = uqtkarray_tools.numpy2uqtk(Xtest)
        self.__ytest_uqtk = uqtkarray.dblArray1D()
        self.__pcmodel_new.EvalPCAtCustPoints(self.__ytest_uqtk,self.__Xtest_uqtk,self.__weights_uqtk)
        self.__ytest = uqtkarray_tools.uqtk2numpy(self.__ytest_uqtk)
        return self.__ytest
    def getsens(self):
        '''
        return sensitivities as dictionary
        '''
        return self.senssum
    def fit(self,X,y,tol=1e-8,sigsq=None,upit=0,nkfolds=1,seed=13):
        '''
        Fit with k-fold cross validation
        '''
        if nkfolds > 1: self.__cv_flag = True

        # get k-fold cross-validation sets
        self.K = kfoldCV(X,y,nfolds=nkfolds,seed=seed)
        self.__E = []
        self.__C = []
        self.__M = []
        for key in self.K:
            self.Xtrain_temp = self.K[key]['xtrain']
            self.ytrain_temp = self.K[key]['ytrain']
            self.Xval_temp = self.K[key]['xval']
            self.yval_temp = self.K[key]['yval']
            # setup, git and predict bcs model
            self.__c_temp, self.__mindex_temp = self.fit0(self.Xtrain_temp,self.ytrain_temp,upit=upit,tol=tol)
            self.__ypred_temp = self.predict(self.Xval_temp,verbose=0)
            self.__err = np.mean((self.__ypred_temp - self.yval_temp)**2)/np.mean(self.yval_temp**2)
            self.__E.append(self.__err)
            self.__C.append(self.__c_temp)
            self.__M.append(self.__mindex_temp)
        # print("avg NMSE is %f" %np.mean(E))
        return self.__E, self.__C, self.__M
