#!/usr/bin/env python
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

"""
Generic tools for evaluation of standard functions
and their integrals
"""

try:
    import numpy as np
except ImportError:
    "Need numpy"

import sys
from math import *
import random as rnd
import itertools

###################################################################################################


def func(xdata,model,func_params):
    """Generic function evaluator.
    Note:
        * Currently only genz functions are implemented.
        * Note that conventional Genz arguments are in [0,1], here the expected input is on [-1,1]
    Arguments:
        * xdata       : Nxd numpy array of input, should be in [-1,1]^d
        * model       : Model name, options are 'genz_osc', 'genz_exp', 'genz_cont', 'genz_gaus', 'genz_cpeak', 'genz_ppeak'
        * func_params : Auxiliary parameters
                      : For genz functions, an array of size d+1, the first entry being the shift, and the rest of the entries are the weights.
                      : See UQTk Manual for Genz formulae.
    Returns:
        * ydata       : An array of outputs of size N.
    """

    # Get the input size
    sam=xdata.shape[0]
    dim=xdata.shape[1]

    # Check the function types and evaluate
    if model == 'genz_osc':
        xdata=0.5*(xdata+1.)
        ydata=np.empty((sam,))
        gcf=func_params[1:]
        xtmp=np.dot(xdata,gcf)
        for j in range(sam):
            ydata[j]=cos(2.*pi*func_params[0]+xtmp[j])

    elif model == 'genz_exp':
        xdata=0.5*(xdata+1.)
        ydata=np.empty((sam,))
        ww=func_params[0]
        gcf=func_params[1:]

        xtmp=np.dot(xdata-ww,gcf)
        for j in range(sam):
            ydata[j]=exp(xtmp[j])

    elif model == 'genz_cont':
        xdata=0.5*(xdata+1.)
        ydata=np.empty((sam,))
        ww=func_params[0]
        gcf=func_params[1:]

        xtmp=np.dot(abs(xdata-ww),gcf)
        for j in range(sam):
            ydata[j]=exp(-xtmp[j])

    elif model == 'genz_gaus':
        xdata=0.5*(xdata+1.)
        ydata=np.empty((sam,))
        ww=func_params[0]
        gcf=func_params[1:]

        xtmp=np.dot((xdata-ww)*(xdata-ww),gcf*gcf)
        for j in range(sam):
            ydata[j]=exp(-xtmp[j])

    elif model == 'genz_cpeak':
        xdata=0.5*(xdata+1.)
        ydata=np.empty((sam,))
        #ww=param[0]
        gcf=func_params[1:]

        xtmp=1.+(np.dot(xdata,gcf)) #use abs if defined on [-1,1]
        for j in range(sam):
            ydata[j]=exp(-(dim+1.)*log(xtmp[j]))

    elif model == 'genz_ppeak':
        xdata=0.5*(xdata+1.)
        ydata=np.empty((sam,))
        ww=func_params[0]
        gcf=func_params[1:]

        for j in range(sam):
            prod=1.
            for i in range(dim):
                prod = prod / (1./(gcf[i]**2.)+(xdata[j,i]-ww)**2.)
            ydata[j]=prod

    elif model == 'ishigami':
        assert(dim==3)
        a=func_params[0]
        b=func_params[1]
        ydata=np.empty((sam,))

        for j in range(sam):
            ydata[j]=np.sin(xdata[j,0])+a*np.sin(xdata[j,1])**2+b*np.sin(xdata[j,0])*xdata[j,2]**4

    elif model == 'sobol':
        assert(dim==func_params.shape[0])
        ydata=np.empty((sam,))
        for j in range(sam):
            val=1.
            for k in range(dim):
                val *= ( (abs(2*xdata[j,k])+func_params[k])/(1.+func_params[k]) )
            ydata[j]=val

    elif model == 'poly_exsens':
        assert(dim==func_params[0])
        ydata=np.empty((sam,))
        for j in range(sam):
            val=1.
            for k in range(dim):
                val *= ( (3./4.)*(xdata[j,k]+1.)**2+1. )/2.
            ydata[j]=val

    return ydata

##################################################################################

def integ_exact(model,func_params):
    """Analytically available function integrals.
    Note:
        * Currently only genz functions are implemented.
        * Note that conventional Genz arguments are in [0,1], here the expected input is on [-1,1]
    Arguments:
        * model       : Model name, options are 'genz_osc', 'genz_exp', 'genz_cont', 'genz_gaus', 'genz_cpeak', 'genz_ppeak'
        * func_params : Auxiliary parameters
                      : For genz functions, an array of size d+1, the first entry being the shift, and the rest of the entries are the weights.
                      : See UQTk Manual for Genz integral formulae.
    Returns:
        * integ_ex    : A real number that is the integral over [-1,1]^d
    """


    if (model=='genz_osc'):
        gcf=func_params
        dim=gcf.shape[0]-1
        integ_ex=cos(2.*pi*gcf[0]+0.5*sum(gcf[1:]))
        for i in range(1,dim+1):
            integ_ex*=(2.*sin(gcf[i]/2.)/gcf[i])
    elif (model=='genz_exp'):
        gcf=func_params
        dim=gcf.shape[0]-1
        integ_ex=1.
        for i in range(1,dim+1):
            at1=exp(-gcf[i]*gcf[0])
            at2=exp(gcf[i]*(1.-gcf[0]))
            integ_ex*=((at2-at1)/(gcf[i]))
    elif (model=='genz_cont'):
        gcf=func_params
        dim=gcf.shape[0]-1
        integ_ex=1.
        for i in range(1,dim+1):
            integ_ex*= ((2.-exp(gcf[i]*(-gcf[0]))-exp(gcf[i]*(gcf[0]-1.)))/gcf[i])
    elif (model=='genz_gaus'):
        gcf=func_params
        dim=gcf.shape[0]-1
        integ_ex=1.
        for i in range(1,dim+1):
            at1=erf(-gcf[i]*gcf[0])
            at2=erf(gcf[i]*(1.-gcf[0]))
            integ_ex*=((at2-at1)*sqrt(pi)/(2.*gcf[i]))
    elif (model=='genz_cpeak'):
        gcf=func_params
        dim=gcf.shape[0]-1
        numer=0.0
        count=1
        denom=1.
        for i in range(1,dim+1):
            comb=list(itertools.combinations(range(1,dim+1),i))
            for j in range(len(comb)):
                assert(i==len(comb[j]))
                #print i,j,pow(-1,i)
                numer+=(pow(-1,i)/(1.+sum(gcf[list(comb[j])])))
                count+=1
            denom*=(i*gcf[i])
        #print count, numer
        integ_ex=(1.+numer)/denom
    elif (model=='genz_ppeak'):
        gcf=func_params
        dim=gcf.shape[0]-1
        integ_ex=1.
        for i in range(1,dim+1):
            at1=np.arctan(-gcf[i]*gcf[0])
            at2=np.arctan(gcf[i]*(1.-gcf[0]))
            integ_ex*=(gcf[i]*(at2-at1))

    return integ_ex

################################################################################
################################################################################

def mainsens_exact(model,func_params):
    """Analytically available main sensitivities for some functions.
    Note:
        * Currently only sobol, ishigami and poly_exsens functions are implemented.
        * Note that conventional sobol arguments are in [0,1], here the expected input is on [-1,1]
    Arguments:
        * model       : Model name, options are 'sobol', 'ishigami', 'poly_exsens'
        * func_params : Auxiliary parameters
    Returns:
        * mainsens    : Main effect Sobol sensitivity index
    """    
    if (model=='sobol'):
        dim=func_params.shape[0]
        mainsens=np.empty((dim,))
        var=1.0
        for i in range(dim):
            mainsens[i]=1./(3.*(1.+func_params[i])**2)
            var*=(mainsens[i]+1.)
        var-=1.0
        mainsens/=var

    elif (model=='ishigami'):
        a=func_params[0]
        b=func_params[1]
        var=a**2/8.+b*np.pi**4/5.+b**2*np.pi**8/18.+0.5
        mainsens=np.empty((3,))
        mainsens[0]=b*np.pi**4/5.+b**2*np.pi**8/50.+0.5
        mainsens[1]=a**2/8.
        mainsens[2]=0.0
        mainsens/=var

    elif (model=='poly_exsens'):
        dim=func_params[0]
        mainsens=(0.2/(1.2**dim-1))*np.ones((dim,))

    else:
        print "No exact sensitivity available for this function. Exiting."
        sys.exit(1)


    return mainsens

##################################################################################
##################################################################################

def main(arg):
    modelname=arg[0]
    input_file=arg[1]
    output_file=arg[2]
    auxparam=[]
    if len(arg)>3:
        auxparam_file=arg[3]
        auxparam=np.loadtxt(auxparam_file,ndmin=1) 

    input=np.loadtxt(input_file,ndmin=2) 

    output=func(input,modelname,auxparam)
    np.savetxt(output_file,output)

if __name__ == "__main__":
    main(sys.argv[1:])
