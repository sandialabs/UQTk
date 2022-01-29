#!/usr/bin/env python
#=====================================================================================
#
#                      The UQ Toolkit (UQTk) version 3.1.2
#                          Copyright (2022) NTESS
#                        https://www.sandia.gov/UQToolkit/
#                        https://github.com/sandialabs/UQTk
#
#     Copyright 2022 National Technology & Engineering Solutions of Sandia, LLC (NTESS).
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

"""
Generic tools for evaluation of standard functions
and their integrals
"""

try:
    import numpy as np
except ImportError:
    "Need numpy"

import sys,os
from math import *
import random as rnd
import itertools


from .mindex_order import getNPC

###################################################################################################


def func(xdata,model,func_params):
    """Generic function evaluator.
    Note:
        * Note that conventional Genz arguments are in [0,1], here the expected input is on [-1,1]
               The shifts are applied to the transformed coordinates in [0,1], to be consistent
               with traditional formulations of Genz functions.
               In future versions of UQTk, we may redefine the Genz functions here to take
               arguments in [0,1]. Need to check if that will create trouble with the other
               function types defined here. Maybe take Genz functions out and make them
               a separate functionality?
    Arguments:
        * xdata       : Nxd numpy array of input, should be in [-1,1]^d
        * model       : Model name, options are 'genz_osc', 'genz_exp', 'genz_cont', 'genz_gaus',
                        'genz_cpeak', 'genz_ppeak', 'ishigami', 'sobol', 'poly_exsens', 'PCmi', 'spp'
        * func_params : Auxiliary parameters
                      : For genz functions, an array of size d+1, the first entry being the shift,
                      : which is the same for all dimensions.
                      : and the rest of the entries are the weights.
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
        print('The shape of ww is: ',ww.shape)
        print('The shape of gcf is: ',gcf.shape)
        print('The shape of x is: ',xdata.shape)
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

    elif model == 'PCmi':
        mindex=func_params[0]
        pccf=func_params[1]
        pctype=func_params[2]

        npc=mindex.shape[0]

        #print "mindex=",mindex, pccf

        np.savetxt('xdata.dat',xdata)

        np.savetxt('pccf.dat',pccf)
        np.savetxt('mi.dat',mindex,fmt='%d')
        assert npc==pccf.shape[0]


        # Evaluate PC to get the training data
        if pctype=='MON':
            ydata=np.zeros((sam,))
            for ipt in range(sam):
                for i in range(npc):
                    cur=pccf[i]
                    for j in range(dim):
                        cur*=(xdata[ipt,j]**mindex[i,j])
                    ydata[ipt]+=cur
        else:
            cmd='pce_eval -x"PC_mi" -f"pccf.dat" -s'+pctype+' -r"mi.dat" > pceval.log'
            print('Running %s'%(cmd))
            os.system(cmd)
            ydata=np.loadtxt('ydata.dat')

    elif model == 'spp':

        ord=func_params[0]
        sp=func_params[1]
        pctype=func_params[2]

        np.savetxt('xdata.dat',xdata)

        npc=getNPC(dim,ord)
        #print "npc = ", npc
        # All coefs set to zero
        sppc=np.zeros((npc))

        # Randomly select some coefficients to be equal to 1
        nzind=np.random.choice(npc, sp,replace=False)

        for i in nzind:
            sppc[i]=1.0
        # Save the 'true' coefficients
        np.savetxt('sppc.dat',sppc)

        # Evaluate PC to get the training data
        cmd='pce_eval -x"PC" -o' + str(ord) + ' -f"sppc.dat" -s'+pctype+' > pceval.log'
        os.system(cmd)
        ydata=np.loadtxt('ydata.dat')

    elif model == 'currinExp':
    	# Currin, C., Mitchell, T., Morris, M., & Ylvisaker, D. (1988).
    	# A Bayesian approach to the design and analysis of computer experiments.
    	# Technical Report 6498. Oak Ridge National Laboratory.
    	# https://www.sfu.ca/~ssurjano/curretal88exp.html (March 2018)
    	xdata=0.5*(xdata+1.)
    	f1 = 1 - np.exp(-1/(2*xdata[:,1]));
    	f2 = 2300*xdata[:,0]**3 + 1900*xdata[:,0]**2 + 2092*xdata[:,0] + 60.0;
    	f3 = 100*xdata[:,0]**3 + 500*xdata[:,0]**2 + 4*xdata[:,0] + 20.0;
    	ydata = f1 * f2 / f3
    	return(ydata)

    elif model == 'currinExpLF':
        # see above comment on 'currinExp', and below
        # Xiong, S., Qian, P. Z., & Wu, C. J. (2013). Sequential design and analysis
        # of high-accuracy and low-accuracy computer codes. Technometrics, 55(1), 37-46
        xdataLF = xdata.copy()+0.1
        y1 = func(xdataLF,'currinExp',[]);
        for i in range(sam):
            xdataLF[i,1] = max(0.0,xdataLF[i,1]-0.2)
        y2 = func(xdataLF,'currinExp',[])
        xdataLF      = xdata.copy()-0.1
        xdataLF[:,1] = xdataLF[:,1]+0.2
        y3 = func(xdataLF,'currinExp',[])
        for i in range(sam):
            xdataLF[i,1] = max(0.0,xdataLF[i,1]-0.2)
        y4 = func(xdataLF,'currinExp',[])
        ydata = 0.25 * (y1 + y2 + y3 + y4)
        return ydata

    elif model == 'park91F1':
        # Park, J.S. (1991). Tuning complex computer codes to data and optimal designs.
        # Ph.D. Thesis, University of Illlinois, Champaign-Urbana.
        # https://www.sfu.ca/~ssurjano/park91a.html (March 2018)
        xdata=0.5*(xdata+1.)
        f1 = 0.5*xdata[:,0]*(np.sqrt(1.0 + (xdata[:,1]+xdata[:,2]**2)*xdata[:,3]/(xdata[:,0]**2)) - 1.0)
        f2 = (xdata[:,0] + 3*xdata[:,3]) * np.exp(1.0 + np.sin(xdata[:,2]))
        ydata = f1 + f2
        return ydata

    elif model == 'park91F1LF':
        # see above comment on 'park91F1', and below
        # Xiong, S., Qian, P. Z., & Wu, C. J. (2013). Sequential design and analysis
        # of high-accuracy and low-accuracy computer codes. Technometrics, 55(1), 37-46
        yHF=func(xdata,'park91F1',[])
        xdata=0.5*(xdata+1.)
        f1 = (1.0+0.1*np.sin(xdata[:,0])) * yHF
        f2 = -2*xdata[:,0] + xdata[:,1]**2 + xdata[:,2]**2
        ydata = f1 + f2 + 0.5;
        return ydata

    elif model == 'park91F2':
        # Park, J.S. (1991). Tuning complex computer codes to data and optimal designs.
        # Ph.D. Thesis, University of Illlinois, Champaign-Urbana.
        # https://www.sfu.ca/~ssurjano/park91a.html (March 2018)
        xdata=0.5*(xdata+1.)
        ydata = 2.0/3.0 * np.exp(xdata[:,0]+xdata[:,1]) - xdata[:,3] * np.sin(xdata[:,2]) + xdata[:,2]
        return ydata

    elif model == 'park91F2LF':
        # see above comment on 'park91F2', and below
        # Xiong, S., Qian, P. Z., & Wu, C. J. (2013). Sequential design and analysis
        # of high-accuracy and low-accuracy computer codes. Technometrics, 55(1), 37-46
        ydata = 1.2*func(xdata,'park91F2',[])-1.0
        return ydata

    else:
        print('Function type is not recognized. Exiting.')
        sys.exit()


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
        print('No exact sensitivity available for this function. Exiting.')
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
