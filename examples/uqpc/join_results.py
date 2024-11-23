#!/usr/bin/env python
#=====================================================================================
#
#                      The UQ Toolkit (UQTk) version 3.1.5
#                          Copyright (2024) NTESS
#                        https://www.sandia.gov/UQToolkit/
#                        https://github.com/sandialabs/UQTk
#
#     Copyright 2024 National Technology & Engineering Solutions of Sandia, LLC (NTESS).
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
#=====================================================================================

################################################################################
# Example script to join results pickle files from various tasks
################################################################################

# Import libraries

import os
import sys
import glob

if (sys.version_info.major==2):
    import cPickle as pick
elif (sys.version_info.major==3):
    import pickle as pick
else:
    print("Only Python 2 or 3 are supported. Exiting.")
    sys.exit()

try:
    import numpy as np
except ImportError:
    print("Numpy module could not be found")



# Searching in task_* directories
if len(glob.glob('task_*/results.pk'))==0:
    print("join_results.py: Files task_*/results.pk do not exist. Exiting.")

# First pass to get the sizes
nout=0
kk=0
for filename in glob.iglob('task_*/results.pk'):
    if (not os.path.isfile(filename)):
        print("join_results.py: %s does not exist. Exiting." % filename)
        sys.exit()
    results=pick.load(open(filename, 'rb'))

    outrange=results['outs']
    if (nout<max(outrange)):
        nout=max(outrange)
    if kk==0:
        inqdp,inpar,ytrain,ytrain_pc,errcheck_pc=results['training']
        npt=inpar.shape[0]
        npar=inpar.shape[1]

        if 'validation' in results.keys():
            qpar_val,inpar_val,yval,yval_pc,errcheck_val_pc=results['validation']
            nval=inpar_val.shape[0]

nout+=1
outrange_all=range(nout)
pccf_all_all=[0]*nout
mindex_all_all=[0]*nout
varfrac_all_all=[0]*nout
ccov_all_all=[0]*nout
allsens_main_all=np.empty((nout,npar))
allsens_total_all=np.empty((nout,npar))
allsens_joint_all=np.empty((nout,npar,npar))
ytrain_pc_all=np.empty((npt,nout))
yval_pc_all=np.empty((nval,nout))
errcheck_pc_all=np.empty((npt,nout))
errcheck_val_pc_all=np.empty((nval,nout))

err_training_all=np.empty((nout,))
err_val_all=np.empty((nout,))


val=False
# Second pass to actually collect data
for filename in glob.iglob('task_*/results.pk'):

    results=pick.load(open(filename, 'rb'))
    print("Collecting from %s : output # %d" % (filename,results['outs']))


    if 'validation' in results.keys():
        val=True

    outrange=results['outs']
    inqdp,inpar,ytrain,ytrain_pc,errcheck_pc=results['training']
    if (val):
        qpar_val,inpar_val,yval,yval_pc,errcheck_val_pc=results['validation']
        err_training,err_val=results['err']
    else:
        err_training=results['err']

    pccf_all,mindex_all,varfrac_all,ccov_all,pc_type=results['pcmi']
    allsens_main,allsens_total,allsens_joint=results['sens']


    i=0
    #print "Adding ",mindex_all,"to",mindex_all_all

    for ind in outrange:
        pccf_all_all[ind]=pccf_all[i]
        mindex_all_all[ind]=mindex_all[i]
        varfrac_all_all[ind]=varfrac_all[i]
        ccov_all_all[ind]=ccov_all[i]

        err_training_all[ind]=err_training[i]
        ytrain_pc_all[:,ind]=ytrain_pc[:,i]
        errcheck_pc_all[:,ind]=errcheck_pc[:,i]
        allsens_main_all[ind,:]=allsens_main[i,:]
        allsens_total_all[ind,:]=allsens_total[i,:]
        allsens_joint_all[ind,:,:]=allsens_joint[i,:,:]

        if (val):
            err_val_all[ind]=err_val[i]
            yval_pc_all[:,ind]=yval_pc[:,i]
            errcheck_val_pc_all[:,ind]=errcheck_val_pc[:,i]

        i+=1



# Write out to a general pickle file
if (val):
    results_all = {'outs':(outrange_all),'training':(inqdp,inpar,ytrain,ytrain_pc_all,errcheck_pc_all),'validation':(qpar_val,inpar_val,yval,yval_pc_all,errcheck_val_pc_all),'pcmi':(pccf_all_all,mindex_all_all,varfrac_all_all,ccov_all_all,pc_type),'sens':(allsens_main_all,allsens_total_all,allsens_joint_all),'err':(err_training_all,err_val_all)}
else:
    results_all = {'outs':(outrange_all),'training':(inqdp,inpar,ytrain,ytrain_pc_all,errcheck_pc_all),'pcmi':(pccf_all_all,mindex_all_all,varfrac_all_all,ccov_all_all,pc_type),'sens':(allsens_main_all,allsens_total_all,allsens_joint_all),'err':(err_training_all)}


pick.dump(results_all,open('results.pk','wb'),-1)



