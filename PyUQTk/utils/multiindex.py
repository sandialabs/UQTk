#!/usr/bin/env python
#=====================================================================================
#
#                      The UQ Toolkit (UQTk) version 3.1.1
#                          Copyright (2021) NTESS
#                        https://www.sandia.gov/UQToolkit/
#                        https://github.com/sandialabs/UQTk
#
#     Copyright 2021 National Technology & Engineering Solutions of Sandia, LLC (NTESS).
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
Scripts for managing multiindices.
"""

import os
import sys

try:
    import numpy as np
except ImportError:
    print('Numpy was not found.')

#############################################################
#############################################################

def gen_mi(mi_type,params):
    """
    Wrapper around the app gen_mi for generating multiindex sets
    Arguments:
        * mi_type : Multiindex tpye, options are 'TO', 'TP', 'HDMR'
        * params  : Parameters, a two-element tuple
                  : First element is the order ('TO'), list of orders per dimension ('TP'), or list of HDMR orders ('HDMR')
                  : Second element is dimensionality
    Returns:
        * mindex  : A 2d array of multiindices.
    """

    # Total-Order truncation
    if mi_type=='TO':
        # Order
        nord=params[0]
        # Dimensionality
        dim=params[1]
        # Command for the app
        cmd='gen_mi -x' + mi_type + ' -p' + str(nord) + ' -q' + str(dim)

    # Tensor-product truncation
    elif mi_type=='TP':
        # A list of orders per dimension
        orders=params[0]
        # Dimensionality
        dim=params[1]
        assert(dim==len(orders))
        # Save the per-dimension orders in a file
        np.savetxt('orders.dat',np.array(orders),fmt='%d')
        # Command for the app
        cmd='gen_mi -x' + mi_type + ' -f orders.dat -q'+str(dim)

    # HDMR trunction
    elif mi_type=='HDMR':
        # A list of per-variate orders
        hdmr_dims=params[0]
        # Dimensionality
        dim=params[1]
        # Save the HDMR dimensions in a file
        np.savetxt('hdmr_dims.dat',np.array(hdmr_dims),fmt='%d')
        # Command for the app
        cmd='gen_mi -x' + mi_type + ' -f hdmr_dims.dat -q'+str(dim)

    else:
        print('Multiindex type is not recognized. Use \'TO\', \'TP\' or \'HDMR\'. Exiting.')
        sys.exit(1)

    # Run the app
    os.system(cmd + ' > gen_mi.out')

    # Load the generated multtindex file
    mindex=np.loadtxt('mindex.dat',dtype=int).reshape(-1,dim)
    return mindex


#############################################################
#############################################################

def mi_addfront_cons(mindex):
    """
    Adding a front to multiindex in a conservative way, i.e.
    a multiindex is added only if *all* parents are in the current set
    """

    print('Adding multiindex front (conservative)')

    npc=mindex.shape[0]
    ndim=mindex.shape[1]
    mindex_f=np.zeros((1,ndim),dtype=int)
    mindex_add=np.zeros((1,ndim),dtype=int)
    mindex_new=np.zeros((1,ndim),dtype=int)
    for i in range(npc):
        cur_mi=mindex[i,:]

        fflag=True
        for j in range(ndim):
            test_mi=np.copy(cur_mi)
            test_mi[j] += 1
            #print "Trying test_mi", test_mi
            fl=True


            if not any(np.equal(mindex,test_mi).all(1)):
                for k in range(ndim):
                    if(test_mi[k]!=0):
                        subt_mi=np.copy(test_mi)
                        subt_mi[k] -= 1

                        if any(np.equal(mindex,subt_mi).all(1)):
                            cfl=True
                            fl=cfl*fl

                        else:
                            fl=False
                            break


                if (fl):
                    if not any(np.equal(mindex_add,test_mi).all(1)):
                        mindex_add=np.vstack((mindex_add,test_mi))
                    if fflag:
                        mindex_f=np.vstack((mindex_f,cur_mi))
                    fflag=False

    mindex_f=mindex_f[1:]
    mindex_add=mindex_add[1:]
    mindex_new=np.vstack((mindex,mindex_add))

    print('Multiindex resized from %d to %d.'%(mindex.shape[0],mindex_new.shape[0]))

    # Returns the new muliindex, the added new multiindices,
    # and the 'front', i.e. multiindices whose children are added
    return [mindex_new,mindex_add,mindex_f]

#############################################################
#############################################################
#############################################################

def mi_addfront(mindex):
    """
    Adding a front to multiindex in a non-conservative way, i.e.
    a multiindex is added only if *any* of the parents is in the current set
    """

    print('Adding multiindex front (non-conservative)')

    npc=mindex.shape[0]
    ndim=mindex.shape[1]

    mindex_f=np.zeros((1,ndim),dtype=int)
    mindex_add=np.zeros((1,ndim),dtype=int)
    mindex_new=np.zeros((1,ndim),dtype=int)
    for i in range(npc):
        cur_mi=mindex[i,:]

        fflag=True
        for j in range(ndim):
            test_mi=np.copy(cur_mi)
            test_mi[j] += 1
            if not any(np.equal(mindex,test_mi).all(1)):
                if not any(np.equal(mindex_add,test_mi).all(1)):
                    mindex_add=np.vstack((mindex_add,test_mi))
                if fflag:
                    mindex_f=np.vstack((mindex_f,cur_mi))
                fflag=False

    mindex_f=mindex_f[1:]
    mindex_add=mindex_add[1:]
    mindex_new=np.vstack((mindex,mindex_add))


    print('Multiindex resized from %d to %d.'%(mindex.shape[0],mindex_new.shape[0]))

    # Returns the new muliindex, the added new multiindices,
    # and the 'front', i.e. multiindices whose children are added
    return [mindex_new,mindex_add,mindex_f]



#####################################################################
