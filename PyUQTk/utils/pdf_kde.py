#!/usr/bin/env python
#=====================================================================================
#                     The UQ Toolkit (UQTk) version @UQTKVERSION@
#                     Copyright (@UQTKYEAR@) Sandia Corporation
#                     http://www.sandia.gov/UQToolkit/
#
#     Copyright (@UQTKYEAR@) Sandia Corporation. Under the terms of Contract DE-AC04-94AL85000
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
KDE PDF computation routines.
"""

import os
try:
    from scipy import stats
except ImportError:
    print('Scipy was not found.')

try:
    import numpy as np
except ImportError:
    print('Numpy was not found.')


#############################################################

def get_pdf(data,target, method='UQTk',verbose=1):
    """
    Compute PDF given data at target points
    with Python built-in method or with the UQTk app

    Arguments:
        * data   : an N x d array of N samples in d dimensions
        * target : an M x d array of target points
                 : can be an integer in method-UQTk case; and is interpreted as
                 : the number of grid points per dimension for a target grid
        * method : 'UQTk' or 'Python'
        * verbose: verbosity on the screen, 0,1, or 2

    Returns:
        * xtarget : target points (same as target, or a grid, if target is an integer)
        * dens    : PDF values at xtarget
    """
    np.savetxt('data',data)

    # Wrapper around the UQTk app
    if (method=='UQTk'):

        if (verbose>1):
            outstr=''
        else:
            outstr=' > pdfcl.log'

        if(type(target)==int):
            cmd='pdf_cl -i data -g '+str(target)+outstr
            if (verbose>0):
                print('Running %s'&(cmd))
            os.system(cmd)

        else:
            np.savetxt('target',target)
            cmd='pdf_cl -i data -x target'+outstr
            if (verbose>0):
                print('Running %s' % (cmd))

            os.system(cmd)

        xtarget=np.loadtxt('dens.dat')[:,:-1]
        dens=np.loadtxt('dens.dat')[:,-1]

    # Python Scipy built-in method of KDE
    elif (method=='Python'):
        assert (type(target)!=int)
        np.savetxt('target',target)

        kde_py=stats.kde.gaussian_kde(data.T)
        dens=kde_py(target.T)
        xtarget=target

    else:
        print('KDE computation method is not recognized (choose \'Python\' or \'UQTk\'). Exiting.')
        sys.exit()

    # Return the target points and the probability density
    return xtarget,dens

