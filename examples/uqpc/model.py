#!/usr/bin/env python
#=========================================================================
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
#=========================================================================

from __future__ import print_function
import argparse
import os
import sys
try:
    import numpy as np
except ImportError:
    print("Numpy was not found.")

sys.path.append(os.environ['UQTK_INS'])
from PyUQTk.utils.func import *


#############################################################
#############################################################

def model(modelPar, xPar=None, modelname='example'):
    if modelname == 'example':
        return model_example(modelPar)
    elif modelname == 'genz':
        return model_genz(modelPar)
    elif modelname == 'ex_xp':
        return model_xp(modelPar, xPar)
    else:
        return model_example(modelPar)


#############################################################

def model_xp(modelPar, xPar):
    npar = modelPar.shape[0]
    mdim = modelPar.shape[1]
    nout = xPar.shape[0]
    xdim = xPar.shape[1]

    output = np.empty((npar, nout))

    output = (1.5 + xPar.sum(axis=1)) * np.exp(-np.outer(modelPar.sum(axis=1)**2, 2. + xPar.sum(axis=1))/3.)

    return output

#############################################################


def model_genz(modelPar):
    npar = modelPar.shape[0]
    mdim = modelPar.shape[1]
    genzparams = 1. / np.arange(1, mdim + 2)
    # genzparams[1:]=1./(genzparams[1:]+1)

    modelnames = ['genz_gaus', 'genz_exp']

    nout = len(modelnames)

    output = np.empty((npar, nout))
    for j in range(nout):
        output[:, j] = func(modelPar, modelnames[j], genzparams)

    return output

# #############################################################


def model_example(modelPar):
    npar = modelPar.shape[0]
    mdim = modelPar.shape[1]

    nout = 7

    # Create the design parameters #TODO maybe should be done outside
    designPar = np.array(range(1, nout + 1)).reshape(-1, 1)  # /float(nout-1)
    np.savetxt('designPar.dat', designPar)

    output = np.empty((npar, nout))
    for i in range(npar):
        print("%d: Running the model with parameter setting " %
              (i + 1), modelPar[i, :])
        for j in range(nout):
            aa = np.dot(modelPar[i, :] + modelPar[i, :]**2,
                        np.exp(np.log(np.arange(1, mdim + 1)) * (designPar[j])))
            output[i, j] = aa * (sum(modelPar[i, :]))

    return output

#############################################################


def main(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", dest="input_file",
                        type=str, default='ptrain.dat', help="Input file name")
    parser.add_argument("-o", "--output", dest="output_file",
                        type=str, default='ytrain.dat', help="Output file name")
    parser.add_argument("-x", "--xinput", dest="xinput_file",
                        type=str, default='xcond.dat',
                        help="X-conditions (design input) file name; used for modelname='ex_xp'")
    parser.add_argument("-m", "--modelname", dest="model_name",
                        type=str, default='example', help="Model name")
    parser.add_argument("-r", "--noise", dest="noise_stdev",
                        type=float, default=0.0, help="Corrupt with Gaussian noise of this size")
    parser.add_argument("-f", "--factor", dest="factor",
                        type=float, default=1.0, help="Multiply output by a factor")
    args = parser.parse_args()

    modelPar = np.loadtxt(args.input_file, ndmin=2)
    xPar = np.loadtxt(args.xinput_file, ndmin=2)

    output = model(modelPar, xPar=xPar, modelname=args.model_name)
    output *= args.factor
    output += args.noise_stdev * np.random.randn(output.shape[0], output.shape[1])
    np.savetxt(args.output_file, output)


if __name__ == "__main__":
    main(sys.argv[1:])
