#!/usr/bin/env python
#=====================================================================================
#
#                      The UQ Toolkit (UQTk) version 3.1.3
#                          Copyright (2023) NTESS
#                        https://www.sandia.gov/UQToolkit/
#                        https://github.com/sandialabs/UQTk
#
#     Copyright 2023 National Technology & Engineering Solutions of Sandia, LLC (NTESS).
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
#=====================================================================================

import os
import sys
import argparse
import numpy as np
from PyUQTk.utils.mindex_order import getNPC

uqtk_bin = os.environ['UQTK_INS'] + "/bin/"

##########################################################
##########################################################
##########################################################

## Sample physical parameter given its PC coefficient samples and xi's samples
def sample_params(parpc, pctype, order, xis, nsam):

    npar = parpc.shape[1]
    paramdatas = []
    for ip in range(npar):
        # if ip in rndind:
        #     np.savetxt('pcf', parpc[:, ip])
        #     np.savetxt('xdata.dat', xis[:nsam, :])
        #     cmd = 'pce_eval -x PC -f pcf -s ' + pctype + ' -o ' + str(order) + ' > pceval.log'
        #     os.system(cmd)
        #     paramdata = np.loadtxt('ydata.dat')
        # else:
        #     paramdata = parpc[0, ip]

        np.savetxt('pcf', parpc[:, ip])
        np.savetxt('xdata.dat', xis[:nsam, :])
        cmd = uqtk_bin + 'pce_eval -x PC -f pcf -s ' + pctype + ' -o ' + str(order) + ' > pceval.log'
        os.system(cmd)
        paramdata = np.loadtxt('ydata.dat')

        paramdatas.append(paramdata)

    # if xis.shape[0]>nsam:
    #     xis=xis[nsam:,:]
    # else:
    #     print "More xis needed!"
    #     sys.exit()

    #print(xis.shape[0], " xis left")
    return np.array(paramdatas).reshape(npar, -1)


######################################################################
######################################################################
######################################################################

## Parse the input arguments
usage_str='Script to sample posterior physical parameters given posterior PC coefficients.'

parser = argparse.ArgumentParser(description=usage_str)
parser.add_argument("-r", "--rdim", dest="rdim", type=int, default=0,
                    help="Number of embedded dimensions")
parser.add_argument("-o", "--order", dest="order", type=int, default=0,
                    help="Order of embedded PC output")
parser.add_argument("-d", "--npar_phys", type=int, dest="npar_phys", default=None,
                    help="Number of physical inputs")
parser.add_argument("-t", "--pctype", dest="pctype", type=str, default="HG",
                    help="Embedded PC type, only relevant if rdim>0")
parser.add_argument("-x", "--nmaxxi", type=int, dest="nmax_xi", default=0,
                    help="Sample xi's or read them from rvar.dat")
parser.add_argument("-s", "--nsamxi", type=int, dest="nsam_xi", default=1,
                    help="Number of xi samples per posterior sample")
parser.add_argument("-e", "--every", type=int, dest="every", default=1,
                    help="Use every n-th sample")
parser.add_argument("-f", "--fixindnom", dest="fixindnom_file", type=str, default=None,
                    help="Indices and nominals of fixed parameters, if any")
# parser.add_argument("-z", "--sigma_pdf", dest="sigma_pdf", action="store_true", default=False,
#                     help="Whether sigma is the last column or not")
# parser.add_argument("-v", "--verbosity", help="increase output verbosity")
args = parser.parse_args()


rdim = args.rdim  # rndind.shape[0]
order = args.order  # double check - this is -t of model_inf,right? Should match with parampccfs.dat
npar_phys = args.npar_phys
pctype = args.pctype
nsam_xi = args.nsam_xi
nmax_xi = args.nmax_xi
every = args.every
fixindnom_file = args.fixindnom_file

# Need better treatment of including sigma
# pchain = np.loadtxt('pchain.dat', ndmin=2)[::every]
# sigma_pdf = args.sigma_pdf
# if sigma_pdf:
#     post_sigma = pchain[:, -1]
#     np.savetxt('post_sigma.dat', post_sigma)
# print("Number of chain parameters       : ", pchain.shape[1])
# if sigma_pdf:
#     print("Including sigma")

# Read PC coefficients' chain samples
parampccfs = np.loadtxt('parampccfs.dat', ndmin=2)
npc = int(parampccfs.shape[0] / npar_phys)
assert(npc == getNPC(rdim, order))

print("Number of PC terms per parameter : ", npc)

tmp = parampccfs.reshape(npc, npar_phys, -1)
parpc_mcmc = tmp[:, :, :-1][:, :, ::every]
parpc_mean = np.average(parpc_mcmc, axis=2)
npost = parpc_mcmc.shape[2]
parpc_map = tmp[:, :, -1]


# Sample xi's, if embedded case, i.e. rdim>0
if rdim > 0:
    if nmax_xi > 0:
        print("Sampling xis")
        cmd = uqtk_bin + '/pce_rv -w PCvar -d ' + str(rdim) + ' -n ' + str(nmax_xi) + \
                         ' -p ' + str(rdim) + ' -x ' + pctype + ' > pcerv.log'
        print("Running ", cmd)
        os.system(cmd)
    else:
        print("Reading pre-sampled rvar.dat")

    xis = np.loadtxt('rvar.dat').reshape(-1, rdim)
    assert(xis.shape[0] >= (npost + 1) * nsam_xi)  # npost is thinned pchain size

    paramdata_map = sample_params(parpc_map, pctype, order, xis, nsam_xi)
    xis = xis[nsam_xi:, :]

    paramdata_mcmc = np.empty((npar_phys, 0))
    for j in range(npost):
        if j % 100 == 0:
            print("Posterior sample ", j, " out of ", npost)

        paramdata_mcmc = np.hstack((paramdata_mcmc,
                                    sample_params(parpc_mcmc[:, :, j], pctype, order, xis, nsam_xi)))
        if xis.shape[0] > nsam_xi:
            xis = xis[nsam_xi:, :]
        else:
            print("Need more xis")
            sys.exit()

else:
    paramdata_map = parpc_map[0, :].reshape(npar_phys, -1)
    paramdata_mcmc = parpc_mcmc[0, :, :]

paramdata_mcmc = paramdata_mcmc.T

# print(paramdata_map.shape)   # npar_phys, nsam_xi
# print(paramdata_mcmc.shape)  # npost*nsam_xi, npar_phys

# If some parameters are fixed at their nominal, take that into account
if fixindnom_file is not None:
    fixindnom = np.loadtxt(fixindnom_file, ndmin=2)

    npar_phys += fixindnom.shape[0]
    map_phparam = np.empty((npar_phys, nsam_xi))
    post_phparam = np.empty((npost * nsam_xi, npar_phys))
    j = 0
    for ipar in range(npar_phys):
        if ipar in fixindnom[:, 0]:
            nom = fixindnom[np.where(fixindnom[:, 0] == ipar), 1]
            post_phparam[:, ipar] = nom
            map_phparam[ipar, :] = nom
        else:
            post_phparam[:, ipar] = paramdata_mcmc[:, j]
            map_phparam[ipar, :] = paramdata_map[j, :]
            j += 1
else:
    map_phparam = paramdata_map.copy()
    post_phparam = paramdata_mcmc.copy()

# Save to outputs
np.savetxt('post_phparam.dat', post_phparam)
np.savetxt('map_phparam.dat', map_phparam)


