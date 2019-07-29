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
#=====================================================================================

################################################################################
# Input PC generation given samples or marginal pc
################################################################################



# Input
# Usage  : prepare_inpc.py <format> <filename> <input_pcorder>
# e.g.   : prepare_inpc.py marg marg_pc.txt 3
#        : prepare_inpc.py sam inp_sam.txt 3
# Output : param_pcf.txt (this is in format that can be used in with the uq_pc.py -c flag)


import os
import sys

try:
	import numpy as np
except ImportError:
	print("Numpy module could not be found")


uqtkbin=os.environ['UQTK_INS']+"/bin/"

input_format=sys.argv[1]
filename=sys.argv[2]
input_pcorder=int(sys.argv[3])


if input_format=="marg":

	with open(filename) as f:
	    margpc = f.readlines()
	dim=len(margpc)



	margpc_all=[]
	maxord=0
	sumord=0
	for i in range(dim):
		margpc_cur=np.array(margpc[i].split(),dtype=float)
		order=margpc_cur.shape[0]-1
		if maxord<order:
			maxord=order
		sumord+=order
		margpc_all.append(margpc_cur)

	# Command for the app
	assert(input_pcorder >= maxord)
	cmd=uqtkbin+'gen_mi -x TO -p' + str(input_pcorder) + ' -q' + str(dim)+' > genmi.log'
	os.system(cmd)
	mindex_totalorder=np.loadtxt('mindex.dat',dtype=int).reshape(-1,dim)


	mindex=np.zeros((1,dim),dtype=int)
	cfs=np.zeros((1,dim))
	for i in range(dim):
		cfs[0,i]=margpc_all[i][0]

	for i in range(dim):

		order_this=margpc_all[i].shape[0]-1
		mindex_this=np.zeros((order_this,dim),dtype=int)
		mindex_this[:,i]=np.arange(1,order_this+1)
		mindex=np.vstack((mindex,mindex_this))

		cfs_this=np.zeros((order_this,dim))
		cfs_this[:,i]=margpc_all[i][1:]
		cfs=np.vstack((cfs,cfs_this))

	k=0
	cfs_totalorder=np.zeros(mindex_totalorder.shape)
	for mi in mindex_totalorder:
		if any(np.equal(mi,mindex).all(1)):
			ind=np.equal(mi,mindex).all(1).tolist().index(True)
			cfs_totalorder[k,:]=cfs[ind,:]
		k+=1

	np.savetxt('param_pcf.txt',cfs_totalorder)


elif input_format=="sam":

	cmd=uqtkbin+'pce_quad -o '+str(input_pcorder)+' -w HG -f '+filename+ ' > pcq.log; mv PCcoeff.dat param_pcf.txt'
	os.system(cmd)


else:
	print("prepare_inpc.py : Input format not recognized. Must be marg or sam.")