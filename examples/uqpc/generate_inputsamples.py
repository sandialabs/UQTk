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
# Example generation of 3d dependent input samples
################################################################################

# Import libraries
try:
	import numpy as np
except ImportError:
	print("Numpy module could not be found.")

# Number of samples
nsam=1000
# Container for samples
rv=np.empty((nsam,3))
# Sample generation
gausrv=np.random.normal(loc=1.0, scale=0.3, size=(nsam,3))
rv[:,0]=np.exp(gausrv[:,0])
rv[:,1]=gausrv[:,0]+np.log(gausrv[:,1]**2)
rv[:,2]=gausrv[:,2]*np.exp(-gausrv[:,1])

# Save to a file
np.savetxt('param_sam.txt',rv)

