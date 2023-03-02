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

