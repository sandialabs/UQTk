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
import sys
import os
uqtksrc=os.environ['UQTK_SRC']
sys.path.append(uqtksrc+"/../UQTk-install")

try:
	from numpy import *
	from matplotlib.pyplot import *
except ImportError:
	print("Need numpy and matplotlib")

from PyUQTk.utils.func import *
from quad_tools import *
import matplotlib.pyplot as plt
import numpy as np
########################################################

###### User Input ######
#Obtain desired model from user
model= input("Please enter desired model from choices:\ngenz_osc\ngenz_exp\ngenz_cont\ngenz_gaus\
	\ngenz_cpeak\ngenz_ppeak\n\n")
#Check that model selected is one listed
model_choices=['genz_osc', 'genz_exp', 'genz_cont','genz_gaus','genz_cpeak', 'genz_ppeak']
while not(model in model_choices):
	model=input("Please choose from listed models: ")

#Obtain number of dimensions from user
ndim=input("Please enter desired dimension: ")
#Ask for valid input if dimension is not a positive integer
while not str.isdigit(ndim) or ndim=="0":
	ndim=input("Please enter a positive integer: ")
#Change type of ndim to integer
ndim=int(ndim)

#Obtain level from user
level=input("Enter the maximum desired level: ")
#Ask for valid input if level is not a positive integer,
while not str.isdigit(level) or level=='0':
	level=input("Please enter a positive integer: ")
#Change type of level to integer
level=int(level)
#Create range of levels to loop through
levels=range(1,level+1)

#Set Genz parameters
if model=="genz_cpeak":
	func_params=ones(ndim+1)*0.1 #All Genz parameters set to 0.1
else:
	func_params=ones(ndim+1) #All Genz parameters set to 1

####### LU sparse Quadrature Integration ######
#List to store total number of quadrature points used
num_pts=[]
#List to store errors
sparse_errors=[]
#Do quadrature integration at varying levels, beginning at level 1 and going up to max level provided
for i in levels:
	#Generate quadrature points
	xpts,wghts=generate_qw(ndim,i, sp='sparse')

	#Evaluate the function
	ypts=func(xpts,model,func_params)
	#Determine number of quadrature points used and add to list
	pts=int(len(ypts))
	num_pts.append(pts)

	#Quadrature integration
	integ=dot(ypts,wghts)
	#Evaluate its exact integral
	integ_ex=integ_exact(model,func_params)
	#Calculate error and add to list of errors
	error=abs(integ_ex-integ)
	sparse_errors.append(error)

###### Monte Carlo Integration ######

#List to store MC errors
mc_errors=[]
#Let number of sampling points vary
#Will take on values of the total number of quad points used in each level
for pts in num_pts:
	#Calculate average error in 10 MC integrations and add to list
	error=find_error(pts,ndim,model,integ_ex, func_params)
	mc_errors.append(error)


###### Create Graph ######
#Create figure
fig, ax = plt.subplots(figsize=(10,10))
#Plot Quadrature Data
plt.plot(num_pts, sparse_errors, color='r', label='Sparse Quadrature')
#Plot MC Data
plt.plot(num_pts, mc_errors, color='b', label= 'Monte Carlo')

#Label Axes
plt.xlabel("Total Number of Sampling points",fontsize=20)
plt.ylabel("Absolute Error in Integral Approximation",fontsize=20)

#Model titles to be displayed on graph
if model=="genz_osc":
    model_title="Genz Oscillatory"
elif model=="genz_gaus":
    model_title="Genz Gaussian"
elif model=="genz_ppeak":
    model_title="Genz Product-Peak"
elif model=="genz_cpeak":
    model_title="Genz Corner-Peak"
elif model=="genz_exp":
    model_title="Genz Exponential"
else:
    model_title="Genz Continuous"

#Add title
fig.suptitle("Comparison of Sparse Quadrature and Monte Carlo Integration\nfor %s Model with\
 Dimension %s"%(model_title, str(ndim)), fontsize=22)
#Make scales logarithmic
plt.yscale('log')
plt.xscale('log')

#Change size of tick labels
plt.tick_params(axis='both', labelsize=16)
#Create legend
plt.legend(loc='lower left')

#Save figure
fig_file_name="sparse_quad.pdf"
plt.savefig(fig_file_name)
print("\nsparse_quad.pdf has been saved")
#Show figure
plt.show()
