#!/usr/bin/env python
#=====================================================================================
#                     The UQ Toolkit (UQTk) version 3.0.4
#                     Copyright (2017) Sandia Corporation
#                     http://www.sandia.gov/UQToolkit/
#
#     Copyright (2017) Sandia Corporation. Under the terms of Contract DE-AC04-94AL85000
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
	"Need numpy and matplotlib"

from PyUQTk.utils.func import *
from quad_tools import *
import matplotlib.pyplot as plt
import numpy as np
####################################################################


###### User Input #######

#Obtain desired model from user
model= raw_input("Please enter desired model from choices:\ngenz_osc\ngenz_exp\ngenz_cont\ngenz_gaus\
	\ngenz_cpeak\ngenz_ppeak\n\n")
#Check that model selected is one listed
model_choices=['genz_osc', 'genz_exp', 'genz_cont','genz_gaus','genz_cpeak', 'genz_ppeak']
while not(model in model_choices):
	model=raw_input("Please choose from listed models: ")

#Obtain number of dimensions from user
ndim=raw_input("Please enter desired dimension: ")
#Ask for valid input if dimension is not a positive integer
while not str.isdigit(ndim) or ndim=="0":
	ndim=raw_input("Please enter a positive integer: ")
#Change type of ndim to integer
ndim=int(ndim)

#Obtain desired maximum number of quadrature points from user
q_pts= raw_input("Enter the desired maximum number of quadrature points per dimension: ")
#Ask for valid input if q_pts is not a positive integer
while not str.isdigit(q_pts) or q_pts=="0":
	q_pts=raw_input("Please enter a positive integer: ")
#Change type of q_pts to integer
q_pts=int(q_pts)


#Let number of quad points per dimension vary from 1 to maximum number provided by user
num_points=range(1,q_pts+1)

#Calculate total number of quadrature points used for each value in num_points
tot_pts=[i**ndim for i in num_points]
#Let the number of sampling points for MC integrations be the same as the total number of quad points
mc_pts=tot_pts

###### Quadrature Integration ######

#list to store quad errors
q_errors=[]


#Loop though different values for number of quad points per dimension
for quad_param in num_points:
	#Generate quadrature points
	xpts,wghts=generate_qw(ndim,quad_param)

	#Evaluate the function
	func_params=ones(ndim+1) #Genz parameters
	ypts=func(xpts,model,func_params)

	#Quadrature integration
	integ=dot(ypts,wghts)
	#Evaluate its exact integral
	integ_ex=integ_exact(model,func_params)
	#Calculate error in quad integration and add to list
	q_error=abs(integ-integ_ex)
	q_errors.append(q_error)


###### Monte Carlo Integration ######

#Empty list to store MC errors
mc_errors=[]
#Let number of sampling points vary
for pts in mc_pts:
	#Calculate average error in 10 MC integrations and add to list
	error=find_error(pts,ndim,model,integ_ex, func_params)
	mc_errors.append(error)

###### Create Graph ######
#Create figure
fig, ax = plt.subplots(figsize=(10,10))
#Plot Quadrature Data
plt.plot(tot_pts, q_errors, color='r', label='Full Quadrature')
#Plot Monte Carlo Data
plt.plot(mc_pts, mc_errors, color='b', label='Monte Carlo')
#Label Axes
plt.xlabel("Total number of Sampling Points",fontsize=20)
plt.ylabel("Absolute Error in Integral Approximation",fontsize=20)
#Create legend
plt.legend(loc='lower left')

#Model titles to be displayed
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
fig.suptitle("Comparison of Full Quadrature and Monte Carlo Integration\nfor %s Model with\
 Dimension %s"%(model_title, str(ndim)), fontsize=22)
#Make scales logarithmic
plt.yscale('log')
plt.xscale('log')
#Change size of tick labels
plt.tick_params(axis='both', labelsize=16)

#Save figure
fig_file_name="quad_vs_mc.pdf"
plt.savefig(fig_file_name)
print "\nquad_vs_mc.pdf has been saved"
#Show figure
plt.show()
