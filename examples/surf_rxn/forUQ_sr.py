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


import os
import shutil
import sys
import numpy as np
import math
import random as rnd
from scipy.stats.mstats import mquantiles

import fileinput
import file_utils

def usage():
    print """
    forUQ_sr.py <pctype> <pcord> <method1> [<method2>] [<method3>]
    
    where: pctype: Type of PC: LU, HG, ...
           pcord : order of PCEs used for output variables
           method? : methods to pursue: ISP, NISP, and/or NISP_MC"""


# Input settings

if len(sys.argv) < 4:
    usage()
    exit(1)

pctype  = str(sys.argv[1])   # PC type
pcord   = int(sys.argv[2])   # PC order of the chain samples 
methods = sys.argv[3:]

nsam = 100                   # Number of samples for NISP_MC


nmeth=len(methods)
print methods
methods_str=' '.join(methods)

## Prepare the xml file
shutil.copyfile('forUQ_surf_rxn.in.xml.templ','surf_rxn.in.xml')
for line in fileinput.input("surf_rxn.in.xml", inplace = 1):
    print line.replace('PCTYPE', pctype),
for line in fileinput.input("surf_rxn.in.xml", inplace = 1):
    print line.replace('PCORDER', str(pcord)),

for im in range(nmeth):
    str_param=''
    if methods[im]=='NISP_MC':
        str_param=str(nsam)
    cmd = './SurfRxn'+methods[im]+'.x '+str_param
    print "Running",cmd
    os.system(cmd)
    cmd = 'mv solution.dat solution_'+methods[im]+'.dat'
    print "Running",cmd
    os.system(cmd)
  

# Postprocessing plots
    cmd = './plSurfRxnMstd.py '+methods[im]
    print "Running",cmd
    os.system(cmd)
 
cmd = './plPDF_method.py u ave '+pctype+' '+str(pcord)+' '+methods_str
print "Running",cmd
os.system(cmd)
cmd = './plPDF_method.py v ave '+pctype+' '+str(pcord)+' '+methods_str
print "Running",cmd
os.system(cmd)
cmd = './plPDF_method.py w ave '+pctype+' '+str(pcord)+' '+methods_str
print "Running",cmd
os.system(cmd)

cmd = './plPDF_method.py u tf '+pctype+' '+str(pcord)+' '+methods_str
print "Running",cmd
os.system(cmd)
cmd = './plPDF_method.py v tf '+pctype+' '+str(pcord)+' '+methods_str
print "Running",cmd
os.system(cmd)
cmd = './plPDF_method.py w tf '+pctype+' '+str(pcord)+' '+methods_str
print "Running",cmd
os.system(cmd)
