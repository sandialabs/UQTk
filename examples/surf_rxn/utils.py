#=====================================================================================
#                     The UQ Toolkit (UQTk) version 3.0.4
#                     Copyright (2017) Sandia Corporation
#                     http://www.sandia.gov/UQToolkit/
#
#     Copyright (2013) Sandia Corporation. Under the terms of Contract DE-AC04-94AL85000
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
import numpy as npy

def get_npc(dims,order):
    if dims < 1 or order < 0:
        print "get_npc() unexpected dims/order values: ",dims,order,", -> Abort !"
        quit()
    nterms=1;
    if order == 0: return nterms;
    for i in range(order):
        nterms = nterms*(dims+order-i)
    for i in range(order):
        nterms = nterms/(order-i)
    return nterms

def compute_err(type,arr1,arr2):
    if ( len(arr1) != len(arr2) ):
        print "compute_err(): lengths of arr1,arr2 do not match:",len(arr1),len(arr2)
        quit()
    if type == "Linf":
        return abs(arr1-arr2).max()
    elif type == "Linfrel":
        return (abs(arr1-arr2).max())/abs(arr1).max()
    elif type == "L1":
        return npy.average(abs(arr1-arr2))
    elif type == "L1rel":
        return npy.average(abs(arr1-arr2))/abs(arr1).max()
    elif type == "L2":
        return npy.sqrt(npy.dot(arr1-arr2,arr1-arr2)/len(arr1))
    elif type == "L2rel":
        return npy.sqrt(npy.dot(arr1-arr2,arr1-arr2)/npy.dot(arr1,arr1))
    else:
        print "compute_err(): unknown method (use Linf, L1, L2, Linfrel, L1rel, L2rel) -> exit"
        quit()

