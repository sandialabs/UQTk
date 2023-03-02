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
import numpy as npy

def get_npc(dims,order):
    if dims < 1 or order < 0:
        print("get_npc() unexpected dims/order values: %d, %d -> Abort !" % (dims,order))
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
        print("compute_err(): lengths of arr1,arr2 do not match: %d, %d" % (len(arr1),len(arr2)))
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
        print("compute_err(): unknown method (use Linf, L1, L2, Linfrel, L1rel, L2rel) -> exit")
        quit()

