#!/usr/bin/env python
#=====================================================================================
#
#                      The UQ Toolkit (UQTk) version 3.1.2
#                          Copyright (2022) NTESS
#                        https://www.sandia.gov/UQToolkit/
#                        https://github.com/sandialabs/UQTk
#
#     Copyright 2022 National Technology & Engineering Solutions of Sandia, LLC (NTESS).
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
import scipy as spy
import os
import shutil
import optparse
import matplotlib.pyplot as plt
from scipy import stats, mgrid, c_, reshape, random, rot90

# local definitions
from utils          import get_npc, compute_err
from forUQ_BB_utils import funcBB, plotKdeBin, plotFunQdpts

# define parser options
parser = optparse.OptionParser()
parser.add_option("--nom", dest="nomvals", default=None,
                  help="Nominal parameter values");
parser.add_option("-s","--std", dest="stdfac", default=0.1, type="float",
                  help="Percentage factor: std/mean")
parser.add_option("-d","--dim", dest="dim", default=1, type="int",
                  help="no. of dimensions (work with first dim parameters)")
parser.add_option("-l","--lev", dest="lev", default=21, type="int",
                  help="Sparsity level or number of quad points per dimension")
parser.add_option("-o","--ord", dest="ord", default=20, type="int",
                  help="Output response surface order")
parser.add_option("-q","--qfs", dest="sp", default="full", type="string",
                  help="full/sparse quadrature")
parser.add_option("--npces", dest="npces", default=100000, type="int",
                  help="No. of output PCE samples for KDE evaluations")
parser.add_option("--npdf", dest="npdf", default=100, type="int",
                  help="No. of grids for KDE evaluations")

# Extract command line options
options, remainder = parser.parse_args()
stdfac = options.stdfac
dim    = options.dim
lev    = options.lev
rord   = options.ord
qfs    = options.sp
npces  = options.npces
npdf   = options.npdf

# TODO: have a usage() command that displays help info when -h flag is given

# command index
icmd=0

if options.nomvals is None:
    #nom=npy.array([1.6,20.75,0.04,1.0,0.36,0.016])
    nom=npy.array([20.75,20.75,0.04,1.0,0.36,0.016])
else:
    dsize=options.nomvals.split(",")
    if ( len(dsize) < dim ):
        print("Error: no. of parameters < no. of dimensions: ",len(dsize),dim)
        quit()
    nom=npy.array([float(x.strip()) for x in dsize])

# define uqtkbin
if os.environ.get("UQTK_INS") is None:
    print("Error: Need to set path to UQTk install direactory as environment variable UQTK_INS -> Abort")
    quit()
else:
    if ( not os.path.isdir(os.environ["UQTK_INS"]) ):
        print("\"",os.environ["UQTK_INS"],"\" is not a valid path -> Abort")
        quit()

print("-----------------------------------------------------------------------------------------------------")
print("   ___                   _         ")
print("  |_ _|_ __  _ __  _   _| |_   _   ")
print("   | || '_ \| '_ \| | | | __| (_)  ")
print("   | || | | | |_) | |_| | |_   _   ")
print("  |___|_| |_| .__/ \__,_|\__| (_)  ")
print("            |_|                    ")
print("-----------------------------------------------------------------------------------------------------")

#--------------------------------------------------------------------------------------------------
# Define utilities
#--------------------------------------------------------------------------------------------------
pcType="HG"
uqtkbin=os.environ["UQTK_INS"]+"/bin"
print("Path to uqtk binaries:")
print("          ",uqtkbin)

genqd  = uqtkbin+"/generate_quad"
pcerv  = uqtkbin+"/pce_rv"
pceval = uqtkbin+"/pce_eval"
pcrsp  = uqtkbin+"/pce_resp"

foundAllBins = True
for fl in [genqd,pcerv,pceval,pcrsp]:
    if ( not os.path.isfile(fl) ):
        print("Could not find:",fl)
        foundAllBins = False

if ( not foundAllBins ):
    print("Could not find all required utilities -> Abort")
    quit()

print("Utilities: ")
print("          ",genqd)
print("          ",pcerv)
print("          ",pceval)
print("          ",pcrsp)
print()

#--------------------------------------------------------------------------------------------------
# Echo command-line parameters
#--------------------------------------------------------------------------------------------------
print("No. of input parameters:..................",dim)
print("Nominal parameter values:.................",nom[:dim])
print("Percentage factor, std/mean:..............",stdfac)
print("Full/sparse quadrature:...................",qfs)
print("Sparsity level or no. of quad points/dim:.",lev)
print("Output response surface order:............",rord)
print("No. of PCE samples for KDE evaluations:...",npces)
print("No. of grids for KDE evaluations:.........",npdf)

#--------------------------------------------------------------------------------------------------
# Prepare file with mean and std - first order Hermite-Gauss chaos
#--------------------------------------------------------------------------------------------------
nomdim=nom[:dim]
stdmat=npy.diag(stdfac*nomdim)
npy.savetxt("pcfile", nomdim, fmt="%.12e", delimiter=" ", newline="\n")
f_handle = open("pcfile", "a+")
npy.savetxt(f_handle, stdmat, fmt="%.12e", delimiter=" ", newline="\n")
f_handle.close()
print("Mean of model parameters:")
print(nomdim)
print("Std.dev. of model parameters:")
print(stdmat)

#-check number of PC terms
inord=1
pccf = npy.genfromtxt('pcfile')
if len(pccf.shape) == 1:
    pccf=npy.array(npy.transpose([pccf]))
indim_par = pccf.shape[1]
npc = pccf.shape[0]
npcc = get_npc(dim,inord)
if npc != npcc:
    print("forUQ_BB(): The number of input PC coefficients does not match to the given dimension and order -> Abort !\n")
    quit()

print("-----------------------------------------------------------------------------------------------------")
print("  ____                    _                 ")
print(" |  _ \ _   _ _ __  _ __ (_)_ __   __ _   _ ")
print(" | |_) | | | | '_ \| '_ \| | '_ \ / _` | (_)")
print(" |  _ <| |_| | | | | | | | | | | | (_| |  _ ")
print(" |_| \_\\\__,_|_| |_|_| |_|_|_| |_|\__, | (_)")
print("                                  |___/     ")
print("-----------------------------------------------------------------------------------------------------")
#--------------------------------------------------------------------------------------------------
# Construct a PCE expansion of model output as a function of input chaos germ
#--------------------------------------------------------------------------------------------------
# 1. generate quadrature points
qdcmd=genqd+" -d"+str(dim)+" -g'HG' -x"+qfs+" -p"+str(lev)+" > logQuad.dat ..."
print("Running:")
print("     ",qdcmd.replace(uqtkbin+"/",""), end=' ')
os.system(qdcmd)
print(" Done !")
shutil.copy("qdpts.dat","xdata.dat")

#--------------------------------------------------------------------------------------------------
# 2. evaluate input PCs at quadrature points
#--------------------------------------------------------------------------------------------------
open("input.dat", 'w').close()
open("logEvalInPC.dat", 'w').close()
#inpceval=pceval+" -x'PC' -p"+str(dim)+" -q"+str(inord)+" -f'pccf.dat' -s"+pcType+" >> logEvalInPC.dat"
inpceval=pceval+" -x'PC' -o"+str(inord)+" -f'pccf.dat' -s"+pcType+" >> logEvalInPC.dat"
print("     ",inpceval.replace(uqtkbin+"/",""),"->", end=' ')
for k in range(indim_par):
    npy.savetxt("pccf.dat",pccf[:,k],fmt='%.18e', delimiter=' ', newline='\n')
    print(k+1, end=' ')
    os.system(inpceval)
    os.system("paste input.dat ydata.dat > tmp; mv tmp input.dat")

print(" Done !")

#--------------------------------------------------------------------------------------------------
# 3. evaluate the forward function
#--------------------------------------------------------------------------------------------------
print("      black-box script ->", end=' ')
funcBB("input.dat","output.dat",xmltpl="surf_rxn.in.xml.tp3",xmlin="surf_rxn.in.xml")
shutil.copy("output.dat","ydata.dat")
print(" Done !")

# generate 1000 samples between min and max of input.dat and evaluate model
idat=npy.genfromtxt("input.dat")
ivals=npy.linspace(idat.min(),idat.max(),1000)
npy.savetxt("input_fcn.dat",ivals,fmt='%.18e', delimiter=' ', newline='\n')
funcBB("input_fcn.dat","output_fcn.dat",xmltpl="surf_rxn.in.xml.tp3",xmlin="surf_rxn.in.xml")
plotFunQdpts("input_fcn.dat","output_fcn.dat","input.dat","output.dat","forUQ_BB_model.pdf")

# 4. build a PC expansion for the output
prspcmd = pcrsp+" -x"+pcType+" -o"+str(rord)+" -d"+str(dim)+" -e > logPCEresp.dat"
print("     ",prspcmd.replace(uqtkbin+"/",""),"->", end=' ')
os.system(prspcmd)
os.system("paste mindex.dat PCcoeff_quad.dat > mipc.dat")
print(" Done !")

#--------------------------------------------------------------------------------------------------
# Compute L2 error between model output and PCE values
#--------------------------------------------------------------------------------------------------
# evaluate input PCs at random values drawn from its type
# No need to specify an order or an input file with PC coefficients as we are
# sampling from a standard random variable that is fully determined by pcType
pcrvcmd=pcerv+" -w'PCvar' -x"+pcType+" -d"+str(dim)+" -n200 -p"+str(dim)+"  > logPCrv.dat"
print("     ",pcrvcmd.replace(uqtkbin+"/",""),"->", end=' ')
os.system(pcrvcmd)
shutil.copy("rvar.dat","xdata.dat")
print(" Done !")

open("input_val.dat", 'w').close()
open("logEvalOutPCrnd.dat", 'w').close()
#outpceval=pceval+" -x'PC' -p"+str(dim)+" -q"+str(inord)+" -f'pccf.dat' -s"+pcType+" >> logEvalInPCrnd.dat"
outpceval=pceval+" -x'PC' -o"+str(inord)+" -f'pccf.dat' -s"+pcType+" >> logEvalInPCrnd.dat"
print("     ",outpceval.replace(uqtkbin+"/",""),"->", end=' ')
for k in range(indim_par):
    npy.savetxt("pccf.dat",pccf[:,k],fmt='%.18e', delimiter=' ', newline='\n')
    os.system(outpceval)
    os.system("paste input_val.dat ydata.dat > tmp; mv tmp input_val.dat")

os.system("rm -rf solution*.dat")
print(" Done !")

# evaluate the forward function and its PC surrogate
print("      black-box script ->", end=' ')
funcBB("input_val.dat","output_val.dat",xmltpl="surf_rxn.in.xml.tp3",xmlin="surf_rxn.in.xml")
print(" Done !")
#outpceval=pceval+" -x'PC' -p"+str(dim)+" -q"+str(rord)+" -f'PCcoeff_quad.dat' -s"+pcType+" > logEvalOutPCrnd.dat"
outpceval=pceval+" -x'PC' -o"+str(rord)+" -f'PCcoeff_quad.dat' -s"+pcType+" > logEvalOutPCrnd.dat"
print("     ",outpceval.replace(uqtkbin+"/",""),"->", end=' ')
os.system(outpceval)
shutil.move("ydata.dat","output_val_pc.dat")
print(" Done !")

#--------------------------------------------------------------------------------------------------
#  cleanup
#--------------------------------------------------------------------------------------------------
os.system("paste input.dat output.dat > xydata.dat")
os.system("paste input.dat ydata_pc.dat > xymodel.dat")
os.system("paste input_val.dat output_val.dat > xydata_val.dat")
os.system("paste input_val.dat output_val_pc.dat > xymodel_val.dat")
os.system("paste input.dat output.dat ydata_pc.dat > xyp.dat")

#--------------------------------------------------------------------------------------------------
#  Extract sensitivity information
#--------------------------------------------------------------------------------------------------
senscmd=uqtkbin+"/pce_sens -m'mindex.dat' -f'PCcoeff_quad.dat' -x'HG' > logSens.dat"
print("     ",senscmd.replace(uqtkbin+"/",""),"->", end=' ')
os.system(senscmd)
print(" Done !")

#--------------------------------------------------------------------------------------------------
#  Sample PCE and plot pdf
#--------------------------------------------------------------------------------------------------
pcerv=uqtkbin+"/pce_rv"
splcmd=pcerv+" -w'PC' -f'PCcoeff_quad.dat' -x'HG' -d1 -n"+str(npces)+" -p"+str(dim)+" -o"+str(rord)+" > logPCEspl.dat"
print("     ",splcmd.replace(uqtkbin+"/",""),"->", end=' ')
os.system(splcmd)
shutil.move("rvar.dat","PCEspls_"+str(npces/1000)+"k.dat")
print(" Done !")
print("      plot pdf")
plotKdeBin("PCEspls_"+str(npces/1000)+"k.dat","forUQ_BB_PCEdens_"+str(npces/1000)+"k.pdf",npdf=npdf)
print("-----------------------------------------------------------------------------------------------------")

# evaluate L2 error between direct function evaluation and the PC values
outval   =npy.genfromtxt("output_val.dat")
outval_pc=npy.genfromtxt("output_val_pc.dat")
print("Relative L2 error between model output and PCE surrogate:",compute_err("L2rel",outval_pc,outval));

# Cleanup
os.system("rm -rf indices.dat n.dat nex.dat norms.dat pccf.dat")
os.system("rm -rf qdpts* wghts* xdata* ydata* mindex.dat")
os.system("rm -rf pcfile slu.dat")
