#!/bin/bash -e
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
#=====================================================================================

#######################################################################
#######################################################################
# Workflow for surrogate-based inference: use this as a readme of steps or run as an example
# Requirement: define UQTK_INS env. variable to point to installed location
# Make an empty folder and run this script there (there will be many files generated!)
# The script calls various other scripts for each step:
#   prep_model*.x and prep_data*.x in principle should be made problem specific
# Plotting utilities are mostly for quick demo: one can copy them and
#   rework for specific needs
#######################################################################
#######################################################################

# Directory forward and inverse UQ scripts
UQPC=${UQTK_INS}/examples/uqpc
IUQ=${UQTK_INS}/examples/iuq

#######################################################################
#######################################################################

# Step 1a: Prepare ensemble of model simulations
echo "Preparing model ensemble"

# Creates the following files in the given formats
#### Mandatory:
# prange.dat   :   d x 2 Parameter physical ranges
# ptrain.dat   :   N x d Parameter samples in physical ranges for training
# qtrain.dat   :   N x d Training samples scaled to [-1,1]^d
# ytrain.dat   :   N x L Model evaluations at training samples, a total of f outputs
#### Optional:
# pval.dat     :   V x d Parameter samples in physical ranges for validation
# qval.dat     :   V x d Validation samples scaled to [-1,1]^d
# yval.dat     :   V x L Model evaluations at validation samples, a total of f outputs
# pnames.txt   :   d x 1 List of parameter names, useful for plotting utilities
# outnames.txt :   L x 1 List of output names, useful for plotting utilities

# 1d example
#$IUQ/prep_model1.x
# 2d example
$IUQ/prep_model2.x > prep_model.out

#######################################################################
#######################################################################

# Step 1b: Build surrogates for each output given ensemble of simulations
echo "Building surrogates"

# Creates results.pk file with all the necessary information on f surrogates
# Number of outputs (f)
NOUT=`awk 'NR==1{print NF}' ytrain.dat`
# Number of training samples (N)
NTRAIN=`awk 'END{print NR}' ptrain.dat`
# Number of validation samples, if any (V)
if [ -f "pval.dat" ]; then
    NVAL=`awk 'END{print NR}' pval.dat`
else
    NVAL=0
fi
echo "Number of output QoIs        : $NOUT "
echo "Number of training samples   : $NTRAIN "
echo "Number of validation samples : $NVAL "

# Forward-UQ, i.e. build surrogate for all outputs
ORDER=3
${UQPC}/uq_pc.py -r offline_post -p prange.dat -m lsq -s rand -n $NTRAIN -v $NVAL -t $ORDER > uqpc.out

#######################################################################
#######################################################################

# Step 1c: Sample plotting examples, relying on results.pk
echo "Plotting surrogate results"

# Ideally, users should use these for quick verification,
# and use them as a guideline to build their own plotting routines

# Plot data-versus-model for surrogate accuracy assessment
${UQPC}/plot.py dm training validation > dm.out

# Plot runId versus data and model for surrogate accuracy assessment
#${UQPC}/plot.py idm training
#${UQPC}/plot.py idm validation

# Plot sensitivities (multi-output bars)
${UQPC}/plot.py sens main > sensmain.out
${UQPC}/plot.py sens total > senstotal.out

# Plot total sensitivities (matrix plots for all outputs and most relevant inputs)
# ${UQPC}/plot.py sensmat total

# Plot main and joint sensitivities for all outputs (circular plots)
#${UQPC}/plot.py senscirc

# Plot high-d representation of multiindex
#${UQPC}/plot.py mindex

#######################################################################
#######################################################################

# Step 1d: Reads results.pk to generate surrogate multiindex and coefficient text files,
# as well as estimated surrogate error file - all to be used for inference
echo "Preparing surrogate coefficents' files"
PMTYPE=pcs
$UQPC/get_modelpc.py -p $PMTYPE -t validation > gemodelpc.out


#######################################################################
#######################################################################

# Step 2a: Prepare observational/experimental data for inference
echo "Preparing observational data"

# Creates the following files in the given formats
#### Mandatory:
# xdata_all.txt    : L x S    L design conditions.
#                             Frequently, S=1 and this is simply counting index of the data
# ydata_all.txt    : L x E    Data values, possible r replicas per x-condition
#                             One can 'fold' the replicas as another column in x-conditions and keep E=1
# dataerr_all.dat  : L x 1    Standard deviation measure of the data. Often using surrogate error as a placeholder
# ind_calib.dat    : L x 1    Masking array of 0's and 1's to indicate which data points are to be used for inference
#### Optional:
# xcond_names.txt  : S x 1    Names of columns of x-condition file, useful for plotting utilities

# 1d example
# $IUQ/prep_data1.x
# 2d example
$IUQ/prep_data2.x > prepdata.out


#######################################################################
#######################################################################

# Step 2b: Select observational/experimental data and run inference
echo "Running inference"

# Prepare data for calibration
$IUQ/prep_calib.x $PMTYPE $NOUT

# Run the calibration
$IUQ/run_infer.x $PMTYPE

#######################################################################
#######################################################################

# Step 2c: Sample plotting examples, relying on outputs of model inference
echo "Postprocessing inference results"

## Currently hardwired for 2d example and embedded case, one should take the routines inside this as an example
$IUQ/postp_infer.x > postp_infer.out

#######################################################################
#######################################################################

# Save surrogate info in a tarball to not overload the directory
tar -cvzf mipcf.tar mindexp* pccfp*

# Cleanup, but be careful with the wildcards
rm -rf rvar.dat designPar.dat regparams.dat lambdas.dat selected.dat
rm -rf errors.dat sigma2.dat Sig.dat
rm -rf mi.dat mindex.dat mindex_new.dat PCcoeff.dat pccf.dat pcf.dat coeff.dat
rm -rf xdata.dat ydata.dat xcheck.dat ycheck.dat ycheck_var.dat data gdata.dat wghts.dat target qdpts.dat
rm -rf *mindex*dat pccf*dat
rm -rf mainsens.dat totsens.dat jointsens.dat
rm -rf *log *tmp
