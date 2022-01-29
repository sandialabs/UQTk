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

# Directory forward and inverse UQ scripts
UQPC=${UQTK_INS}/examples/uqpc
IUQ=${UQTK_INS}/examples/iuq


## Generate posterior samples of model parameters
## (i.e. mix posterior with model error dimension, if any)
# Two embedded dimensions -r 2
# Second order embedded NISP -o 2
$IUQ/get_postsam.py -r 2 -o 2 -d 3 -t LU -x 100000 -s 1

## Scale to physical range
$UQPC/scale.x post_phparam.dat to prange.dat post_phparam_phys.dat

## Plot posterior PDFs of inputs
$IUQ/plot_pdfs.py -p post_phparam_phys.dat -n pnames.txt -g prange.dat -t tri
cp pdf_tri.eps pdf_tri_inputs.eps

## Prior predictive samples of surrogate
$IUQ/get_postpred.py -p qtrain.dat -f ytrain_surr.dat -n 70

## Posterior predictive samples of surrogate
$IUQ/get_postpred.py -p post_phparam.dat -f post_pred.dat -n 70

## Plot posterior PDFs of select outputs
# example of 0th, 10th and 33th outputs
$IUQ/plot_pdfs.py 0 10 33 -p post_pred.dat -n outnames.txt -t tri

## Plot prior and posterior together with data
$IUQ/plot_prpost.py 0 1 2 -y ytrain_surr.dat -z post_pred.dat -d ydata_all.txt -n outnames.txt

# $UQPC/transpose_file.x mapparam.dat > mapparamt.dat
# $UQPC/scale.x mapparamt.dat to prange.dat mapparamt_phys.dat
# $UQPC/transpose_file.x mapparamt_phys.dat > mapparam_phys.dat

## Select indices for plotting 1d slices, e.g. pick the second column =-0.11111
paste xdata_all.txt ydata_all.txt > xydata_all.txt
awk '($2+0.111)**2<0.001{print NR-1}' xdata_all.txt > ind_plot.dat
awk '($2+0.111)**2<0.001{print $3}' xydata_all.txt > ydata_plot.dat

## Sensitivities of selected
$IUQ/plot_sens.py -i ind_plot.dat -s allsens_total.dat -p pnames.txt -o outnames.txt

## Reorganize means and variances
awk '{print $1}' fmeans.dat > means.dat
paste fvars.dat surr_errors.dat | awk '{print $4*$4, $2, $1}' > vars.dat # SE, PU, ME

## Plot variance-decomposed fits
## -c 0 means plotting with respect to 0th column of xdata_all.txt
$IUQ/plot_fit1d.py -i ind_plot.dat -x xdata_all.txt -y ytrain.dat -d ydata_plot.dat -c 0 -r linear -m means.dat -v vars.dat

## Plot shaded quantiles of prior and posterior
## -c 0 means plotting with respect to 0th column of xdata_all.txt
$IUQ/plot_shade.py -i ind_plot.dat -x xdata_all.txt -y ytrain.dat -z post_pred.dat -d ydata_plot.dat -c 0
