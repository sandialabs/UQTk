/* =====================================================================================
                     The UQ Toolkit (UQTk) version @UQTKVERSION@
                     Copyright (@UQTKYEAR@) Sandia Corporation
                     http://www.sandia.gov/UQToolkit/

     Copyright (@UQTKYEAR@) Sandia Corporation. Under the terms of Contract DE-AC04-94AL85000
     with Sandia Corporation, the U.S. Government retains certain rights in this software.

     This file is part of The UQ Toolkit (UQTk)

     UQTk is free software: you can redistribute it and/or modify
     it under the terms of the GNU Lesser General Public License as published by
     the Free Software Foundation, either version 3 of the License, or
     (at your option) any later version.

     UQTk is distributed in the hope that it will be useful,
     but WITHOUT ANY WARRANTY; without even the implied warranty of
     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
     GNU Lesser General Public License for more details.

     You should have received a copy of the GNU Lesser General Public License
     along with UQTk.  If not, see <http://www.gnu.org/licenses/>.

     Questions? Contact Bert Debusschere <bjdebus@sandia.gov>
     Sandia National Laboratories, Livermore, CA, USA
===================================================================================== */
/// \file inference.h
/// \author K. Sargsyan  2016 -
/// \brief Header for the model inference tools

#include <string>

#include "Array1D.h"
#include "Array2D.h"


/// \brief Main function for inferring model parameters
/// \note This is written in a fortran style, i.e. some arguments are inputs, and the rest are output
/// \param[in] *forwardFuncs          : an array of y=f(p,x) functions that take np-by-pdim and nx-by-xdim input arrays
///                                     and returns an np-by-nx output, see tools/func.h for several examples
/// \param[in] funcinfo               : auxiliary function-specific information (can be 0)
/// \param[in] likType                : likelihood type: options are
///                                    'full', 'marg', 'mvn', 'gausmarg', 'abc', 'abcm', 'classical', 'koh' (see UQTk Manual)
/// \param[in] priorType              : prior type: options are 'uniform', 'normal', 'inverse', etc ... (see UQTk Manual)
/// \param[in] priora, priorb         : prior parameters (todo: need to be made dimension-specific)
///                                     for uniform, it is the range, for normal it is the moments
/// \param[in] xdata                  : x-values of data, nx-by-xdim
/// \param[in] ydata                  : y-values of data, nx-by-neach
/// \param[in] xgrid                  : x-values where predictive moments are computed after the inference
///                                   : can be the same as xdata
/// \param[in] dataNoiseInference     : indicator, data noise stdev is fixed(0), inferred (1), or log-inferred (2)
/// \param[in] datanoise_array        : data noise stdev array, if fixed, otherwise merely an MCMC starting point
/// \param[in] pdim                   : model parameter dimensionality
/// \param[in] order                  : order of output PC that is computed via NISP in the likelihood
/// \param[in] rndInd                 : array of indices of parameters to be randomized
/// \param[in] fixIndNom              : array of indices and nominal values of parameters to be fixed
/// \param[in] pdfType                : type of PDF PC parameterization, options are 'pct','pci' and 'full' (see UQTk Manual)
/// \param[in] pcType                 : type of PC for the PDF parameterization, options are all common PC types, e.g. 'HG','LU'
/// \param[in] seed                   : integer seed for MCMC
/// \param[in] nmcmc                  : number of MCMC steps to follow optimization; if 0, then only optimization is performed if optimflag is True
/// \param[in] mcmcgamma              : gamma (scaling) parameter for adaptive MCMC
/// \param[in] optimflag              : indicates if optimization is prepended to MCMC
/// \param[in] chstart                : initial chain state
/// \param[in] chsig                  : initial non-adaptive, dimensionwise proposal jump size
/// \param[in] likParam, likParam_int : likelihood parameters (currently, only needed for KDE-based likelihoods,
///                                    'full' and 'marg', to pass the KDE bandwidth and number of samples)
/// \param[in] pgrid                  : parameter grid, if requested, to compute exact posterior
///                                     (can be empty)
/// \param[in] nburn                  : burn-in for MCMC to write to pchain
/// \param[in] nstep                  : thinning of MCMC to write to pchain
/// \param[in,out] pchain             : thinned chain file with nburn and nstep applied, to be used for postprocessing
///                                   : if given as non-empty array, the MCMC is skipped, and only postprocessing is performed
/// \param[out] mapparam              : MAP parameters
/// \param[out] datavar_map           : MAP value of data variances
/// \param[out] pmean_map, pvar_map   : MAP values of parameter means and variances
/// \param[out] fmean_map, fvar_map   : MAP values of function predition (at xgrid) means and variances
/// \param[out] postave_datavar       : posterior average of data variances
/// \param[out] p_postave_mean        : posterior average of parameter mean
/// \param[out] p_postave_var         : posterior average of parameter variance
/// \param[out] p_postvar_mean        : posterior variance of parameter mean
/// \param[out] f_postave_mean        : posterior average of function prediction mean
/// \param[out] f_postave_var         : posterior average of function prediction variance
/// \param[out] f_postvar_mean        : posterior variance of function prediction mean
/// \param[out] paramPCcfs            : each column is a vector of parameter PC coefficients corresponding to an MCMC sample from pchain
///                                     the last column is the MAP value of parameter PC coefficients
void infer_model(Array1D< Array2D<double> (*)(Array2D<double>&, Array2D<double>&, Array2D<double>&, void *) > forwardFuncs, void* funcInfo,
	string likType,
	string priorType, double priora, double priorb,
	Array2D<double>& xdata,Array2D<double>& ydata, Array2D<double>& xgrid,
	int dataNoiseInference, Array1D<double>& datanoise_array,
	int pdim,int order,Array1D<int>& rndInd,Array2D<double>& fixIndNom,string pdfType,string pcType,
	int seed, int nmcmc, double mcmcgamma, bool optimflag, Array1D<double>& chstart, Array1D<double>& chsig,
	double likParam, int likParam_int,
	Array2D<double>& pgrid,Array2D<double>& pchain, int nburn, int nstep,
	Array1D<double>& mapparam, Array1D<double>& datavar_map,
	Array1D<double>& pmean_map, Array1D<double>& pvar_map,
	Array1D<double>& fmean_map, Array1D<double>& fvar_map,
	Array1D<double>& postave_datavar,
	Array1D<double>& p_postave_mean, Array1D<double>& p_postave_var, Array1D<double>& p_postvar_mean,
   	Array2D<double>& f_postsam_mean, Array1D<double>& f_postave_mean, Array1D<double>& f_postave_var, Array1D<double>& f_postvar_mean,
   	Array2D<double>& paramPCcfs);

/// \brief Log-posterior function given a vector of parameters (chain state) and a void* set of auxiliary variables
double LogPosterior(Array1D<double>& m, void* mypost_void);


