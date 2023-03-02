/* =====================================================================================

                      The UQ Toolkit (UQTk) version 3.1.3
                          Copyright (2023) NTESS
                        https://www.sandia.gov/UQToolkit/
                        https://github.com/sandialabs/UQTk

     Copyright 2023 National Technology & Engineering Solutions of Sandia, LLC (NTESS).
     Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government
     retains certain rights in this software.

     This file is part of The UQ Toolkit (UQTk)

     UQTk is open source software: you can redistribute it and/or modify
     it under the terms of BSD 3-Clause License

     UQTk is distributed in the hope that it will be useful,
     but WITHOUT ANY WARRANTY; without even the implied warranty of
     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
     BSD 3 Clause License for more details.

     You should have received a copy of the BSD 3 Clause License
     along with UQTk. If not, see https://choosealicense.com/licenses/bsd-3-clause/.

     Questions? Contact the UQTk Developers at <uqtk-developers@software.sandia.gov>
     Sandia National Laboratories, Livermore, CA, USA
===================================================================================== */

#ifndef LIKELIHOOD_H
#define LIKELIHOOD_H

// Include UQTk header files
#include "arrayio.h"
#include "utils.h"

/**
 * @brief Evaluate the log of the likelihood function at the given parameters. 
 * 
 * @param parameters The parameter values where to evaluate the likelihood.
 * @param info Struct that contains problem-specific parameters.
 * @return The log of the likelihood evaluated at the given parameter values.
 */
double eval_log_likelihood(Array2D<double> &parameters, Info &info);

/**
 * @brief Evaluate the log of the prior at the given parameters. 
 * 
 * @param parameters The parameter values where to evaluate the prior.
 * @param info Struct that contains problem-specific parameters.
 * @return The log of the prior evaluated at the given parameter values.
 */
double eval_log_prior(Array2D<double> &parameters, Info &info);

/**
 * @brief Evaluate the log of the posterior at the given parameters. 
 * 
 * @param params Parameters where the posterior should be evaluated.
 * @param info Struct that contains problem-specific parameters.
 * @return The log of the posterior evaluated at the given parameter values.
 */
double eval_log_posterior(Array1D<double> &params, void *info);

#endif