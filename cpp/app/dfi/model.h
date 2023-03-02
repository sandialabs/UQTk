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

#ifndef MODEL_H
#define MODEL_H

// Include header files
#include "utils.h"

/**
 * @brief Evaluate the model. 
 * 
 * @param parameters Point where to evaluate the surrogate model.
 * @param d Data set number.
 * @param n Measurement station number.
 * @param info Struct that contains problem-specific parameters.
 * @return The value of the surrogate model evaluated at the given parameter values.
 */
double eval_surrogate_model(Array2D<double> &parameters, int d, int n, Info &info);

/**
 * @brief Evaluate the gradient of the model. 
 * 
 * @param grad array to store the gradient.
 * @param parameters Point where to evaluate the gradient of the surrogate model.
 * @param d Data set number.
 * @param n Measurement station number.
 * @param info Struct that contains problem-specific parameters.
 */
void eval_grad_surrogate_model(Array2D<double> &grad, Array2D<double> &parameters, int d, int n, Info &info);

#endif