/* =====================================================================================

                      The UQ Toolkit (UQTk) version 3.1.0
                          Copyright (2020) NTESS
                        https://www.sandia.gov/UQToolkit/
                        https://github.com/sandialabs/UQTk

     Copyright 2020 National Technology & Engineering Solutions of Sandia, LLC (NTESS).
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
#include "Array1D.h"

/// \brief Write Array data to a file with name filename
/// 
/// \note opens and closes the file
void WriteToFile(Array1D<double>& data, char* filename);

/// \brief Write all PC modes to a file
///
/// Arguments:
///   \li const double tym: current time
///   \li const double u*, v*, w*: vector with PC coefficients for three solution components u, v, w
///   \li const int n: number of PC terms
///   \li FILE* f_dump: C file pointer to write to
/// \note Assumes the file to be open already
void WriteModesToFilePtr(const double tym, const double* u, const double* v, const double* w, const int n, FILE* f_dump);


/// \brief Write mean and std. dev. to a file
///
/// Arguments:
///   \li const double tym: current time
///   \li const double u0, v0, w0: means of three solution components u, v, w
///   \li const double u_std, v_std, w_std: std. dev. of three solution components u, v, w
///   \li FILE* f_dump: C file pointer to write to
void WriteMeanStdDevToFilePtr(const double tym, const double u0, const double v0, const double w0, 
                                                const double u_std, const double v_std, const double w_std, FILE* f_dump);

/// \brief Write mean and std. dev. to screen
///
/// Arguments:
///   \li const int step: current time step
///   \li const double tym: current time
///   \li const double u0, v0, w0: means of three solution components u, v, w
///   \li const double u_std, v_std, w_std: std. dev. of three solution components u, v, w
void WriteMeanStdDevToStdOut(const int step, const double tym, const double u0, const double v0, const double w0, 
                                                       const double u_std, const double v_std, const double w_std);
