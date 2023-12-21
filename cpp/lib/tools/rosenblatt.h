/* =====================================================================================

                      The UQ Toolkit (UQTk) version 3.1.4
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

     Questions? Contact the UQTk Developers at https://github.com/sandialabs/UQTk/discussions
     Sandia National Laboratories, Livermore, CA, USA
===================================================================================== */
/// \file rosenblatt.h
/// \brief Header for tools related to Rosenblatt transformation
#ifndef ROSENBLATT_H
#define ROSENBLATT_H


/// \brief Generates a new sample by inverse Rosenblatt defined by a given sample set (xi) and a 'uniform' sample from the unit hypercube
/// \note Performs inverse Rosenblatt map \f$\xi=R^{-1}(u)\in\mathbf{R}^d\f$ given a single sample of \f$u\in[0,1]^d\f$
/// \note The command line utility app/pce_quad uses this map; see UQTk Manual entry for the app for technical description of the map
/// \param[in] unif    : 1-dimensional array of size \f$d\f$ corresponding to a sample \f$u\in[0,1]^d\f$
/// \param[in] xi      : 2-dimensional array of size \f$d\times M\f$ corresponding to samples that define the target distribution
/// \param[out] newXi  : 1-dimensional array of size \f$d\f$ corresponding to a new sample \f$\xi=R^{-1}(u)\in\mathbf{R}^d\f$
/// \param[in] sig     : 1-dimensional array of size \f$d\f$ for dimension-specific KDE bandwidths
void invRos(Array1D<double>& unif, Array2D<double>& xi, Array1D<double>& newXi, Array1D<double>& sig);

/// \brief This is a version of invRos() with the same bandwidth (bw) for all dimensions
void invRos(Array1D<double>& unif, Array2D<double>& xi, Array1D<double>& newXi, double bw);

/// \brief This is a version of invRos() with an automatic bandwidth selection based on a rule of thumb.
/// \note The rule of thumb is not always reliable. It is recommended to test various bandwidths.
void invRos(Array1D<double>& unif, Array2D<double>& xi, Array1D<double>& newXi);

/// \brief This is a version of invRos() with an automatic bandwidth selection based on a rule of thumb, and
///        operating on set of uniform samples (rather than one at a time)
/// \param[in] unif    : 2-dimensional array of size \f$N\times d\f$ corresponding to N samples \f$u\in[0,1]^d\f$ (uniform)
/// \param[in] xi      : 2-dimensional array of size \f$d\times M\f$ corresponding to samples that define the arbitrary target distribution
/// \param[out] newXi  : 2-dimensional array of size \f$N\times d\f$ corresponding to a set of N new samples \f$\xi=R^{-1}(u)\in\mathbf{R}^d\f$
/// \note The rule of thumb is not always reliable. It is recommended to test various bandwidths.
void invRos(Array2D<double>& unif, Array2D<double>& xi, Array2D<double>& newXi);

/// \brief Calculates 'rule of thumb' optimal KDE bandwidths for a multi-dimensional data
/// \note Employs Silverman's rule-of-thumb, with a homemade factor adjustment accounting for samples that are near boundaries
/// \note This rule-of-thumb is quite heuristic; use at your own risk
/// \param[in] data    : 2-dimensional array of size \f$N\times d\f$ corresponding to samples
/// \param[out] bdwth  : 1-dimensional array of size \f$d\f$ for dimension-specific KDE bandwidths
void get_opt_KDEbdwth(const Array2D<double>& data, Array1D<double>& bdwth);

/// \brief Given input samples (xi), generates uniform samples by Rosenblatt map whose PDF is defined by a given sample set (xi_data)
/// \note Performs Rosenblatt map \f$u=R(\xi)\f$ given multiple samples of \f$\xi\f$
/// \param[in] xi      : 2-dimensional array of size \f$N\times d\f$ for input samples \f$\xi\f$ of the Rosenblatt map
/// \param[in] xi_data : 2-dimensional array of size \f$M\times d\f$ corresponding to samples \f$\xi_d\f$ from PDF defining the Rosenblatt map
/// \param[out] unif   : 2-dimensional array of size \f$N\times d\f$ for output samples of the Rosenblatt map; will land in \f$[0,1]^d\f$ if the input \f$\xi\f$ is sampled from the PDF of \f$\xi_d\f$
/// \param[in] sig     : 1-dimensional array of size \f$d\f$ for dimension-specific KDE bandwidths
void Rosen(Array2D<double>& xi, Array2D<double>& xi_data, Array2D<double>& unif, Array1D<double>& sig);

/// \brief This is a version of Rosen() with the same bandwidth (bw) for all dimensions
void Rosen(Array2D<double>& xi, Array2D<double>& xi_data, Array2D<double>& unif, double bw);

/// \brief This is a version of Rosen() with an automatic bandwidth selection based on a rule of thumb.
/// \note The rule of thumb is not always reliable. It is recommended to test various bandwidths.
void Rosen(Array2D<double>& xi, Array2D<double>& xi_data, Array2D<double>& unif);


#endif // ROSENBLATT_H
