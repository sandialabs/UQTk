/* =====================================================================================

                      The UQ Toolkit (UQTk) version @UQTKVERSION@
                          Copyright (@UQTKYEAR@) NTESS
                        https://www.sandia.gov/UQToolkit/
                        https://github.com/sandialabs/UQTk

     Copyright @UQTKYEAR@ National Technology & Engineering Solutions of Sandia, LLC (NTESS).
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

#ifndef UTILS_H
#define UTILS_H

// Include std header files
#include <iomanip>
#include <fstream>

// Include UQTk header files
#include "arrayio.h"
#include "PCSet.h"

/**
 * @brief Struct to hold problem-specific parameters.
 */
struct Info
{
    // Number of uncertain parameters to calibrate
    int nb_of_parameters;
    // List of prior types and prior parameters for each parameter (length nb_of_parameters)
    Array1D<Array1D<double> > prior;
    // Number of (experimental) data sets
    int nb_of_data_sets;
    // Number of synthetic data sets for each (experimental) data set (length nb_of_data_sets)
    Array1D<int> nb_of_synthetic_data_sets;
    // Number of measurement stations for each (experimental) data set (length nb_of_data_sets)
    Array1D<int> nb_of_measurement_stations;
    // (Experimental) data sets (length nb_of_data_sets, size nb_of_measurement_stations x 2 each)
    Array1D< Array2D<double> > data_sets;
    // The PCE multi-indices (length nb_of_data_sets, length nb_of_measurement_stations, size nb_of_terms x nb_of_parameters each)
    Array1D< Array1D< Array2D<int> > > mindices;
    // The PCE coefficients for each data set (length nb_of_data_sets, length nb_of_measurement_stations, length nb_of_terms each)
    Array1D< Array1D< Array1D<double> > > pccfs;
    // List of (pointers to) PCE surrogate models (length nb_of_data_sets)
    Array1D< Array1D<PCSet* > > pces;
    // Data standard deviation scaling factors (length nb_of_data_sets)
    Array1D<double> betas;
    // Synthetic data sets (length nb_of_data_sets, size nb_of_measurement_stations x nb_of_synthetic_data_sets)
    Array1D< Array2D<double> > synthetic_data_sets;
    // Data set weights to be used in the likelihood function (length nb_of_data_sets)
    Array1D<double> weights;
    // Integer used to indicate the chain burnin factor
    int burnin;
    // Integer used to indicate the chain subsampling rate
    int every;
    // Variables used for incremental standard deviation updates
    bool compute_pushforward_variance;
    int c1;
    int c2;
    Array1D<Array1D<double> > m;
    Array1D<Array1D<double> > s;
    // Variables used for Fisher information matrix
    bool compute_fisher_information;
    Array1D<double> grad_log_likelihood;
};

/**
 * @brief Print current status of MCMC chain.
 *
 * @param parameters Point where to evaluate the surrogate model.
 * @param info Struct that contains problem-specific parameters.
 */
void print_status(Array1D<double> &parameters, Info &info);

/**
 * @brief Check if a file with the given name exists
 * 
 * @param name The file name to check
 * @return True if a file with the given name exists, false otherwise
 */
// Note: this can be replaced with std::filesystem::exists once UQTk switches to c++17
inline bool file_exists(const std::string& name) {
    ifstream f(name.c_str());
    return f.good();
}

#endif
