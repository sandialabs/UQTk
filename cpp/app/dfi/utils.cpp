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

// Include header files
#include "utils.h"
#include "model.h"

void print_status(Array1D<double> &parameters, Info &info)
{
    // Print parameter values
    for (int i = 0; i < info.nb_of_parameters; i++)
    {
        std::cout << "  => p" << std::left << std::setw(4) << i + 1 << ":"
                  << " " << std::setw(10) << parameters(i) << std::endl;
    }

    // Print comparison between pushforward posterior and measurements
    for (int d = 0; d < info.nb_of_data_sets; d++)
    {
        std::cout << "  => Data set " << d + 1 << "/" << info.nb_of_data_sets << ":" << std::endl;
        for (int n = 0; n < info.nb_of_measurement_stations(d); n++)
        {
            Array2D<double> parameters_2d(0, parameters.XSize());
            parameters_2d.insertRow(parameters, 0);
            std::cout << "       -> x" << std::left << std::setw(4) << n + 1 << ":"
                      << " measurement = " << std::setw(10) << info.data_sets(d)(n, 0)
                      << ", prediction = " << std::setw(10) << eval_surrogate_model(parameters_2d, d, n, info) << std::endl;
        }
    }
}