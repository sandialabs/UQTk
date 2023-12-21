
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

// Include std header files
#include <math.h>
#include <iomanip>
#include <cmath>

// Include other header files
#include "likelihood.h"
#include "model.h"

double eval_log_likelihood(Array2D<double> &parameters, Info &info)
{
    // TODO: maybe use constexpr here?
    const double pi = 4 * std::atan(1);

    // Check if we should update the pushforward posterior variance estimate and/or fisher information
    // TODO: this could be faster if we use templates
    bool compute_pushforward_variance = false;
    bool compute_fisher_information = false;
    if (info.c1 > info.burnin)
    {
        if (((info.c1 - 1) % info.every) == 0)
        {
            info.c2 += 1;
            compute_pushforward_variance = info.compute_pushforward_variance;
            compute_fisher_information = info.compute_fisher_information;
        }
    }

    // Allocate array to store the gradient (for computing the Fisher information matrix)
    Array2D<double> grad_g_d(1, info.nb_of_parameters);

    // Compute log likelihood
    double log_likelihood = 0;
    for (int d = 0; d < info.nb_of_data_sets; d++)
    {
        double alpha_d = info.weights(d);
        double K_d = info.nb_of_synthetic_data_sets(d);
        for (int n = 0; n < info.nb_of_measurement_stations(d); n++)
        {
            // Compute log likelihood constribution
            double s_d = info.data_sets(d)(n, 1);
            double beta_s_d_squared = info.betas(d) * s_d * s_d;
            double g_d = eval_surrogate_model(parameters, d, n, info);
            double numerator = 0;
            double numerator_sq = 0;
            for (int k = 0; k < K_d; k++)
            {
                double delta = info.synthetic_data_sets(d)(n, k) - g_d;
                numerator += delta;
                numerator_sq += delta * delta;
            }
            log_likelihood += alpha_d * (std::log(2 * pi * beta_s_d_squared) + numerator_sq / (K_d * beta_s_d_squared));

            // Update pushforward standard deviations
            if (compute_pushforward_variance)
            {
                double d1 = g_d - info.m(d)(n);
                info.m(d)(n) += d1 / info.c2;
                double d2 = g_d - info.m(d)(n);
                info.s(d)(n) += d1 * d2;
            }

            // Update fisher information matrix
            // WARNING: gradient computation in UQTk is really slow at the moment
            // Could be faster if we manually obtain the derivative PCE
            if (compute_fisher_information)
            {
                eval_grad_surrogate_model(grad_g_d, parameters, d, n, info);
                for (int i = 0; i < info.nb_of_parameters; i++)
                {
                    double grad_log_likelihood = alpha_d * numerator * grad_g_d(0, i) / (K_d * beta_s_d_squared);
                    double d = grad_log_likelihood * grad_log_likelihood - info.grad_log_likelihood(i);
                    info.grad_log_likelihood(i) += d / info.c2;
                }
            }
        }
    }

    // Return log likelihood
    return -0.5 * log_likelihood;
}

double eval_log_prior(Array2D<double> &parameters, Info &info)
{
    // TODO: maybe use constexpr here?
    const double pi = 4 * std::atan(1);

    // Compute the log prior
    double log_prior = 0;
    for (int i = 0; i < info.nb_of_parameters; i++)
    {
        int prior_type = (int) info.prior(i)(0);
        if (prior_type == 0) // Uniform prior
        {
            double a = info.prior(i)(1);
            double b = info.prior(i)(2);
            if (parameters(0, i) < b && parameters(0, i) > a)
            {
                log_prior -= std::log(b - a);
            }
            else
            {
                return -1e80; // return some large value
            }
        }
        else if (prior_type == 1) // Gaussian prior
        {
            double mu = parameters(0, i) - info.prior(i)(1);
            double sigma = info.prior(i)(2);
            log_prior -= (mu * mu) / (2 * sigma * sigma) + 0.5 * log(2 * pi) + log(sigma);
        }
    }

    // Return log prior
    return log_prior;
}

double eval_log_posterior(Array1D<double> &parameters, void *myinfo)
{
    // Cast to Info object
    Info *info = (Info*) myinfo;

    // Update sample counter
    info->c1 += 1;

    // Parameters to 2D (needed for surrogate model evaluation)
    Array2D<double> parameters_2d(0, parameters.XSize());
    parameters_2d.insertRow(parameters, 0);

    // Return result
    return eval_log_likelihood(parameters_2d, *info) + eval_log_prior(parameters_2d, *info);
}