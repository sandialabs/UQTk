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
/// \file amcmc.h
/// \author K. Sargsyan, C. Safta, B. Debusschere, 2012 -
/// \brief Header file for the Adaptive Markov Chain Monte Carlo class

#ifndef UQTKAMCMC_H_SEEN
#define UQTKAMCMC_H_SEEN

#include "dsfmt_add.h"
#include "mcmc.h"
#include "arrayio.h"
#include "arraytools.h"

#include <iostream>
#include <string.h>
#include <stdio.h>
#include <sstream>

//*****************************************

/// \class Adaptive MCMC
/// \brief Adaptive Markov Chain Monte Carlo class. Derived from the base class for MCMC
///        Implemented the algorithms for aMCMC
class AMCMC:public MCMC{
public:
    ///\brief Delegating Constructor
    AMCMC(double (*logposterior)(Array1D<double>&, void *), void *postinfo):MCMC(logposterior,postinfo){};
    AMCMC(LogPosteriorBase& L):MCMC(L){};


    // Initialization and set functions for private variables that are necessary to the aMCMC algorithms

    /// \brief Initialize adaptivity step parameters for aMCMC
    void initAdaptSteps(int adaptstart,int adaptstep, int adaptend);
    /// \brief Initialize the scaling factor gamma for aMCMC
    void initAMGamma(double gamma_);
    /// \brief Initialize the covariance 'nugget' for aMCMC
    void initEpsCov(double eps_cov_);

    // Get functions for any variables

    /// \brief Get function for the size=3 vector (t_start,t_step,t_end) that indicates when the adaptivity starts, how often the proposal covariance is updated and when the adaptivity ends, respectively
    void getAdaptSteps(Array1D<int> adaptstep_);
    /// \brief Get function for the coefficient behind the covariance scaling factor
    double getGamma();
    /// \brief Get function for the offset epsilon for Cholesky to be computationally feasible
    double getEpsCov();

    // Print functions:

    /// \brief Function to print the chain information on the screen
    void printChainSetup();

    virtual void runChain(int ncalls, Array1D<double>& chstart) override;
    virtual void runChain(int ncalls) override;

private:
    double gamma; // The coefficient behind the covariance scaling factor
    double eps_cov; // The offset epsilon for Cholesky to be computationally feasible
    int nSubSteps_ = 1;
    Array2D<double> curcov; // Covariance of the chain values so far
    Array1D<double> curmean; // Mean of the chain values sampled so far
    Array1D<int> adaptstep_; // a size=3 vector (t_start,t_step,t_end) that indicates when the adaptivity starts, how often the proposal covariance is updated and when the adaptivity ends, respectively

    // Flags to indicate if corresponding values are initated or not
    bool adaptstepInit_ = false;
    bool gammaInit_ = false;
    bool epscovInit_ = false;

    //Default values for Gamma and EPS_Cov
    double default_gamma_ = 0.01;
    double default_eps_cov_ = 1e-8;

    // Proposal Function
    void proposal(Array1D<double>& m_t,Array1D<double>& m_cand,int t);

    Array2D<double> propLCov_; // The Cholesky factor(square-root) of proposal covariance

    double probOldNew(Array1D<double>& a, Array1D<double>& b) override {return 0.0;};
};

#endif
