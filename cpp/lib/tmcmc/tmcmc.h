/* =====================================================================================

                      The UQ Toolkit (UQTk) version 3.1.2
                          Copyright (2022) NTESS
                        https://www.sandia.gov/UQToolkit/
                        https://github.com/sandialabs/UQTk

     Copyright 2022 National Technology & Engineering Solutions of Sandia, LLC (NTESS).
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
/// \file mcmc.h
/// \author K. Sargsyan, C. Safta, B. Debusschere, 2012 -
/// \brief Header file for the Transitional Markov chain Monte Carlo class

#ifndef UQTKTMCMCClass_H_SEEN
#define UQTKTMCMCClass_H_SEEN

#include "dsfmt_add.h"
#include "mcmc.h"
#include "arrayio.h"
#include "arraytools.h"

#include <iostream>
#include <string.h>
#include <stdio.h>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <math.h>
#include <vector>
#include <algorithm>
#include <random>
#include <numeric>
#include <functional>
#include <string>

//*****************************************
/// \class Transitional MCMC
/// \brief Transitional Markov Chain Monte Carlo class. Derived from the base class for MCMC
///        Implemented the algorithms for TMCMC
class TMCMC:public MCMC{
public:

    // Delegating constructor
    TMCMC():MCMC(){};

    // Initalization and set functions:

    /// \brief Set the default values for TMCMC function
    void initDefaults();
    /// \brief Initialize the number of processes for TMCMC
    void initTMCMCNprocs(int tmcmc_nprocs);
    /// \brief Initialize the coefficient behind the covariance scaling
    ///        factor for TMCMC
    void initTMCMCGamma(double tmcmc_gamma);
    /// \brief Initialize the maximum allowed coefficient of variation for
    ///        the weights in TMCMC
    void initTMCMCCv(double tmcmc_cv);
    /// \brief Initialize the multiplicative factor for chain length to
    ///        encourage mixing in TMCMC
    void initTMCMCMFactor(int tmcmc_MFactor);
    /// \brief Initialize the choice to resample according to BASIS and
    ///        CATMIPs in TMCMC
    void initTMCMCBasis(bool tmcmc_basis);
    /// \brief Initialize the CATMIPs resampling parameter for TMCMC
    void initTMCMCCATSteps(int tmcmc_CATSteps);

    // Get functions:

    /// \brief Get the number of processes for TMCMC
    int getTMCMCNprocs();
    /// \brief Get the coefficient behind the covariance scaling
    ///        factor for TMCMC
    double getTMCMCGamma();
    /// \brief Get the maximum allowed coefficient of variation for
    ///        the weights in TMCMC
    double getTMCMCCv();
    /// \brief Get the multiplicative factor for chain length to
    ///        encourage mixing in TMCMC
    int getTMCMCMFactor();
    /// \brief Get the choice to resample according to BASIS and
    ///        CATMIPs in TMCMC
    bool getTMCMCBasis();
    /// \brief Get the CATMIPs resampling parameter for TMCMC
    int getTMCMCCATSteps();

    virtual void runChain(int ncalls, Array1D<double>& chstart) override;
    virtual void runChain(int ncalls) override;

private:
    int TMCMCNprocs; // The number of processes to use for parallel likelihood evaluation
    double TMCMCGamma; // The coefficient behind the covariance scaling factor
    double TMCMCCv; // The maximum allowed coefficient of variation for the weights
    int TMCMCMFactor; // The multiplicative factor for chain length to encourage mixing
    bool TMCMCBasis; // Resample according to BASIS and CATMIPs
    int TMCMCCATSteps; // CATMIPs resampling parameter
    std::vector<std::vector<double>> tmcmc_rngs; // The ranges for all samples

    // Initalization Flags
    bool tmcmcNprocsInit_ = false;
    bool tmcmcGammaInit_ = false;
    bool tmcmcCvInit_ = false;
    bool tmcmcMFactorInit_ = false;
    bool tmcmcBasisInit_ = false;
    bool tmcmcCATStepsInit_ = false;
    bool tmcmcRngsInit_ = false;

    // TMCMC Defaults
    int default_tmcmc_nprocs_ = 4;
    double default_tmcmc_gamma_ = -1.0;
    double default_tmcmc_cv_ = 0.1;
    int default_tmcmc_MFactor_ = 1;
    bool default_tmcmc_basis_ = true;
    int default_tmcmc_CATSteps_ = 1;
    ///\todo Write in the default TMCMC Rngs Vector
    std::vector<std::vector<double>> default_tmcmc_rngs_;

    double probOldNew(Array1D<double>& a, Array1D<double>& b) override {return 0.0;};
};

typedef std::vector<double>  RealVector;
typedef std::vector<int>     IntVector;
typedef std::vector<bool>    BoolVector;
typedef std::vector<char>    CharVector;

double tmcmc (RealVector &spls, RealVector &lprior, RealVector &llik,
           double gm, int nspl, int iseed,
           int nProcs, int ndim, double cv, int MFactor,
           bool basis, int CATSteps, int write_flag) ;

#endif
