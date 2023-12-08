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

     Questions? Contact the UQTk Developers at https://github.com/sandialabs/UQTk/discussions
     Sandia National Laboratories, Livermore, CA, USA
===================================================================================== */
/// \file mcmc.h
/// \author K. Sargsyan, C. Safta, B. Debusschere, 2012 -
/// \brief Header file for the Markov chain Monte Carlo class

#ifndef UQTKMALA_H_SEEN
#define UQTKMALA_H_SEEN

#include "dsfmt_add.h"
#include "mcmc.h"
#include "arrayio.h"
#include "arraytools.h"

#include <iostream>
#include <string.h>
#include <stdio.h>
#include <sstream>

//*****************************************

/// \class MALA or Langevin Sampling
/// \brief MALA Markov Chain Monte Carlo class. Derived from the MCMC base class
class MALA:public MCMC{
public:
    // Delegating constructor
    MALA(double (*logposterior)(Array1D<double>&, void *), void *postinfo):MCMC(logposterior,postinfo){};
    MALA(LogPosteriorBase& L):MCMC(L){};

    // Initialization and set functions for private variables that are necessary to the MALA algorithms

    /// \brief Initialize epsilon for MALA
    void initEpsMALA(double eps_mala_);
    /// \brief Set the gradient function given a pointer to a gradient of a logPosterior function
    void setGradient(void (*gradlogPosterior)(Array1D<double>&, Array1D<double>&, void *));

    // Get functions:
    /// \brief Get epsilon for MALA
    double getEpsMALA();
    /// \brief Get gradient function by passing in log posterior pointer
    void getGradient(void (*gradlogPosterior)(Array1D<double>&, Array1D<double>&, void *));
    /// \brief Gets if the gradient function is set as a bool 
    bool getGradientFlag();

    /// \brief Function to evaluate the gradient of log-posterior
    void evalGradLogPosterior(Array1D<double>& m, Array1D<double>& grads);

    virtual void runChain(int ncalls, Array1D<double>& chstart) override;
    virtual void runChain(int ncalls) override;
    virtual void runOptim(Array1D<double>& start) override;

private:
    double eps_mala; // Epsilon for MALA algorithm
    int nSubSteps_ = 1;

    ///\brief Proposal Function
    void proposal(Array1D<double>& m_t,Array1D<double>& m_cand);

    // Flags to see if corresponding values are initalized or not
    bool epsMalaInit_ = false;

    double default_eps_mala_ = 0.1; // Default epsilon value for mala

    virtual double probOldNew(Array1D<double>& a, Array1D<double>& b) override;

    double evallogMVN_diag(Array1D<double>& x,Array1D<double>& mu,Array1D<double>& sig2); // Evaluate MVN

    void (*gradlogPosterior_)(Array1D<double>&, Array1D<double>&, void *);

    bool gradflag_ = false; // Flag to check whether gradient information is given or not
};

#endif
