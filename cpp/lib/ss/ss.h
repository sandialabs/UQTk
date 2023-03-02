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
/// \file mcmc.h
/// \author K. Sargsyan, C. Safta, B. Debusschere, 2012 -
/// \brief Header file for the Single Site Markov Chain Monte Carlo class

#ifndef UQTKSS_H_SEEN
#define UQTKSS_H_SEEN

#include "dsfmt_add.h"
#include "mcmc.h"
#include "arrayio.h"
#include "arraytools.h"

#include <iostream>
#include <string.h>
#include <stdio.h>
#include <sstream>

//*****************************************

/// \class Single-Site MCMC
/// \brief Single-Site Markov Chain Monte Carlo class. Derived from the base class for MCMC
///        Implemented the algorithms for single-site (Metropolis-within-Gibbs)
class SS:public MCMC{
public:
    // Delegating Constructors
    SS(double (*logposterior)(Array1D<double>&, void *), void *postinfo):MCMC(logposterior,postinfo){};
    SS(LogPosteriorBase& L):MCMC(L){};

    virtual void runChain(int ncalls, Array1D<double>& chstart) override;
    virtual void runChain(int ncalls) override;
    void setNSubSteps(){
      nSubSteps_ = this -> GetChainDim();
    };
    int getNSubSteps() override;
private:
    int nSubSteps_;

    // Proposal Function
    void proposal(Array1D<double>& m_t,Array1D<double>& m_cand,int dim);

    double probOldNew(Array1D<double>& a, Array1D<double>& b) override {return 0.0;};
};

#endif
