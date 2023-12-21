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
/// \file mcmc.h
/// \author K. Sargsyan, C. Safta, B. Debusschere, 2012 -
/// \brief Header file for the Markov chain Monte Carlo class

#ifndef UQTKMMALA_H_SEEN
#define UQTKMMALA_H_SEEN

#include "dsfmt_add.h"
#include "mcmc.h"
#include "arrayio.h"
#include "arraytools.h"
#include "mala.h"

#include <iostream>
#include <string.h>
#include <stdio.h>
#include <sstream>

//*****************************************

/// \class MMALA or Manifold variant of Langevian Sampling
/// \brief MMALA Markov Chain Monte Carlo class. Derived from the base class for MALA
///        Implemented the MMALA algorithms
class MMALA:public MALA{
public:
  //Delegating constructors
  MMALA(double (*logposterior)(Array1D<double>&, void *), void *postinfo):MALA(logposterior,postinfo){};
  MMALA(LogPosteriorBase& L):MALA(L){};
  /// \brief Get metric tensor function
  void getMetricTensor(void (*metricTensor)(Array1D<double>&, Array2D<double>&, void *));

  /// \brief Set the metric tensor function given a pointer to a metric tensor function
  void setMetricTensor(void (*metricTensor)(Array1D<double>&, Array2D<double>&, void *));


private:
  int nSubSteps_ = 1;
  // Proposal function
  void proposal(Array1D<double>& m_t,Array1D<double>& m_cand);

  double probOldNew(Array1D<double>& a, Array1D<double>& b) override;

  void (*metricTensor_)(Array1D<double>&, Array2D<double>&, void *); // Pointer to metric tensorr function

  void (*gradlogPosterior_)(Array1D<double>&, Array1D<double>&, void *);

  bool tensflag_ = false; // Flag that indicates whether tensor information is given or not

};

#endif
