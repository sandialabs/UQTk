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

     Questions? Contact the UQTk Developers at <uqtk-developers@software.sandia.gov>
     Sandia National Laboratories, Livermore, CA, USA
===================================================================================== */
/// \file Sampling.hpp
/// \author B. Debusschere, C. Safta, K. Sargsyan, K. Chowdhary 2016 -
/// \brief Header file for the Sampling class

#ifndef SAMPLING_HPP_SEEN
#define SAMPLING_HPP_SEEN

#include <climits>
#include <math.h>

#include "gen_defs.h"
#include "probability.h"
#include "Array2D.h"
#include "dsfmt_add.h"


/// \class  Sampling
/// \brief  Provides functions for Monte-Carlo sampling, including
/// Latin Hypercube Sampling and Improved Hypercube Sampling; also
/// provides Quasi Monte Carlo methods

class Sampling {
public:

  /// \brief Constructor: initializes the Sampling class
  Sampling(const std::string sp_type, int ndim) {
    assert(ndim>0);
    if (sp_type==std::string("qmc")) {
      hBase_.resize(ndim,0.0);
      int ierr = getPrimes(ndim, hBase_);
      hSkip_.resize(ndim,0);
      hJump_.resize(ndim,1);
      hStep_ = 0;
    }
  }

  /// \brief Destructor: cleans up all memory and destroys object
  ~Sampling() {
    // do nothing
  }

  /// \brief Latin Hypercube Sampling, uniform distribution in all directions
  void unifLHS(const int nsample, const int ndim, const int zSeed, double *rvar);
  /// \brief Latin Hypercube Sampling, uniform distribution in all directions
  void unifLHS(const int nsample, const int ndim, dsfmt_t *rnstate, double *rvar);
  /// \brief Latin Hypercube Sampling, uniform distribution in all directions
  void unifLHS(const int zSeed, Array2D<double> &rvar);
  /// \brief Latin Hypercube Sampling, uniform distribution in all directions
  void unifLHS(dsfmt_t *rnstate, Array2D<double> &rvar);

  /// \brief Latin Hypercube Sampling, normal distribution in all directions
  void normLHS(const int nsample, const int ndim, const int zSeed, double *rvar);
  /// \brief Latin Hypercube Sampling, normal distribution in all directions
  void normLHS(const int zSeed, Array2D<double> &rvar);

  /// \brief Improved Hypercube Sampling, uniform distribution in all directions
  void unifIHS(const int dfac, dsfmt_t *rnstate, const int ndim, const int ns, double *rndnos);
  /// \brief Improved Hypercube Sampling, uniform distribution in all directions
  void unifIHS(const int dfac, dsfmt_t *rnstate, Array2D<double> &rndnos);

  /// \brief Quasi Monte Carlo, Halton sequence
  void getHaltonSeq ( const int nelem, const int dim, std::vector<double> &seq );
  /// \brief Quasi Monte Carlo, Halton sequence
  void getHaltonSeq ( const int nelem, const int dim, const int step,
  		                const std::vector<int> &skip, const std::vector<int> &jump, const std::vector<int> &base,
  		                std::vector<double> &seq);
  /// \brief Quasi Monte Carlo, Hammersley sequence
  void getHammersleySeq ( const int nelem, const int dim, double *seq );

 private:

  /// \brief Permutes a given array in place
  void getPerm(const int nn, const int seed, int* perm);
  /// \brief Retrieves and Improved Hypercube Sampling permutation
  void getIHSperm(const int ndim, const int ns, int *x, const int dupl, dsfmt_t *rnstate);
  /// \brief Returns first nPrime prime numbers
  int getPrimes(const int nPrime, std::vector<int> &primes);

  int hStep_;
  std::vector<int> hSkip_;
  std::vector<int> hJump_;
  std::vector<int> hBase_;

};

#endif /* !SAMPLING_HPP_SEEN */
