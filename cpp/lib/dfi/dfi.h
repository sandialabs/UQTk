/* =====================================================================================
                     The UQ Toolkit (UQTk) version 3.0.4
                     Copyright (2017) Sandia Corporation
                     http://www.sandia.gov/UQToolkit/

     Copyright (2013) Sandia Corporation. Under the terms of Contract DE-AC04-94AL85000
     with Sandia Corporation, the U.S. Government retains certain rights in this software.

     This file is part of The UQ Toolkit (UQTk)

     UQTk is free software: you can redistribute it and/or modify
     it under the terms of the GNU Lesser General Public License as published by
     the Free Software Foundation, either version 3 of the License, or
     (at your option) any later version.

     UQTk is distributed in the hope that it will be useful,
     but WITHOUT ANY WARRANTY; without even the implied warranty of
     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
     GNU Lesser General Public License for more details.

     You should have received a copy of the GNU Lesser General Public License
     along with UQTk.  If not, see <http://www.gnu.org/licenses/>.

     Questions? Contact Bert Debusschere <bjdebus@sandia.gov>
     Sandia National Laboratories, Livermore, CA, USA
===================================================================================== */
#ifndef UQTKDFI_H_SEEN
#define UQTKDFI_H_SEEN


#include "dsfmt_add.h"
#include "mcmc.h"

class DFISetupBase{
public:

  virtual double f(Array1D<double>&, Array1D<double>&, Array1D<double>&){};
  virtual double S(Array1D<MCMC::chainstate> inner_samples){};
};

class DFIInner: public LikelihoodBase{
public:
  DFIInner(DFISetupBase& d);
  DFIInner(){};
  ~DFIInner(){};

  // params for running the inner chain
  int ndim, nBurn, nCalls; // dim of inner chain
  Array1D<double> beta0_;  // initial start
  Array1D<double> gammas_; // initial proposal covariance
  Array1D<MCMC::chainstate> samples_; // holds beta samples from inner chain

  DFISetupBase* d_;        // class which holds logposterior
  Array1D<double> z_;     // internal z data (change in outer loop)
  Array1D<double> sigma_; // internal sigma (change in outer loop)

  // variables for holding mean and variance
  Array1D<double> means_; 
  Array1D<double> stds_;
  Array1D<double> quants;
  
  double eval(Array1D<double>&); // evaluate logposterior
  void getSamples(); // get Array1D of samples

  // params for summary statistics
  Array1D<double> means0; 
  Array1D<double> stds0; 
  double delta1, delta2; 
  double S(); // comparison of summary statistics

};

class DFI: public LikelihoodBase{
public:
  int zdim; 
  int sdim; 
  int nBeta;  

  DFIInner* dfi_inner_;
  int nCalls; 

  DFI(DFIInner&);
  ~DFI(){};

  // params for eval function
  Array1D<double> z_; 
  Array1D<double> sigma_; 

  double eval(Array1D<double>&);
  void runChain(int nCalls, Array1D<double> gammas, Array1D<double> start, int seed, int node);

  void getMLE(Array1D<double>& xstart);
};

#endif /* UQTKDFI_H_SEEN */
