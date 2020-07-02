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
// \file mcmc.cpp
// \author K. Sargsyan, C. Safta, B. Debusschere, 2012 -
// \brief Markov chain Monte Carlo class

#include <math.h>
#include <float.h>
#include "error_handlers.h"
#include "deplapack.h"

#include "tools.h"
#include "arrayio.h"
#include "arraytools.h"
#include "mcmc-refactor.h"
#include "tmcmc.h"
#include "gen_defs.h"
#include "lbfgs_routines.h"

double neg_logposteriorproxy(int chaindim, double* m, void* classpointer);
void grad_neg_logposteriorproxy(int chaindim, double* m, double* grads, void* classpointer);

MCMC::MCMC(double (*logposterior)(Array1D<double>&, void *), void *postinfo){
  // Set Flag
  FLAG = 0;

  // Setting entering the pointers to the Log Posterior function as well as the data function
  postInfo_ = postinfo;
  logPosterior_ = logPosterior;

  // Initiate the random number generator seed
  // \todo This needs to be made more generic
  seed_ = 13;
  dsfmt_init_gen_rand(&RandomState,seed_);

  // Set the location of the last chain state written (-1 means nothing is written to files yet)
  lastwrite_ = -1;

  // Set Write Flag
  WRITE_FLAG = 1;

  return;
}

MCMC::MCMC(LikelihoodBase& L){
  FLAG = 1;
  L_ = &L;

  // Initiate the random number generator seed
  // \todo This needs to be made more generic
  seed_ = 13;
  dsfmt_init_gen_rand(&RandomState,seed_);

  // Set the location of the last chain state written (-1 means nothing is written to files yet)
  lastwrite_=-1;

  // Set Write Flag
  WRITE_FLAG = 1;

  return;
}

MCMC::MCMC(){
  FLAG = 2;

  // Initiate the random number generator seed
  // \todo This needs to be made more generic
  seed_ = 13;
  dsfmt_init_gen_rand(&RandomState,seed_);

  // Set the location of the last chain state written (-1 means nothing is written to files yet)
  lastwrite_=-1;

  // Set Write Flag
  WRITE_FLAG = 1;

  return;
}

void MCMC::setWriteFlag(int I){
  WRITE_FLAG = I;
}

void MCMC::setGradient(void (*gradlogPosterior)(Array1D<double>&, Array1D<double>&, void *)){
  gradlogPosterior_ = gradlogPosterior;
  gradflag_ = true;
  return;
}

void MCMC::setMetricTensor(void (*metricTensor)(Array1D<double>&, Array2D<double>&, void *)){
  metricTensor_ = metricTensor;
  tensflag_ = true;
  return;
}

void MCMC::setFcnAccept(void (*fcnAccept)(void *))
{
  fcnAccept_ = fcnAccept;
  fcnAcceptFlag_ = true;
  return;
}

void MCMC::setFcnReject(void (*fcnReject)(void *))
{
  fcnReject_ = fcnReject;
  fcnRejectFlag_ = true;
  return;
}

void MCMC::setChainDim(int chdim){
  chainDim_ = chdim;
  chaindimInit_ = true;
  return;
}

void MCMC::initChainPropCov(Array2D<double>& propcov){
  // Initialize the proposal covariance matrix
  methodinfo_.chcov=propcov;
  // Set the initialization flag to True
  propcovInit_=true;

  return;
}

