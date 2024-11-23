/* =====================================================================================

                      The UQ Toolkit (UQTk) version 3.1.5
                          Copyright (2024) NTESS
                        https://www.sandia.gov/UQToolkit/
                        https://github.com/sandialabs/UQTk

     Copyright 2024 National Technology & Engineering Solutions of Sandia, LLC (NTESS).
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
// \file amcmc.cpp
// \author K. Sargsyan, C. Safta, B. Debusschere, 2012
// \brief Adaptive Markov Chain Monte Carlo class

#include <math.h>
#include <float.h>
#include "error_handlers.h"
#include "deplapack.h"

#include "tools.h"
#include "arrayio.h"
#include "arraytools.h"
#include "mcmc.h"
#include "gen_defs.h"
#include "lbfgs_routines.h"
#include "amcmc.h"

void AMCMC::printChainSetup(){
  if (this->gammaInit_)
    cout << "Gamma            : " << this->gamma << endl;
  else
    cout << "Gamma (default)  : " << this->default_gamma_ << endl;
  if (this->epscovInit_)
    cout << "Eps_Cov          : " << this->eps_cov << endl;
  else
    cout << "Eps_Cov (default): " << this->default_eps_cov_ << endl;

  return;
}

void AMCMC::runChain(int ncalls, Array1D<double>& chstart){
  // Check mandatory information
  if(!(this -> getDimInit())){
    throw Tantrum((string) "Chain dimensionality needs to be initialized");
  }

  // Check what is not initialized and use defaults instead
  // \todo Specify defaults somewhere more transparently

  // Set defaults proposal covariance
  if(!(this -> getPropCovInit())){
    Array1D<double> chsig(this -> GetChainDim(),0.e0);
    for(int i=0;i<this -> GetChainDim();i++) chsig(i)=MAX(fabs(0.1*chstart(i)),0.001);
    this->initChainPropCovDiag(chsig);
  }

  // Set defaults output format
  if (!(this -> getOutputInit())){
    this->setOutputInfo("txt","chain.dat", max(1,(int) ncalls/100), max(1,(int) ncalls/20));
  }

  // Set the default parameters for aMCMC
  if(!adaptstepInit_){
    this->initAdaptSteps((int) ncalls/10,10,ncalls);
  }
  if(!gammaInit_){
    this->initAMGamma(this->default_gamma_);
  }
  if(!epscovInit_){
    this->initEpsCov(this->default_eps_cov_);
  }

  // Work variables for simplicity
  string output = this -> getOutputType();

  // Initial chain state
  this -> setCurrentStateStep(0);
  this -> setCurrentStateState(chstart);
  this -> setCurrentStateAlfa(0.0);
  this -> setCurrentStatePost(this->evalLogPosterior(chstart));
  this->updateMode();
  this -> addCurrentState();
  //fullChain_.PushBack(currState_);

  // Number of accepted steps and number of all trials
  int nacc=0;
  int nall=0;

  // No new mode found yet
  this -> setNewMode(false);

  for(int t = 1; t <= ncalls; ++t){
    this -> setCurrentStateStep(t);
    double sum_alpha=0.0;

    // Create a new proposed sample
    for(int is = 0; is < nSubSteps_; ++is){
      Array1D<double> m_cand;
      Array1D<double> state;
      this -> getCurrentStateState(state);
      this -> proposal(state, m_cand, t);

      // Evaluate the posterior at the new sample point
      double eval_cand = this->evalLogPosterior(m_cand);

      // Evaluate the new|old and old|new proposals
      Array1D<double> state1;
      this -> getCurrentStateState(state1);
      double old_given_new = this->probOldNew(state1, m_cand);
      double new_given_old = this->probOldNew(m_cand,state1);

      // Accept or reject it
      double post = this -> getCurrentStatePost();
      double alpha = exp(eval_cand - post + old_given_new - new_given_old);
      sum_alpha += alpha;
      if (this->inDomain(m_cand) && (alpha>=1 || alpha > dsfmt_genrand_urv(&RandomState))) { // Accept and update the state
        nacc++;
        this -> setCurrentStateState(m_cand);
        this -> setCurrentStatePost(eval_cand);
        if (this->getFcnAcceptInit())
          this -> runAcceptFcn();
      } // If state not accepted, keep previous state as the current state
      else{
        if (this->getFcnRejectInit())
          this -> runRejectFcn();
      }

      ++nall;
    }

    double alfa = sum_alpha / nSubSteps_;
    this -> setCurrentStateAlfa(alfa);

    // Append the current state to the array of all past states
    this -> addCurrentState();

    // Keep track of the mode (among the locations visited so far)
    // \todo maybe only store tmode_(we save the full chain anyway)
    if (this -> getCurrentStatePost() > this -> getModeStatePost()){
      this->updateMode();
      this -> setNewMode(true);
    }

    this -> setAcceptRatio((double) nacc/nall);

    if(this -> getWriteFlag() == 1){
      // Output to Screen
      if( t % this -> getScreenFreq() == 0 || t==ncalls){

        printf("%lg %% completed; Chain step %d\n", 100.*t/ncalls,t);
        printf("================= Current logpost:%f, Max logpost:%f, Accept rate:%f\n",this -> getCurrentStatePost(),this -> getModeStatePost(),this -> getAcceptRatio());
        printf("================= Current MAP params: ");
        Array1D<double> state;
        this -> getModeStateState(state);
        for(int ic=0;ic<this -> GetChainDim();ic++){
          printf("par(%d)=%f ",ic,state(ic));
        }
        cout << endl;

      }

        // Output to File
      if( t % this -> getFileFreq() == 0 || t==ncalls){

        if(!strcmp(output.c_str(),"txt"))
          this->writeChainTxt(this -> getFilename());
        else  if(!strcmp(output.c_str(),"bin"))
          this->writeChainBin(this -> getFilename());
        else
          throw Tantrum((string) "Chain output type is not recognized");
        this -> setLastWrite(t);
      }
    }
  }

  return;
}

void AMCMC::runChain(int ncalls){
  Array1D<double> chstart(this -> GetChainDim(),0.e0);

  this->runChain(ncalls, chstart);
}

void AMCMC::proposal(Array1D<double>& m_t,Array1D<double>& m_cand,int t){
  int chol_info=0;
  char lu='L';

  // xm[] is the mean of m_t[] over all previous states, X_0,...,X_{t-1}
  // at this stage, index t, we know X_0,X_1,...,X_{t-1}
  // and we're seeking to find X_t, the new state of the chain
  // also evaluate covt, the covariance matrix

  if(t == 1){
    // at the first iteration, the mean is easy and the sample covariance is 0
    curmean = m_t;
    curcov.Resize(this -> GetChainDim(),this -> GetChainDim(),0.e0);
  }
  else if(t > 1 && t < adaptstep_(2)){
    for (int i = 0; i < this -> GetChainDim(); ++i) {
      curmean(i)  = ( curmean(i)*(t-1.) + m_t(i) )/t;
    }

    for(int i = 0; i < this -> GetChainDim(); ++i){
      for(int j = 0; j < i + 1; ++j){
        curcov(i,j) = ( (t-2.)/(t-1.) ) * curcov(i,j) + ( t/((t-1.)*(t-1.)) ) * ( m_t(i) - curmean(i) )*( m_t(j) - curmean(j) );
      }
    }

    for (int i = 0; i < this -> GetChainDim(); i++){
      for(int j = i + 1; j < this -> GetChainDim(); ++j){
        curcov(i,j) = curcov(j,i);
      }
    }
  }

  // Jump size
  double sigma = gamma * 2.4 * 2.4 / (double) this -> GetChainDim();

  if(t == 1){
    //propLCov_=chcov;
    this -> getChainPropCov(propLCov_);

    // Cholesky factorization of the proposal covariance propLCov_, done in-place
    // Note, for diagonal covariances, this is an overkill
    int chdim = this -> GetChainDim();
    FTN_NAME(dpotrf)(&lu,&chdim, propLCov_.GetArrayPointer(),&chdim,&chol_info);
  }

  if ( ( t > adaptstep_(0) ) && ( (t % adaptstep_(1) ) ==  0 ) && t <= adaptstep_(2) ){
    for (int i = 0; i < this -> GetChainDim(); ++i){
      for(int j = 0; j < this -> GetChainDim(); ++j){
        propLCov_(i,j) = sigma*(curcov(i,j) + (i==j) * eps_cov ) ;
      }
    }

    //chcov = propLCov_;
    this -> initChainPropCov(propLCov_);

    // Cholsky factorization of the proposal covariance propLCov_, done in-place
    int chdim = this -> GetChainDim();
    FTN_NAME(dpotrf)(&lu,&chdim, propLCov_.GetArrayPointer(),&chdim,&chol_info);

    // Catch the error in Cholesky factorization
    if (chol_info != 0 ) {
      printf("Error in Cholesky factorization, info=%d, printing the matrix below:\n", chol_info);

      for(int i=0;i< this -> GetChainDim();i++){
        for(int j=0;j< this -> GetChainDim();j++)
          printf("%lg ",propLCov_(i,j));
        printf("\n");
      }

      exit(1);
    }
  }

  // Candidate state is a multivariate normal sample away from the current state
  m_cand=m_t;
  Array1D<double> xi(this -> GetChainDim(),0.e0);
  for (int i=0; i < this -> GetChainDim(); i++) {
    xi(i)=dsfmt_genrand_nrv(&RandomState);
    double Lnrv=0.0;
    for (int j=0; j < i+1; j++) {
      Lnrv += propLCov_(i,j)*xi(j);
    }
    m_cand(i) += Lnrv;
  }

  return;

}

void AMCMC::initAdaptSteps(int adaptstart,int adaptstep, int adaptend){
  this -> adaptstep_.PushBack(adaptstart);
  this -> adaptstep_.PushBack(adaptstep);
  this -> adaptstep_.PushBack(adaptend);

  adaptstepInit_ = true;

  return;
}

void AMCMC::initAMGamma(double gamma_){
  gamma = gamma_;
  gammaInit_ = true;

  return;
}

void AMCMC::initEpsCov(double eps_cov_){
  eps_cov = eps_cov_;
  epscovInit_ = true;

  return;
}

void AMCMC::getAdaptSteps(Array1D<int> adaptstep){
  this -> adaptstep_ = adaptstep;

  return;
}

double AMCMC::getGamma(){
  return gamma;
}

double AMCMC::getEpsCov(){
  return eps_cov;
}
