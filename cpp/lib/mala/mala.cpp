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
// \file mala.cpp
// \author K. Sargsyan, C. Safta, B. Debusschere, 2012
// \brief MALA Markov chain Monte Carlo class

#include <math.h>
#include <float.h>
#include "error_handlers.h"
#include "deplapack.h"
#include "mala.h"
#include "tools.h"
#include "arrayio.h"
#include "arraytools.h"
#include "mcmc.h"
#include "gen_defs.h"
#include "lbfgs_routines.h"

double neg_logposteriorproxy(int chaindim, double* m, void* classpointer);
void grad_neg_logposteriorproxy(int chaindim, double* m, double* grads, void* classpointer);

void MALA::setGradient(void (*gradlogPosterior)(Array1D<double>&, Array1D<double>&, void *)){
  gradlogPosterior_ = gradlogPosterior;
  gradflag_ = true;
  return;
}

void MALA::getGradient(void (*gradlogPosterior)(Array1D<double>&, Array1D<double>&, void *)){
  gradlogPosterior = gradlogPosterior_;
  return;
}

bool MALA::getGradientFlag(){
  return gradflag_;
}

void MALA::runOptim(Array1D<double>& start){
  int n=start.Length();
  int m=5;
  Array1D<int> nbd(n,0);
  Array1D<double> l(n,0.e0);
  Array1D<double> u(n,0.e0);

  for (int i=0;i<n;i++){
    if (this -> getLowerFlag(i) && !(this -> getUpperFlag(i))){
      nbd(i)=1;
      l(i)=this->getLower(i);
    }
    if (this -> getUpperFlag(i) && !(this -> getLowerFlag(i))){
      nbd(i)=3;
      u(i)=this->getUpper(i);
    }
    if (this -> getUpperFlag(i) && this -> getLowerFlag(i)){
      nbd(i)=2;
      l(i)=this->getLower(i);
      u(i)=this->getUpper(i);
    }
  }

  void* info=this;

  if (gradflag_)
    lbfgsDR(n,m,start.GetArrayPointer(),nbd.GetArrayPointer(),l.GetArrayPointer(),u.GetArrayPointer(),neg_logposteriorproxy,grad_neg_logposteriorproxy,info) ;
  else
    lbfgsDR(n,m,start.GetArrayPointer(),nbd.GetArrayPointer(),l.GetArrayPointer(),u.GetArrayPointer(),neg_logposteriorproxy,NULL,info) ;

  this -> setCurrentStateStep(0);
  this -> setCurrentStateState(start);
  this -> setCurrentStateAlfa(0.0);
  this -> setCurrentStatePost(this->evalLogPosterior(start));
  this->updateMode();

  return;
}

void MALA::runChain(int ncalls, Array1D<double>& chstart){
  // Check mandatory information
  if(!(this -> getDimInit())){
    throw Tantrum((string) "Chain dimensionality needs to be initialized");
  }

  // Check mandatory information specific to mala
  if(!(this -> getGradientFlag())){
    throw Tantrum((string) "Gradient function needs to be initialized");
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

  if(!epsMalaInit_){
    this->initEpsMALA(default_eps_mala_);
  }

  // Work variables for simplicity
  string output = this -> getOutputType();

  // Initial chain state
  chainstate state;
  this -> setCurrentStateStep(0);
  this -> setCurrentStateState(chstart);
  this -> setCurrentStateAlfa(0.0);
  this -> setCurrentStatePost(this->evalLogPosterior(chstart));
  this->updateMode();
  this -> addCurrentState();

  // Number of accepted steps and number of all trials
  int nacc=0;
  int nall=0;

  // No new mode found yet
  this -> setNewMode(false);

  for(int t = 1; t < ncalls; ++t){
    this -> setCurrentStateStep(t);
    double sum_alpha=0.0;

    // Create a new proposed sample
    for(int is = 0; is  < nSubSteps_; ++is){
      Array1D<double> m_cand;
      Array1D<double> state;
      this -> getCurrentStateState(state);
      this -> proposal(state, m_cand);

      // Evaluate the posterior at the new sample point
      double eval_cand = this->evalLogPosterior(m_cand);

      // Evaluate the new|old and old|new proposals
      /// \todo See the previous note in proposal for AM
      double old_given_new = this->probOldNew(state, m_cand);
      double new_given_old = this->probOldNew(m_cand,state);

      // Accept or reject it
      double alpha = exp(eval_cand - this -> getCurrentStatePost() + old_given_new - new_given_old);
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
        for(int ic=0;ic<this -> GetChainDim();ic++)
          printf("par(%d)=%f ",ic,state(ic));
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

void MALA::runChain(int ncalls){
  Array1D<double> chstart(this -> GetChainDim(),0.e0);

  this->runChain(ncalls, chstart);
}

void MALA::evalGradLogPosterior(Array1D<double>& m, Array1D<double>& grads){
  // Evaluate given the log-posterior function defined by the user in the constructor
  void *post;
  this -> getPostInfo(post);
  gradlogPosterior_(m,grads,post);

  return;
}

void MALA::proposal(Array1D<double>& m_t,Array1D<double>& m_cand){
  Array1D<double> grads;
  gradlogPosterior_(m_t,grads,NULL);
  cout << "grads= " << grads(0) << " " << grads(1) << endl;
  m_cand = m_t;
  for (int i=0; i < this -> GetChainDim(); i++) {
    m_cand(i) += eps_mala * eps_mala *grads(i)/2.;
    m_cand(i) += eps_mala * dsfmt_genrand_nrv(&RandomState);
  }

  return;
}

double MALA::probOldNew(Array1D<double>& a, Array1D<double>& b){
  double logprob;
  Array1D<double> gradb;

  ///\todo figure out this because the grad log posterior is a private variable in MCMC, might need to put it in the MALA class
  gradlogPosterior_(b,gradb,NULL);
  double eps2 = eps_mala * eps_mala;
  Array1D<double> bmean(this -> GetChainDim(),0.e0);
  Array1D<double> diagcov(this -> GetChainDim(),0.e0);

  for (int i=0;i<this -> GetChainDim();i++){
    bmean(i)=b(i)+eps2*gradb(i)/2.0;
    diagcov(i)=eps2;
  }

  logprob=evallogMVN_diag(a,bmean,diagcov);

  return logprob;
}

double MALA::evallogMVN_diag(Array1D<double>& x,Array1D<double>& mu,Array1D<double>& sig2){
  double pi=4.0*atan(1.0);

  double value=0.e0;

  // \todo Put sanity checks on dimensions

  for (int i=0;i<this->GetChainDim();i++){
    value -= 0.5*log(2.*pi*sig2(i));
    value -= (x(i)-mu(i))*(x(i)-mu(i))/(2.0*sig2(i));
  }
  return value;
}

void MALA::initEpsMALA(double eps_mala_){
  this -> eps_mala = eps_mala_;
  epsMalaInit_ = true;

  return;
}

double MALA::getEpsMALA(){
  return this -> eps_mala;
}

void grad_neg_logposteriorproxy(int chaindim, double* m, double* grads, void* classpointer){
  MALA* thisClass=(MALA*) classpointer;

  if(chaindim != thisClass->GetChainDim()){
    throw Tantrum(std::string("neg_logposteriorproxy: The passed in MCMC chain dimension does not match the  dimension of the MChain class instance"));
  }

  Array1D<double> mm(chaindim,0.e0);

  for(int i = 0; i < chaindim; ++i){
    mm(i) = m[i];
  }

  // Call the posterior function and return its result
  Array1D<double> grads_arr;
  thisClass->evalGradLogPosterior(mm, grads_arr);

  for(int i=0;i<chaindim;i++)
    grads[i]=-grads_arr(i);

  return;
}

double neg_logposteriorproxy(int chaindim, double* m, void* classpointer){
  MCMC* thisClass = (MCMC*) classpointer;

  if(chaindim != thisClass -> GetChainDim()){
    throw Tantrum(std::string("neg_logposteriorproxy: The passed in MCMC chain dimension does not match the  dimension of the MChain class instance"));
  }

  Array1D<double> mm(chaindim,0.e0);

  for(int i = 0; i < chaindim; ++i){
    mm(i) = m[i];
  }

  return -thisClass -> evalLogPosterior(mm);
}
