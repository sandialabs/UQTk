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
// \author K. Sargsyan, C. Safta, B. Debusschere, 2012
// \brief Single Site Markov Chain Monte Carlo class

#include <math.h>
#include <float.h>
#include "error_handlers.h"
#include "deplapack.h"

#include "tools.h"
#include "arrayio.h"
#include "arraytools.h"
#include "mcmc.h"
#include "tmcmc.h"
#include "gen_defs.h"
#include "lbfgs_routines.h"
#include "ss.h"

void SS::runChain(int ncalls, Array1D<double>& chstart){
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

  // Work variables for simplicity
  string output = this -> getOutputType();

  // Initial chain state
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
      this -> proposal(state, m_cand, is);

      // Evaluate the posterior at the new sample point
      double eval_cand = this->evalLogPosterior(m_cand);

      // Evaluate the new|old and old|new proposals
      /// \todo See notes in AM and MALA
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

void SS::runChain(int ncalls){
  Array1D<double> chstart(this -> GetChainDim(),0.e0);

  this->runChain(ncalls, chstart);
}

void SS::proposal(Array1D<double>& m_t,Array1D<double>& m_cand,int dim)
{
  // Single-site proposal
  m_cand = m_t;

  Array2D<double> chcov;

  this -> getChainPropCov(chcov);

  m_cand(dim) += ( sqrt(chcov(dim,dim))*dsfmt_genrand_nrv(&RandomState) );

  return;
}

int SS::getNSubSteps(){
  return this -> nSubSteps_;
}
