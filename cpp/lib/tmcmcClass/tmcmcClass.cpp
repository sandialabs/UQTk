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
// \brief Transitional Markov chain Monte Carlo class

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
#include "tmcmcClass.h"

void TMCMC::runChain(int ncalls, Array1D<double>& chstart){
  // Check mandatory information
  if(!(this -> getDimInit())){
    throw Tantrum((string) "Chain dimensionality needs to be initialized");
  }

  // Check what is not initialized and use defaults instead
  // \todo Specify defaults somewhere more transparently
  if(!tmcmcNprocsInit_)
    this->initTMCMCNprocs(this->default_tmcmc_nprocs_);

  if(!tmcmcGammaInit_)
    this->initTMCMCGamma(this->default_tmcmc_gamma_);

  if(!tmcmcCvInit_)
    this->initTMCMCCv(this->default_tmcmc_cv_);

  if(!tmcmcMFactorInit_)
    this->initTMCMCMFactor(this->default_tmcmc_MFactor_);

  if(!tmcmcBasisInit_)
    this->initTMCMCBasis(this->default_tmcmc_basis_);

  if(!tmcmcCATStepsInit_)
    this->initTMCMCCATSteps(this->default_tmcmc_CATSteps_);

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

  double logevid;

  std::vector<double> samples;
  std::vector<double> logpriors;
  std::vector<double> logliks;
  std::ofstream evidFile("Evidence.dat");

  // Run TMCMC, get evidence
  logevid = tmcmc(samples, logpriors, logliks,
  TMCMCGamma, ncalls,
  this -> getSeed(), TMCMCNprocs, this-> GetChainDim(),
  TMCMCCv, TMCMCMFactor,
  TMCMCBasis, TMCMCCATSteps, this -> getWriteFlag());

  // clean up
  std::ifstream moveFile("tmcmc_moveIntermediateFiles.sh");
  if (moveFile.is_open()) {
    std::string moveStr = "./tmcmc_moveIntermediateFiles.sh TMCMCIntermediates";
    system(moveStr.c_str());
  }

  // Convert samples to chain format
  Array1D<double> sample(this -> GetChainDim());
  // std::cout << "hero1" << std::endl;

  for (int i=0; i < ncalls; i++) {
    this -> setCurrentStateStep(i+1);
    for (int j=0; j<this -> GetChainDim(); j++ ){
      // std::cout << i << " " << j << " " << samples[i*this->chainDim_+j] << std::endl;
      sample(j) = samples[i*this -> GetChainDim()+j];
    }
    this -> setCurrentStateState(sample);
    this -> setCurrentStateAlfa(0.0);
    this -> setCurrentStatePost(logpriors[i]+logliks[i]);
    this -> addCurrentState();
  }

  evidFile << std::setprecision(18) << logevid << std::endl;

  if (this -> getWriteFlag() == 1){
    // Output to file
    if(!strcmp(output.c_str(),"txt")){
      string name = this -> getFilename();
      this->writeChainTxt(name);
    }
    else if(!strcmp(output.c_str(),"bin")){
      string name = this -> getFilename();
      this->writeChainBin(name);
    }
    else
      throw Tantrum((string) "Chain output type is not recognized");
  }

  evidFile.close();
}

void TMCMC::runChain(int ncalls){
  Array1D<double> chstart(this -> GetChainDim(),0.e0);

  this->runChain(ncalls, chstart);
}

void TMCMC::initTMCMCNprocs(int tmcmc_nprocs){
  TMCMCNprocs = tmcmc_nprocs;
  tmcmcNprocsInit_ = true;

  return;
}

void TMCMC::initTMCMCGamma(double tmcmc_gamma){
  TMCMCGamma = tmcmc_gamma;
  tmcmcGammaInit_ = true;

  return;
}

void TMCMC::initTMCMCCv(double tmcmc_cv){
  TMCMCCv = tmcmc_cv;
  tmcmcCvInit_ = true;

  return;
}

void TMCMC::initTMCMCMFactor(int tmcmc_MFactor){
  TMCMCMFactor = tmcmc_MFactor;
  tmcmcMFactorInit_ = true;

  return;
}

void TMCMC::initTMCMCBasis(bool tmcmc_basis){
  TMCMCBasis = tmcmc_basis;
  tmcmcBasisInit_ = true;

  return;
}

void TMCMC::initTMCMCCATSteps(int tmcmc_CATSteps){
  tmcmcCATStepsInit_ = true;
  TMCMCCATSteps = tmcmc_CATSteps;

  return;
}

int TMCMC::getTMCMCNprocs(){
  return TMCMCNprocs;
}

double TMCMC::getTMCMCGamma(){
  return TMCMCGamma;
}

double TMCMC::getTMCMCCv(){
  return TMCMCCv;
}

int TMCMC::getTMCMCMFactor(){
  return TMCMCMFactor;
}

bool TMCMC::getTMCMCBasis(){
  return TMCMCBasis;
}

int TMCMC::getTMCMCCATSteps(){
  return TMCMCCATSteps;
}
