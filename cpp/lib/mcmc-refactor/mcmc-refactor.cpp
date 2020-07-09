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
// \file mcmc-refactor.cpp
// \author K. Sargsyan, C. Safta, B. Debusschere, 2012 - L. Boll, 2020
// \brief Markov chain Monte Carlo class and its subsequent derived classes

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

/*MCMC::MCMC(LikelihoodBase& L){
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
}*/

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
  chcov=propcov;
  // Set the initialization flag to True
  propcovInit_=true;

  return;
}

void MCMC::initChainPropCovDiag(Array1D<double>& sig){
  // Create a diagonal matrix and fill in the diagonal terms
  chcov.Resize(chainDim_,chainDim_,0.e0);
  for(int i = 0; i < chainDim_; ++i){
    chcov(i,i)=sig(i)*sig(i);
  }

  // Set the initialization flag to True
  propcovInit_=true;

  return;
}

void MCMC::setOutputInfo(string outtype, string file,int freq_file, int freq_screen){
  outputinfo_.type = outtype;
  outputinfo_.filename = file;
  outputinfo_.freq_chainfile = freq_file;
  outputinfo_.freq_outscreen = freq_screen;
  // Set the initialization flag to True
  outputInit_ = true;

  return;
}

void MCMC::namesPrepended(){
  namesPrepend = true;

  return;
}

void MCMC::setSeed(int seed){
  seed_ = seed;
  dsfmt_init_gen_rand(&RandomState,seed);
  return;
}

void MCMC::setLower(double lower, int i){
  Lower_(i) = lower;
  lower_flag_(i) = 1;

  return;
}

void MCMC::setUpper(double upper, int i){
  Upper_(i) = lower;
  upper_flag_(i) = 1;

  return;
}

void MCMC::setDefaultDomain(){
  Lower_.Resize(chainDim_,-DBL_MAX);
  Upper_.Resize(chainDim_,DBL_MAX);
  lower_flag_.Resize(chainDim_,0);
  upper_flag_.Resize(chainDim_,0);

  return;
}

void MCMC::getChainPropCov(Array2D<double>& propcov){
  // Get the proposal covariance matrix
  propcov = chcov;

  return;
}

string MCMC::getFilename(){
  return outputinfo_.filename;
}

int MCMC::getWriteFlag(){
  return WRITE_FLAG;
}

void MCMC::getSamples(int burnin, int every, Array2D<double>& samples){
  int nCalls = fullChain_.Length();
  samples.Resize(chainDim_,0); // initialize sample array
  int j=0;
  for (int i = burnin; i < nCalls; i+=every){
    samples.insertCol(fullChain_(i).state,j);
    j++;
  }

  return;
}

void MCMC::getSamples(Array2D<double>& samples){
  getSamples(0,1,samples);
}

void MCMC::getGradient(void (*gradlogPosterior)(Array1D<double>&, Array1D<double>&, void *)){
  gradlogPosterior = gradlogPosterior_;
  return;
}

void MCMC::getMetricTensor(void (*metricTensor)(Array1D<double>&, Array2D<double>&, void *)){
  metricTensor = metricTensor_;
  return;
}

void MCMC::getFcnAccept(void (*fcnAccept)(void *)){
  fcnAccept = fcnAccept_;
  return;
}

void MCMC::getFcnAccept(void (*fcnReject)(void *)){
  fcnReject = fcnReject_;
  return;
}

string MCMC::getOutputType(){
  return outputinfo_.outtype;
}

int MCMC::getFileFreq(){
  return outputinfo_.freq_file;
}

int MCMC::getScreenFreq(){
  return outputinfo_.freq_screen;
}

bool MCMC::getNamesPrepended(){
  return namesPrepend;
}

int MCMC::getSeed(){
  return seed_;
}

double MCMC::getLower(int i){
  return Lower_(i)
}

double MCMC::getUpper(int i){
  return Upper_(i);
}

bool MCMC::getDimInit(){
  return chaindimInit_;
}

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

void MCMC::resetChainState(){
  fullChain_.Clear();

  return;
}

void MCMC::resetChainFilename(string filename){
  this -> resetChainState();
  lastwrite_ = -1;
  outputinfo_.filename = filename;

  return;
}

void MCMC::parseBinChain(string filename, Array1D<chainstate>& readchain){
  double tmp;
  int readstep,i=0;

  FILE *fb = fopen(filename.c_str(),"rb");

  chainstate curchain;
  // Read the binary file and write to an array of chain states
  while( fread(&readstep,sizeof(int),1,fb) ) {

    curchain.step=readstep;
    assert(readstep==i);

    curchain.state.Resize(chainDim_);

    fread(curchain.state.GetArrayPointer(),chainDim_*sizeof(double),1,fb);
    fread(&tmp,sizeof(double),1,fb);
    curchain.alfa=tmp;
    fread(&tmp,sizeof(double),1,fb);
    curchain.post=tmp;
    readchain.PushBack(curchain);
    i++;
  }

  fclose(fb);

  return;
}

void MCMC::writeFullChainTxt(string filename, Array1D<chainstate> fullchain){
  // Open the text file in a write mode
  char* writemode = "w";

  // Append if the names already prepended
  if (namePrepend_)
    writemode = "a";

  FILE* f_out;
  if(!(f_out = fopen(filename.c_str(),writemode))){
    printf("writeChain: could not open file '%s'\n",filename.c_str());
    exit(1);
  }

  // Write the full chain
  for(int i=0;i<fullchain.XSize();i++){
    fprintf(f_out, "%d ", fullchain(i).step);
    for(int ic=0;ic<chainDim_;ic++)
      fprintf(f_out, "%24.16lg ", fullchain(i).state(ic));
    fprintf(f_out, "%24.16lg %24.16lg \n", fullchain(i).alfa,fullchain(i).post);
  }

  // Closte the file
 if(fclose(f_out)){
   printf("writeChain: could not close file '%s'\n",filename.c_str());
   exit(1);
 }

 // Report
 printf("Written the full chain out to text file %s\n",filename.c_str());

 return;
}

void MCMC::getFullChain(Array1D<chainstate>& readchain){
  readchain = fullChain_;
  return;
}

void MCMC::appendMAP(){
  fullChain_.PushBack(modeState_);
  return;
}

double MCMC::getMode(Array1D<double>& MAPparams){
  MAPparams = modeState_.state;
  return modeState_.post;
}

void MCMC::runOptim(Array1D<double>& start){
  int n=start.Length();
  int m=5;
  Array1D<int> nbd(n,0);
  Array1D<double> l(n,0.e0);
  Array1D<double> u(n,0.e0);

  for (int i=0;i<n;i++){
    if (lower_flag_(i) && !upper_flag_(i)){
      nbd(i)=1;
      l(i)=this->Lower_(i);
    }
    if (upper_flag_(i) && !lower_flag_(i)){
      nbd(i)=3;
      u(i)=this->Upper_(i);
    }
    if (upper_flag_(i) && lower_flag_(i)){
      nbd(i)=2;
      l(i)=this->Lower_(i);
      u(i)=this->Upper_(i);
    }
  }

  void* info=this;

  if (gradflag_)
    lbfgsDR(n,m,start.GetArrayPointer(),nbd.GetArrayPointer(),l.GetArrayPointer(),u.GetArrayPointer(),neg_logposteriorproxy,grad_neg_logposteriorproxy,info) ;
  else
    lbfgsDR(n,m,start.GetArrayPointer(),nbd.GetArrayPointer(),l.GetArrayPointer(),u.GetArrayPointer(),neg_logposteriorproxy,NULL,info) ;

  currState_.step=0;
  currState_.state=start;
  currState_.alfa=0.0;
  currState_.post=this->evalLogPosterior(start);
  this->updateMode();

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
  if(!propcovInit_){
    Array1D<double> chsig(this -> GetChainDim(),0.e0);
    for(int i=0;i<this -> GetChainDim();i++) chsig(i)=MAX(fabs(0.1*chstart(i)),0.001);
    this->initChainPropCovDiag(chsig);
  }

  // Set defaults output format
  if (!outputInit_){
    this->setOutputInfo("txt","chain.dat", max(1,(int) ncalls/100), max(1,(int) ncalls/20));
  }

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
  string output=outputinfo_.type;

  // Initial chain state
  currState_.step=0;
  currState_.state=chstart;
  currState_.alfa=0.0;
  currState_.post=this->evalLogPosterior(chstart);
  this->updateMode();
  fullChain_.PushBack(currState_);

  // Number of accepted steps and number of all trials
  int nacc=0;
  int nall=0;

  // No new mode found yet
  newMode_=false;

  for(int t = 1; t < nCalls, ++t){
    currState_.step=t;
    double sum_alpha=0.0;

    // Create a new proposed sample
    for(int is = 0; is  < nSubSteps_; ++is){
      Array1D<double> m_cand;
      this -> proposal(currState_.state, m_cand, t);

      // Evaluate the posterior at the new sample point
      double eval_cand = this->evalLogPosterior(m_cand);

      // Evaluate the new|old and old|new proposals
      double old_given_new = this->probOldNew(currState_.state, m_cand);
      double new_given_old = this->probOldNew(m_cand,currState_.state);

      // Accept or reject it
      double alpha = exp(eval_cand - currState_.post + old_given_new - new_given_old);
      sum_alpha += alpha;
      if (this->inDomain(m_cand) && (alpha>=1 || alpha > dsfmt_genrand_urv(&RandomState))) { // Accept and update the state
        nacc++;
        currState_.state = m_cand;
        currState_.post = eval_cand;
        if (this->fcnAcceptFlag_)
          this->fcnAccept_(this->postInfo_);
      } // If state not accepted, keep previous state as the current state
      else{
        if (this->fcnRejectFlag_)
          this->fcnReject_(this->postInfo_);
      }

      ++nall;
    }

    currState_.alfa = sum_alpha / nSubSteps_;

    // Append the current state to the array of all past states
    fullChain_.PushBack(currState_);

    // Keep track of the mode (among the locations visited so far)
    // \todo maybe only store tmode_(we save the full chain anyway)
    if (currState_.post > modeState_.post){
      this->updateMode();
      newMode_=true;
    }

    accRatio_ = (double) nacc/nall;

    if(WRITE_FLAG == 1){
      // Output to Screen
      if( t % outputinfo_.freq_outscreen == 0 || t==ncalls){

        printf("%lg %% completed; Chain step %d\n", 100.*t/ncalls,t);
        printf("================= Current logpost:%f, Max logpost:%f, Accept rate:%f\n",currState_.post,modeState_.post,accRatio_);
        printf("================= Current MAP params: ");
        for(int ic=0;ic<this -> GetChainDim();ic++)
          printf("par(%d)=%f ",ic,modeState_.state(ic));
        cout << endl;

      }

        // Output to File
      if( t % outputinfo_.freq_chainfile == 0 || t==ncalls){

        if(!strcmp(output.c_str(),"txt"))
          this->writeChainTxt(outputinfo_.filename);
        else  if(!strcmp(output.c_str(),"bin"))
          this->writeChainBin(outputinfo_.filename);
        else
          throw Tantrum((string) "Chain output type is not recognized");
        lastwrite_ = t;
      }
    }
  }

  return;
}

void MALA::runChain(int ncalls, Array1D<double>& chstart){
  // Check mandatory information
  if(!(this -> getDimInit())){
    throw Tantrum((string) "Chain dimensionality needs to be initialized");
  }

  // Check what is not initialized and use defaults instead
  // \todo Specify defaults somewhere more transparently

  // Set defaults proposal covariance
  if(!propcovInit_){
    Array1D<double> chsig(this -> GetChainDim(),0.e0);
    for(int i=0;i<this -> GetChainDim();i++) chsig(i)=MAX(fabs(0.1*chstart(i)),0.001);
    this->initChainPropCovDiag(chsig);
  }

  // Set defaults output format
  if (!outputInit_){
    this->setOutputInfo("txt","chain.dat", max(1,(int) ncalls/100), max(1,(int) ncalls/20));
  }

  if(!epsMalaInit_){
    this->initEpsMALA(default_eps_mala_);
  }

  // Work variables for simplicity
  string output=outputinfo_.type;

  // Initial chain state
  currState_.step=0;
  currState_.state=chstart;
  currState_.alfa=0.0;
  currState_.post=this->evalLogPosterior(chstart);
  this->updateMode();
  fullChain_.PushBack(currState_);

  // Number of accepted steps and number of all trials
  int nacc=0;
  int nall=0;

  // No new mode found yet
  newMode_=false;

  for(int t = 1; t < nCalls, ++t){
    currState_.step=t;
    double sum_alpha=0.0;

    // Create a new proposed sample
    for(int is = 0; is  < nSubSteps_; ++is){
      Array1D<double> m_cand;
      this -> proposal(currState_.state, m_cand);

      // Evaluate the posterior at the new sample point
      double eval_cand = this->evalLogPosterior(m_cand);

      // Evaluate the new|old and old|new proposals
      double old_given_new = this->probOldNew(currState_.state, m_cand);
      double new_given_old = this->probOldNew(m_cand,currState_.state);

      // Accept or reject it
      double alpha = exp(eval_cand - currState_.post + old_given_new - new_given_old);
      sum_alpha += alpha;
      if (this->inDomain(m_cand) && (alpha>=1 || alpha > dsfmt_genrand_urv(&RandomState))) { // Accept and update the state
        nacc++;
        currState_.state = m_cand;
        currState_.post = eval_cand;
        if (this->fcnAcceptFlag_)
          this->fcnAccept_(this->postInfo_);
      } // If state not accepted, keep previous state as the current state
      else{
        if (this->fcnRejectFlag_)
          this->fcnReject_(this->postInfo_);
      }

      ++nall;
    }

    currState_.alfa = sum_alpha / nSubSteps_;

    // Append the current state to the array of all past states
    fullChain_.PushBack(currState_);

    // Keep track of the mode (among the locations visited so far)
    // \todo maybe only store tmode_(we save the full chain anyway)
    if (currState_.post > modeState_.post){
      this->updateMode();
      newMode_=true;
    }

    accRatio_ = (double) nacc/nall;

    if(WRITE_FLAG == 1){
      // Output to Screen
      if( t % outputinfo_.freq_outscreen == 0 || t==ncalls){

        printf("%lg %% completed; Chain step %d\n", 100.*t/ncalls,t);
        printf("================= Current logpost:%f, Max logpost:%f, Accept rate:%f\n",currState_.post,modeState_.post,accRatio_);
        printf("================= Current MAP params: ");
        for(int ic=0;ic<this -> GetChainDim();ic++)
          printf("par(%d)=%f ",ic,modeState_.state(ic));
        cout << endl;

      }

        // Output to File
      if( t % outputinfo_.freq_chainfile == 0 || t==ncalls){

        if(!strcmp(output.c_str(),"txt"))
          this->writeChainTxt(outputinfo_.filename);
        else  if(!strcmp(output.c_str(),"bin"))
          this->writeChainBin(outputinfo_.filename);
        else
          throw Tantrum((string) "Chain output type is not recognized");
        lastwrite_ = t;
      }
    }
  }

  return;
}

void SS::runChain(int ncalls, Array1D<double>& chstart){
  // Check mandatory information
  if(!(this -> getDimInit())){
    throw Tantrum((string) "Chain dimensionality needs to be initialized");
  }

  // Check what is not initialized and use defaults instead
  // \todo Specify defaults somewhere more transparently

  // Set defaults proposal covariance
  if(!propcovInit_){
    Array1D<double> chsig(this -> GetChainDim(),0.e0);
    for(int i=0;i<this -> GetChainDim();i++) chsig(i)=MAX(fabs(0.1*chstart(i)),0.001);
    this->initChainPropCovDiag(chsig);
  }

  // Set defaults output format
  if (!outputInit_){
    this->setOutputInfo("txt","chain.dat", max(1,(int) ncalls/100), max(1,(int) ncalls/20));
  }

  // Work variables for simplicity
  string output=outputinfo_.type;

  // Initial chain state
  currState_.step=0;
  currState_.state=chstart;
  currState_.alfa=0.0;
  currState_.post=this->evalLogPosterior(chstart);
  this->updateMode();
  fullChain_.PushBack(currState_);

  // Number of accepted steps and number of all trials
  int nacc=0;
  int nall=0;

  // No new mode found yet
  newMode_=false;

  for(int t = 1; t < nCalls, ++t){
    currState_.step=t;
    double sum_alpha=0.0;

    // Create a new proposed sample
    for(int is = 0; is  < nSubSteps_; ++is){
      Array1D<double> m_cand;
      this -> proposal(currState_.state, m_cand, is);

      // Evaluate the posterior at the new sample point
      double eval_cand = this->evalLogPosterior(m_cand);

      // Evaluate the new|old and old|new proposals
      double old_given_new = this->probOldNew(currState_.state, m_cand);
      double new_given_old = this->probOldNew(m_cand,currState_.state);

      // Accept or reject it
      double alpha = exp(eval_cand - currState_.post + old_given_new - new_given_old);
      sum_alpha += alpha;
      if (this->inDomain(m_cand) && (alpha>=1 || alpha > dsfmt_genrand_urv(&RandomState))) { // Accept and update the state
        nacc++;
        currState_.state = m_cand;
        currState_.post = eval_cand;
        if (this->fcnAcceptFlag_)
          this->fcnAccept_(this->postInfo_);
      } // If state not accepted, keep previous state as the current state
      else{
        if (this->fcnRejectFlag_)
          this->fcnReject_(this->postInfo_);
      }

      ++nall;
    }

    currState_.alfa = sum_alpha / nSubSteps_;

    // Append the current state to the array of all past states
    fullChain_.PushBack(currState_);

    // Keep track of the mode (among the locations visited so far)
    // \todo maybe only store tmode_(we save the full chain anyway)
    if (currState_.post > modeState_.post){
      this->updateMode();
      newMode_=true;
    }

    accRatio_ = (double) nacc/nall;

    if(WRITE_FLAG == 1){
      // Output to Screen
      if( t % outputinfo_.freq_outscreen == 0 || t==ncalls){

        printf("%lg %% completed; Chain step %d\n", 100.*t/ncalls,t);
        printf("================= Current logpost:%f, Max logpost:%f, Accept rate:%f\n",currState_.post,modeState_.post,accRatio_);
        printf("================= Current MAP params: ");
        for(int ic=0;ic<this -> GetChainDim();ic++)
          printf("par(%d)=%f ",ic,modeState_.state(ic));
        cout << endl;

      }

        // Output to File
      if( t % outputinfo_.freq_chainfile == 0 || t==ncalls){

        if(!strcmp(output.c_str(),"txt"))
          this->writeChainTxt(outputinfo_.filename);
        else  if(!strcmp(output.c_str(),"bin"))
          this->writeChainBin(outputinfo_.filename);
        else
          throw Tantrum((string) "Chain output type is not recognized");
        lastwrite_ = t;
      }
    }
  }

  return;
}

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
  if(!propcovInit_){
    Array1D<double> chsig(this -> GetChainDim(),0.e0);
    for(int i=0;i<this -> GetChainDim();i++) chsig(i)=MAX(fabs(0.1*chstart(i)),0.001);
    this->initChainPropCovDiag(chsig);
  }

  // Set defaults output format
  if (!outputInit_){
    this->setOutputInfo("txt","chain.dat", max(1,(int) ncalls/100), max(1,(int) ncalls/20));
  }

  // Work variables for simplicity
  string output=outputinfo_.type;

  double logevid;

  std::vector<double> samples;
  std::vector<double> logpriors;
  std::vector<double> logliks;
  std::ofstream evidFile("Evidence.dat");

  // Run TMCMC, get evidence
  logevid = tmcmc(samples, logpriors, logliks,
  tmcmc_gamma, ncalls,
  seed_, tmcmc_nprocs, this-> GetChainDim(),
  tmcmc_cv, tmcmc_MFactor,
  tmcmc_basis, tmcmc_CATSteps, WRITE_FLAG);

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
    currState_.step=i+1;
    for (int j=0; j<this -> GetChainDim(); j++ ){
      // std::cout << i << " " << j << " " << samples[i*this->chainDim_+j] << std::endl;
      sample(j) = samples[i*this -> GetChainDim()+j];
    }
    currState_.state=sample;
    currState_.alfa=0.0;
    currState_.post=logpriors[i]+logliks[i];
    fullChain_.PushBack(currState_);
  }

  evidFile << std::setprecision(18) << logevid << std::endl;

  if (WRITE_FLAG == 1){
    // Output to file
    if(!strcmp(output.c_str(),"txt"))
      this->writeChainTxt(outputinfo_.filename);
    else if(!strcmp(output.c_str(),"bin"))
      this->writeChainBin(outputinfo_.filename);
    else
      throw Tantrum((string) "Chain output type is not recognized");
  }

  evidFile.close();
}

void MCMC::runChain(int ncalls){
  Array1D<double> chstart(chainDim_,0.e0);

  this->runChain(ncalls, chstart);
}

bool MCMC::newModeFound(){
  return newMode_;
}

void MCMC::getAcceptRatio(double * accrat){
  *accrat = accRatio_;

  return;
}

double MCMC::getAcceptRatio(){
  return accRatio_;
}

int MCMC::GetChainDim() const{
  return chainDim_;
}

double MCMC::evalLogPosterior(Array1D<double>& m){
  return logPosterior_(m,postInfo_);
}

void MCMC::evalGradLogPosterior(Array1D<double>& m, Array1D<double>& grads){
  // Evaluate given the log-posterior function defined by the user in the constructor
  gradlogPosterior_(m,grads,postInfo_);

  return;
}

bool MCMC::inDomain(Array1D<double>& m){
  int nd = m.XSize();

  for (int id=0;id<nd;id++){
    if (!lower_flag_(id))
      return true;
    else if (m(id)<Lower_(id))
      return false;

    if (!upper_flag_(id))
      return true;
    else if (m(id)>Upper_(id))
      return false;
    }

  return true;
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
  else if(t > 1 && t < adaptstep(2)){
    for (int i = 0; i < this -> GetChainDim(); ++i) {
      curmean(i)  = ( curmean(i)*(t-1.) + m_t(i) )/t;
    }

    for(int i = 0; i < this -> GetChainDim(); ++i){
      for(int j = 0; j < i - 1; ++j){
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
    propLCov_=chcov;

    // Cholesky factorization of the proposal covariance propLCov_, done in-place
    // Note, for diagonal covariances, this is an overkill
    FTN_NAME(dpotrf)(&lu,&(this -> GetChainDim()), propLCov_.GetArrayPointer(),&(this -> GetChainDim()),&chol_info);
  }

  if ( ( t > adaptstep(0) ) && ( (t % adaptstep(1) ) ==  0 ) && t <= adaptstep(2) ){
    for (int i = 0; i < this -> GetChainDim(); ++i){
      for(int j = 0; j < this -> GetChainDim(); ++j){
        propLCov_(i,j) = sigma*(curcov(i,j) + (i==j) * eps_cov ) ;
      }
    }

    chcov = propLCov_;

    // Cholsky factorization of the proposal covariance propLCov_, done in-place
    FTN_NAME(dpotrf)(&lu,&(this -> GetChainDim()), propLCov_.GetArrayPointer(),&(this -> GetChainDim()),&chol_info);

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

void MALA::proposal(Array1D<double>& m_t,Array1D<double>& m_cand){
  Array1D<double> grads;
  gradlogPosterior_(m_t,grads,NULL);
  cout << "grads= " << grads(0) << " " << grads(1) << endl;
  m_cand = m_t;
  for (int i=0; i < this -> GetChainDim(); i++) {
    m_cand(i) += epsMALA_*epsMALA_*grads(i)/2.;
    m_cand(i) += epsMALA_*dsfmt_genrand_nrv(&RandomState);
  }

  return;
}

void MMALA::proposal(Array1D<double>& m_t,Array1D<double>& m_cand){
  int chol_info=0;
  char lu='L';

  Array1D<double> grads;
  this->gradlogPosterior_(m_t,grads,NULL);
  Array2D<double> mtensorinv;
  this->metricTensor_(m_t,mtensorinv,NULL);
  m_cand=m_t;

  Array1D<double> mtggrads;
  prodAlphaMatVec(mtensorinv, grads, 1.0, mtggrads) ;

  Array2D<double> sqrt_mtensorinv;
  sqrt_mtensorinv = mtensorinv;
  FTN_NAME(dpotrf)(&lu,&(this -> GetChainDim()), sqrt_mtensorinv.GetArrayPointer(),&(this -> GetChainDim()),&chol_info);
  // Catch the error in Cholesky factorization
  if (chol_info != 0 )
    printf("Error in Cholesky factorization, info=%d\n", chol_info);

  for (int i=0; i < this -> GetChainDim(); i++) {
    m_cand(i) += epsMALA_*epsMALA_*mtggrads(i)/2.;
    for (int j=0; j < i+1; j++) {
      m_cand(i) += epsMALA_*sqrt_mtensorinv(i,j)*dsfmt_genrand_nrv(&RandomState);
    }
  }

  return;
}

void SS::proposal(Array1D<double>& m_t,Array1D<double>& m_cand,int dim)
{
  // Single-site proposal
  m_cand = m_t;
  m_cand(dim) += ( sqrt(chcov(dim,dim))*dsfmt_genrand_nrv(&RandomState) );

  return;
}

double MALA::probOldNew(Array1D<double>& a, Array1D<double>& b){
  double logprob;
  Array1D<double> gradb;

  gradlogPosterior_(b,gradb,NULL);
  double eps2=this->epsMALA_*this->epsMALA_;
  Array1D<double> bmean(this->chainDim_,0.e0);
  Array1D<double> diagcov(this->chainDim_,0.e0);

  for (int i=0;i<this -> GetChainDim();i++){
    bmean(i)=b(i)+eps2*gradb(i)/2.0;
    diagcov(i)=eps2;
  }

  logprob=evallogMVN_diag(a,bmean,diagcov);

  return logprob;
}

double MMALA::probOldNew(Array1D<double>& a, Array1D<double>& b){
  ///\todo In the original code there is a seperate branch for MMALA in the function probOldNew, but it is blank and has nothing in it. It would need to be defined in this setup.
}

double MCMC::evallogMVN_diag(Array1D<double>& x,Array1D<double>& mu,Array1D<double>& sig2){

}

void MCMC::updateMode(){

}

void MCMC::writeChainTxt(string filename){

}

void MCMC::writeChainBin(string filename){

}
