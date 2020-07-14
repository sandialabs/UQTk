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
#include "mcmc.h"
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

void MALA::setGradient(void (*gradlogPosterior)(Array1D<double>&, Array1D<double>&, void *)){
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
  // Set the initialization flag to true
  propcovInit_=true;

  return;
}

void MCMC::initChainPropCovDiag(Array1D<double>& sig){
  // Create a diagonal matrix and fill in the diagonal terms
  chcov.Resize(chainDim_,chainDim_,0.e0);
  for(int i = 0; i < chainDim_; ++i){
    chcov(i,i)=sig(i)*sig(i);
  }

  // Set the initialization flag to true
  propcovInit_=true;

  return;
}

void MCMC::setOutputInfo(string outtype, string file,int freq_file, int freq_screen){
  outputinfo_.outtype = outtype;
  outputinfo_.filename = file;
  outputinfo_.freq_file = freq_file;
  outputinfo_.freq_screen = freq_screen;
  // Set the initialization flag to true
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
  Upper_(i) = upper;
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

void MALA::getGradient(void (*gradlogPosterior)(Array1D<double>&, Array1D<double>&, void *)){
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

void MCMC::getFcnReject(void (*fcnReject)(void *)){
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
  return Lower_(i);
}

double MCMC::getUpper(int i){
  return Upper_(i);
}

bool MCMC::getDimInit(){
  return chaindimInit_;
}

bool MALA::getGradientFlag(){
  return gradflag_;
}

void MCMC::getPropLCov(Array2D<double>& LCov){
  LCov = propLCov_;
  return;
}

void MCMC::getPostInfo(void *post){
  post = postInfo_;
  return;
}

bool MCMC::getPropCovInit(){
  return propcovInit_;
}

bool MCMC::getOutputInit(){
  return outputInit_;
}

int MCMC::getLastWrite(){
  return lastwrite_;
}

void MCMC::setLastWrite(int i){
  lastwrite_ = i;
  return;
}

void MCMC::setAcceptRatio(double d){
  accRatio_ = d;
  return;
}

void MCMC::runAcceptFcn(){
  this->fcnAccept_(this->postInfo_);
  return;
}

void MCMC::runRejectFcn(){
  this->fcnReject_(this->postInfo_);
  return;
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

int MCMC::getFullChainSize(){
  return this -> fullChain_.XSize();
}

void MCMC::addCurrentState(){
  fullChain_.PushBack(currState_);
}

void MCMC::setCurrentStateStep(int i){
  currState_.step = i;
  return;
}

void MCMC::getCurrentStateState(Array1D<double>& state){
  state = currState_.state;
  return;
}

double MCMC::getCurrentStatePost(){
  return currState_.post;
}

double MCMC::getModeStatePost(){
  return modeState_.post;
}

void MCMC::getModeStateState(Array1D<double>& state){
  state = modeState_.state;
  return;
}

/// \todo This might not be a MALA exclusive function but no other method needs a gradient
void MALA::runOptim(Array1D<double>& start){
  int n=start.Length();
  int m=5;
  Array1D<int> nbd(n,0);
  Array1D<double> l(n,0.e0);
  Array1D<double> u(n,0.e0);

  for (int i=0;i<n;i++){
    if (lower_flag_(i) && !upper_flag_(i)){
      nbd(i)=1;
      l(i)=this->getLower(i);
    }
    if (upper_flag_(i) && !lower_flag_(i)){
      nbd(i)=3;
      u(i)=this->getUpper(i);
    }
    if (upper_flag_(i) && lower_flag_(i)){
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

  for(int t = 1; t < ncalls; ++t){
    this -> setCurrentStateStep(t);
    double sum_alpha=0.0;

    // Create a new proposed sample
    for(int is = 0; is  < nSubSteps_; ++is){
      Array1D<double> m_cand;
      Array1D<double> state;
      this -> getCurrentStateState(state);
      this -> proposal(state, m_cand, t);

      // Evaluate the posterior at the new sample point
      double eval_cand = this->evalLogPosterior(m_cand);

      // Evaluate the new|old and old|new proposals
      /// \Note If there is an issue with the code check here first because if proposal changes the current state then the state object might change
      double old_given_new = this->probOldNew(state, m_cand);
      double new_given_old = this->probOldNew(m_cand,state);

      // Accept or reject it
      double post = this -> getCurrentStatePost();
      double alpha = exp(eval_cand - post + old_given_new - new_given_old);
      sum_alpha += alpha;
      if (this->inDomain(m_cand) && (alpha>=1 || alpha > dsfmt_genrand_urv(&RandomState))) { // Accept and update the state
        nacc++;
        /*currState_.state = m_cand;
        currState_.post = eval_cand;*/
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

void MALA::evalGradLogPosterior(Array1D<double>& m, Array1D<double>& grads){
  // Evaluate given the log-posterior function defined by the user in the constructor
  void *post;
  this -> getPostInfo(post);
  gradlogPosterior_(m,grads,post);

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
  else if(t > 1 && t < adaptstep_(2)){
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
  int chdim = this -> GetChainDim();
  FTN_NAME(dpotrf)(&lu,&chdim, sqrt_mtensorinv.GetArrayPointer(),&chdim,&chol_info);
  // Catch the error in Cholesky factorization
  if (chol_info != 0 )
    printf("Error in Cholesky factorization, info=%d\n", chol_info);

  for (int i=0; i < this -> GetChainDim(); i++) {
    m_cand(i) += this -> getEpsMALA() * this -> getEpsMALA() * mtggrads(i)/2.;
    for (int j=0; j < i+1; j++) {
      m_cand(i) += this -> getEpsMALA() *sqrt_mtensorinv(i,j)*dsfmt_genrand_nrv(&RandomState);
    }
  }

  return;
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

double MMALA::probOldNew(Array1D<double>& a, Array1D<double>& b){
  return 0.0;
  ///\todo In the original code there is a seperate branch for MMALA in the function probOldNew, but it is blank and has nothing in it. It would need to be defined in this setup.
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

void MCMC::updateMode(){
  // Update the chain mode (MAP state)
  modeState_.step=currState_.step;
  modeState_.state=currState_.state;
  modeState_.post=currState_.post;
  modeState_.alfa=-1.0;

  return;
}

void MCMC::setNewMode(bool value){
  newMode_ = value;
  return;
}

bool MCMC::getFcnAcceptInit(){
  return fcnAcceptFlag_;
}

bool MCMC::getFcnRejectInit(){
  return fcnRejectFlag_;
}

void MCMC::writeChainTxt(string filename){
  // Choose whether write or append
  char* writemode="w";
  if (this -> getLastWrite() >= 0 || namesPrepend)
    writemode="a";

  // Open the text file
  FILE* f_out;
  if(!(f_out = fopen(filename.c_str(),writemode))){
    printf("writeChain: could not open file '%s'\n",filename.c_str());
    exit(1);
  }

  // Write to the text file
  for(int i = this -> getLastWrite() + 1; i<this->fullChain_.XSize();i++){
    fprintf(f_out, "%d ", this->fullChain_(i).step);
    for(int ic=0;ic<this->chainDim_;ic++)
      fprintf(f_out, "%24.16lg ", this->fullChain_(i).state(ic));
    fprintf(f_out, "%24.16lg %24.16lg \n", this->fullChain_(i).alfa, this->fullChain_(i).post);

  }

  // Close the text file
  if(fclose(f_out)){
    printf("writeChain: could not close file '%s'\n",filename.c_str());
    exit(1);
  }

  // Report
  printf("Written the states %d - %d to the text file %s\n", this->lastwrite_+1,  this->fullChain_.XSize()-1, filename.c_str());

  return;
}

void MCMC::writeChainBin(string filename){
  // Choose whether write or append
  char* writemode="wb";
  if (lastwrite_>=0)
    writemode="ab";

  // Open the binary file
  FILE* f_out;
  if(!(f_out = fopen(filename.c_str(),writemode))){
    printf("writeChain: could not open file '%s'\n",filename.c_str());
    exit(1);
  }


  // Write to the binary file
  for(int i=this->lastwrite_+1;i<this->fullChain_.XSize();i++){
    fwrite(&(this->fullChain_(i).step), sizeof(int), 1, f_out);
    fwrite(fullChain_(i).state.GetArrayPointer(),this->chainDim_*sizeof(double),1, f_out);
    fwrite(&(this->fullChain_(i).alfa), sizeof(double), 1, f_out);
    fwrite(&(this->fullChain_(i).post), sizeof(double), 1, f_out);
  }

  // CLose the binary file
  if(fclose(f_out)){
    printf("writeChain: could not close file '%s'\n",filename.c_str());
    exit(1);
  }

  // Report
  printf("Written the states %d - %d to the binary file %s\n",this->lastwrite_+1, this->fullChain_.XSize()-1, filename.c_str());

  return;
}

void MALA::initEpsMALA(double eps_mala_){
  this -> eps_mala = eps_mala_;
  epsMalaInit_ = true;

  return;
}

double MALA::getEpsMALA(){
  return this -> eps_mala;
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

int SS::getNSubSteps(){
  return this -> nSubSteps_;
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
