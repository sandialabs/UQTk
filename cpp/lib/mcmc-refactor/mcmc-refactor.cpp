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

void MCMC::initChainPropCovDiag(Array1D<double>& sig){
  // Create a diagonal matrix and fill in the diagonal terms
  methodinfo_.chcov.Resize(this->chainDim_,this->chainDim_,0.e0);
  for(int i = 0; i < chainDim_; ++i){
    methodinfo_.chcov(i,i)=sig(i)*sig(i);
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
  Lower_.Resize(this->chainDim_,-DBL_MAX);
  Upper_.Resize(this->chainDim_,DBL_MAX);
  lower_flag_.Resize(this->chainDim_,0);
  upper_flag_.Resize(this->chainDim_,0);
  
  return;
}

void MCMC::getChainPropCov(Array2D<double>& propcov){
  // Get the proposal covariance matrix
  propcov=methodinfo_.chcov;
  
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

void AMCMC::printChainSetup(){
  if (this->gammaInit_)
    cout << "Gamma            : " << this->methodinfo_.gamma << endl;
  else
    cout << "Gamma (default)  : " << this->default_gamma_ << endl;
  if (this->epscovInit_)
    cout << "Eps_Cov          : " << this->methodinfo_.eps_cov << endl;
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

    curchain.state.Resize(this->chainDim_);

    fread(curchain.state.GetArrayPointer(),this->chainDim_*sizeof(double),1,fb);
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
    for(int ic=0;ic<this->chainDim_;ic++)
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
  if(!chaindimInit_){
    throw Tantrum((string) "Chain dimensionality needs to be initialized");
  }

  // Check what is not initialized and use defaults instead
  // \todo Specify defaults somewhere more transparently

  // Set defaults proposal covariance
  if(!propcovInit_){
    Array1D<double> chsig(this->chainDim_,0.e0);
    for(int i=0;i<this->chainDim_;i++) chsig(i)=MAX(fabs(0.1*chstart(i)),0.001);
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
      proposal(currState_.state, m_cand, t);

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
        for(int ic=0;ic<this->chainDim_;ic++)
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