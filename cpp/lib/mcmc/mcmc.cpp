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

MCMC::MCMC(double (*logposterior)(Array1D<double>&, void *), void *postinfo){
  // Set Flag
  FLAG = 0;

  // Setting entering the pointers to the Log Posterior function as well as the data function
  postInfo_ = postinfo;
  logPosterior_ = logposterior;

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

MCMC::MCMC(LogPosteriorBase& L){
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
  this -> setDefaultDomain();
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
  if (namesPrepend)
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

void MCMC::setCurrentStateState(Array1D<double>& newState){
  currState_.state = newState;
  return;
}

double MCMC::getCurrentStatePost(){
  return currState_.post;
}

void MCMC::setCurrentStatePost(double newPost){
  currState_.post = newPost;
  return;
}

double MCMC::getModeStatePost(){
  return modeState_.post;
}

void MCMC::setCurrentStateAlfa(double newAlfa){
  currState_.alfa = newAlfa;
  return;
}

void MCMC::getModeStateState(Array1D<double>& state){
  state = modeState_.state;
  return;
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

  lbfgsDR(n,m,start.GetArrayPointer(),nbd.GetArrayPointer(),l.GetArrayPointer(),u.GetArrayPointer(),neg_logposteriorproxy,NULL,info) ;

  this -> setCurrentStateStep(0);
  this -> setCurrentStateState(start);
  this -> setCurrentStateAlfa(0.0);
  this -> setCurrentStatePost(this->evalLogPosterior(start));
  this->updateMode();

  return;
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
  if (FLAG == 0){
    return logPosterior_(m,postInfo_);
  }
  if (FLAG == 1){
     return L_->eval(m);
  }
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

int MCMC::getLowerFlag(int i){
  return lower_flag_(i);
}

int MCMC::getUpperFlag(int i){
  return upper_flag_(i);
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
