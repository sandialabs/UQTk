/* =====================================================================================
                     The UQ Toolkit (UQTk) version 3.0.4
                     Copyright (2017) Sandia Corporation
                     http://www.sandia.gov/UQToolkit/

     Copyright (2017) Sandia Corporation. Under the terms of Contract DE-AC04-94AL85000
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
#include "mcmc.h"
#include "gen_defs.h"
#include "lbfgs_routines.h"


double neg_logposteriorproxy(int chaindim, double* m, void* classpointer);
void grad_neg_logposteriorproxy(int chaindim, double* m, double* grads, void* classpointer);

//******************NEW ROUTINES*************************//
void MCMC::setSeed(int seed){
  dsfmt_init_gen_rand(&RandomState,seed);
}
void MCMC::setWriteFlag(int I){
  WRITE_FLAG = I;
}

void MCMC::resetChainState(){
  fullChain_.Clear();
}
//******************NEW ROUTINES*************************//

MCMC::MCMC(double (*logPosterior)(Array1D<double>&, void *), void *postinfo)
{

  //*****************************************
  FLAG = 0;
  //*****************************************

  postInfo_=postinfo;
  logPosterior_=logPosterior;

  // / Set the initialization flags to false
  gradflag_=false;
  tensflag_=false;

  chaindimInit_=false;
  propcovInit_=false;
  methodInit_=false;
  outputInit_=false;
  adaptstepInit_=false;
  gammaInit_=false;
  epscovInit_=false;
  epsMalaInit_=false;


  // Initiate the random number generator seed
  // \todo This needs to be made more generic
  seed_=13;
  dsfmt_init_gen_rand(&RandomState,seed_);

  // Set the location of the last chain state written (-1 means nothing is written to files yet)
  lastwrite_=-1;

  // By default the names are not prepended
  namePrepend_=false;

  // Set defaults
  this->initDefaults();

  WRITE_FLAG = 1;

  return;
}

//*****************************************
MCMC::MCMC(LikelihoodBase& L)
{
  FLAG = 1;
  L_ = &L;

  // / Set the initialization flags to false
  gradflag_=false;
  tensflag_=false;

  chaindimInit_=false;
  propcovInit_=false;
  methodInit_=false;
  outputInit_=false;
  adaptstepInit_=false;
  gammaInit_=false;
  epscovInit_=false;
  epsMalaInit_=false;

  // Set the location of the last chain state written (-1 means nothing is written to files yet)
  lastwrite_=-1;

  // Initiate the random number generator seed
  // \todo This needs to be made more generic
  seed_=13;
  dsfmt_init_gen_rand(&RandomState,seed_); // you can override with setseed

  // By default the names are not prepended
  namePrepend_=false;

  // Set defaults
  this->initDefaults();

  WRITE_FLAG = 1;


  return;
}

// create samples where each column is a sample
void MCMC::getSamples(int burnin, int every,Array2D<double>& samples)
{
    int nCalls = fullChain_.Length();
    samples.Resize(chainDim_,0); // initialize sample array
    int j=0;
    for (int i = burnin; i < nCalls; i+=every){
        samples.insertCol(fullChain_(i).state,j);
        j++;
    }
}

// create samples where each column is a sample
void MCMC::getSamples(Array2D<double>& samples)
{
    getSamples(0,1,samples);
}
//*****************************************


void MCMC::setGradient(void (*gradlogPosterior)(Array1D<double>&, Array1D<double>&, void *))
{
  gradlogPosterior_ = gradlogPosterior;
  gradflag_ = true;
  return;
}

void MCMC::setMetricTensor(void (*metricTensor)(Array1D<double>&, Array2D<double>&, void *))
{
  metricTensor_ = metricTensor;
  tensflag_ = true;
  return;
}

void MCMC::initDefaults()
{
  this->default_method_="am";
  this->default_gamma_=0.01;
  this->default_eps_cov_=1e-8;

  this->default_eps_mala_=0.1;

  this->newMode_=false;
  this->accRatio_ = -1.0;

  return;
}

void MCMC::printChainSetup()
{
  if (this->methodInit_)
    cout << "Method           : " << this->methodinfo_.type << endl;
  else
    cout << "Method (default) : " << this->default_method_ << endl;
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



double MCMC::evalLogPosterior(Array1D<double>& m){
  // Evaluate given the log-posterior function defined by the user in the constructor
  //*****************************************
  if (FLAG == 0){
    return logPosterior_(m,postInfo_);
  }
  if (FLAG == 1){
     return L_->eval(m);
  }

  //*****************************************
}

void MCMC::evalGradLogPosterior(Array1D<double>& m, Array1D<double>& grads)
{
  // Evaluate given the log-posterior function defined by the user in the constructor
  gradlogPosterior_(m,grads,postInfo_);

  return;

}

void MCMC::initMethod(string method)
{
  // Set the method type
  methodinfo_.type=method;
  // Set the initialization flag to True
  methodInit_=true;
  return;
}

void MCMC::initAdaptSteps(int adaptstart,int adaptstep, int adaptend)
{
  // Initialize the vector containing the adaptivity information: when to start, how often and when to stop.
  methodinfo_.adaptstep.Resize(3);
  methodinfo_.adaptstep(0)=adaptstart;
  methodinfo_.adaptstep(1)=adaptstep;
  methodinfo_.adaptstep(2)=adaptend;

  // Set the initialization flag to True
  adaptstepInit_=true;

  return;
}

void MCMC::initAMGamma(double gamma)
{
  // Initialize the scale factor gamma
  methodinfo_.gamma=gamma;
  // Set the initialization flag to True
  gammaInit_=true;

  return;
}

void MCMC::initEpsCov(double eps_cov)
{
  // Initialize the covariance 'nugget'
  methodinfo_.eps_cov=eps_cov;
  // Set the initialization flag to True
  epscovInit_=true;

  return;
}

void MCMC::initEpsMALA(double eps_mala)
{
  // Initialize the epsilon parameter for MALA algorithm
  epsMALA_=eps_mala;
  // Set the initialization flag to True
  epsMalaInit_=true;

  return;
}

void MCMC::setOutputInfo(string outtype, string file,int freq_file, int freq_screen)
{
  outputinfo_.type=outtype;
  outputinfo_.filename=file;
  outputinfo_.freq_chainfile=freq_file;
  outputinfo_.freq_outscreen=freq_screen;
  // Set the initialization flag to True
  outputInit_=true;
  return;
}

void MCMC::initChainPropCov(Array2D<double>& propcov)
{
  // Initialize the proposal covariance matrix
  methodinfo_.chcov=propcov;
  // Set the initialization flag to True
  propcovInit_=true;
  return;
}

void MCMC::initChainPropCovDiag(Array1D<double>& sig)
{

  // Create a diagonal matrix and fill in the diagonal terms
  methodinfo_.chcov.Resize(this->chainDim_,this->chainDim_,0.e0);
  for(int i=0;i<this->chainDim_;i++) methodinfo_.chcov(i,i)=sig(i)*sig(i);
  // Set the initialization flag to True
  propcovInit_=true;
  return;
}

void MCMC::getChainPropCov(Array2D<double>& propcov)
{
  // Get the proposal covariance matrix
  propcov=methodinfo_.chcov;
  return;
}

double neg_logposteriorproxy(int chaindim, double* m, void* classpointer)
{

  MCMC* thisClass=(MCMC*) classpointer;

  int aa=thisClass->GetChainDim();

// Double check chain dimensionality
  if(chaindim != thisClass->GetChainDim()){

    throw Tantrum(std::string("neg_logposteriorproxy: The passed in MCMC chain dimension") +
                  " does not match the dimension of the MChain class instance");
}

  Array1D<double> mm(chaindim,0.e0);

  for(int i=0;i<chaindim;i++)
    mm(i)=m[i];




// Call the posterior function and return its result
  return -thisClass->evalLogPosterior(mm);




}

void grad_neg_logposteriorproxy(int chaindim, double* m, double* grads, void* classpointer)
{

  MCMC* thisClass=(MCMC*) classpointer;

  int aa=thisClass->GetChainDim();

// Double check chain dimensionality
  if(chaindim != thisClass->GetChainDim()){

    throw Tantrum(std::string("grad_neg_logposteriorproxy: The passed in MCMC chain dimension") +
                  " does not match the dimension of the MChain class instance");
}

  Array1D<double> mm(chaindim,0.e0);

  for(int i=0;i<chaindim;i++)
    mm(i)=m[i];




// Call the posterior function and return its result
  Array1D<double> grads_arr;
  thisClass->evalGradLogPosterior(mm, grads_arr);


  for(int i=0;i<chaindim;i++)
    grads[i]=-grads_arr(i);


  return;


}

void MCMC::runOptim(Array1D<double>& start)
{
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



void MCMC::runChain(int ncalls, Array1D<double>& chstart)
{
  // Check the mandatory initialization
  if(!chaindimInit_)
    throw Tantrum((string) "Chain dimensionality needs to be initialized");

  // Check what is not initialized and use defaults instead
  // \todo Specify defaults somewhere more transparently

  // Set defaults proposal covariance
  if(!propcovInit_){
    Array1D<double> chsig(this->chainDim_,0.e0);
    for(int i=0;i<this->chainDim_;i++) chsig(i)=MAX(fabs(0.1*chstart(i)),0.001);
    this->initChainPropCovDiag(chsig);
  }

  // Set defaults output format
  if (!outputInit_)
    this->setOutputInfo("txt","chain.dat", max(1,(int) ncalls/100), max(1,(int) ncalls/20));


  // Set the default method
  if (!methodInit_)
    methodinfo_.type=default_method_;

  // Set the default parameters for aMCMC
  if(!strcmp(this->methodinfo_.type.c_str(),"am")){

    if(!adaptstepInit_)
      this->initAdaptSteps((int) ncalls/10,10,ncalls);

    if(!gammaInit_)
      this->initAMGamma(this->default_gamma_);

    if(!epscovInit_)
      this->initEpsCov(this->default_eps_cov_);

  }
  else if(!strcmp(this->methodinfo_.type.c_str(),"mala")){
    if(!epsMalaInit_)
      this->initEpsMALA(this->default_eps_mala_);
  }


  // For simplicity, work variables
  string method=methodinfo_.type;
  string output=outputinfo_.type;


  // Set the number of substeps per one chain step
  if(!strcmp(method.c_str(),"ss"))
    nSubSteps_=chainDim_;
  else //if(!strcmp(method.c_str(),"am")) or if(!strcmp(method.c_str(),"mala"))
    nSubSteps_=1;


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

  // Main loop
  for (int t=1; t <= ncalls; t++) {
    currState_.step=t;
    double sum_alpha=0.0;

    // Create a new proposed sample
    // If the method is 'ss'(Single-Site), one actually searches all dimensions before recording the new state, i.e.
    // one chain step has d substeps, where d is the chain dimensionality.
    // For 'am'(Adaptive), d is set to 1.
    for (int is=0; is<nSubSteps_; is++){
      Array1D<double> m_cand;
      if(!strcmp(method.c_str(),"ss"))
        this->proposalSingleSite(currState_.state, m_cand, is);
      else if(!strcmp(method.c_str(),"am"))
        this->proposalAdaptive(currState_.state, m_cand, t);
      else if(!strcmp(method.c_str(),"mala"))
        this->proposalMALA(currState_.state, m_cand);
      else if(!strcmp(method.c_str(),"mmala"))
        this->proposalMMALA(currState_.state, m_cand);
      else
        throw Tantrum((string) "Chain running method is not recognized");

      // Evaluate the posterior at the new sample point
      double eval_cand = this->evalLogPosterior(m_cand);

      // Evaluate the new|old and old|new proposals
      double old_given_new = this->probOldNew(currState_.state, m_cand);
      double new_given_old = this->probOldNew(m_cand,currState_.state);



      // Accept or reject it
      double alpha = exp(eval_cand - currState_.post + old_given_new - new_given_old);
        // cout << t << " " << eval_cand << " " << currState_.post << " " << old_given_new << " " << new_given_old << " " << alpha << endl;
       // cout << "alpha = " << alpha << endl;
      sum_alpha+=alpha;
      if (alpha>=1 || alpha > dsfmt_genrand_urv(&RandomState)){ // Accept and update the state
          if (this->inDomain(m_cand))
          {
              nacc++;
              currState_.state = m_cand;
              currState_.post = eval_cand;
          }
      } // If state not accepted, keep previous state as the current state

      nall++;
    }
    currState_.alfa=sum_alpha/nSubSteps_;

    // Append the current state to the array of all past states
    fullChain_.PushBack(currState_);

    // Keep track of the mode (among the locations visited so far)
    // \todo maybe only store tmode_(we save the full chain anyway)
    if (currState_.post > modeState_.post){
      this->updateMode();
      newMode_=true;
    }

    accRatio_ = (double) nacc/nall;



    if (WRITE_FLAG == 1){

      // Output to the screen
      if( t % outputinfo_.freq_outscreen == 0 || t==ncalls){

        printf("%lg %% completed; Chain step %d\n", 100.*t/ncalls,t);
        printf("================= Current logpost:%f, Max logpost:%f, Accept rate:%f\n",currState_.post,modeState_.post,accRatio_);
        printf("================= Current MAP params: ");
        for(int ic=0;ic<this->chainDim_;ic++)
          printf("par(%d)=%f ",ic,modeState_.state(ic));
        cout << endl;

      }

      // Output to file
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

  }  // End of main loop

  return;
}

void MCMC::updateMode()
{
  // Update the chain mode (MAP state)
  modeState_.step=currState_.step;
  modeState_.state=currState_.state;
  modeState_.post=currState_.post;
  modeState_.alfa=-1.0;

  return;
}

bool MCMC::newModeFound()
{
  // Check to see if a new mode was found during last call to runChain
  return newMode_;
}

void MCMC::getAcceptRatio(double * accrat)
{
  // Returns the acceptance ratio
  *accrat = accRatio_;
  return;
}

void MCMC::proposalSingleSite(Array1D<double>& m_t,Array1D<double>& m_cand,int dim)
{
  // Single-site proposal
  m_cand=m_t;
  m_cand(dim) += ( sqrt(methodinfo_.chcov(dim,dim))*dsfmt_genrand_nrv(&RandomState) );

  return;
}


void MCMC::proposalAdaptive(Array1D<double>& m_t,Array1D<double>& m_cand,int t)
{
  int chol_info=0;
  char lu='L';

  // xm[] is the mean of m_t[] over all previous states, X_0,...,X_{t-1}
  // at this stage, index t, we know X_0,X_1,...,X_{t-1}
  // and we're seeking to find X_t, the new state of the chain
  // also evaluate covt, the covariance matrix

  if (t == 1) {
 // at the first iteration, the mean is easy and the sample covariance is 0
     methodinfo_.curmean=m_t;
     methodinfo_.curcov.Resize(this->chainDim_,this->chainDim_,0.e0);
  } else if( t > 1 && t < methodinfo_.adaptstep(2) ){
// after the first iteration, start keeping track of the sample mean
     for (int i=0; i < chainDim_; i++) {
          methodinfo_.curmean(i)  = ( methodinfo_.curmean(i)*(t-1.) + m_t(i) )/t;
     }



     for (int i=0; i < this->chainDim_; i++)
       for (int j=0; j < i+1; j++)
   methodinfo_.curcov(i,j) = ( (t-2.)/(t-1.) )*methodinfo_.curcov(i,j) +
     ( t/((t-1.)*(t-1.)) ) * ( m_t(i) - methodinfo_.curmean(i) )*( m_t(j) - methodinfo_.curmean(j) );

     //transpose
     for (int i=0; i < chainDim_; i++)
       for (int j=i+1; j < chainDim_ ; j++)
   methodinfo_.curcov(i,j) = methodinfo_.curcov(j,i) ;

  }



  // Jump size
  double sigma = methodinfo_.gamma * 2.4 * 2.4 / (double)this->chainDim_;

  if(t ==1) {
    propLCov_=methodinfo_.chcov;

    // Cholesky factorization of the proposal covariance propLCov_, done in-place
    // Note, for diagonal covariances, this is an overkill
    FTN_NAME(dpotrf)(&lu,&chainDim_, propLCov_.GetArrayPointer(),&chainDim_,&chol_info);

  }

  if ( ( t > methodinfo_.adaptstep(0) ) && ( (t % methodinfo_.adaptstep(1) ) ==  0 ) && t <= methodinfo_.adaptstep(2) ) {

    for (int i=0; i < chainDim_; i++)
      for (int j=0; j < chainDim_; j++)
  propLCov_(i,j) = sigma*( methodinfo_.curcov(i,j) + (i==j)*methodinfo_.eps_cov ) ;

    methodinfo_.chcov=propLCov_;


    // Cholsky factorization of the proposal covariance propLCov_, done in-place
    FTN_NAME(dpotrf)(&lu,&chainDim_, propLCov_.GetArrayPointer(),&chainDim_,&chol_info);

    // Catch the error in Cholesky factorization
    if (chol_info != 0 ) {
      printf("Error in Cholesky factorization, info=%d, printing the matrix below:\n", chol_info);

      for(int i=0;i<chainDim_;i++){
        for(int j=0;j<chainDim_;j++)
          printf("%lg ",propLCov_(i,j));
        printf("\n");
      }

      exit(1);
    }
  }


  // Candidate state is a multivariate normal sample away from the current state
    m_cand=m_t;
    Array1D<double> xi(chainDim_,0.e0);
  for (int i=0; i < chainDim_; i++) {
      xi(i)=dsfmt_genrand_nrv(&RandomState);
    double Lnrv=0.0;
    for (int j=0; j < i+1; j++) {
      Lnrv += propLCov_(i,j)*xi(j);
    }
    m_cand(i) += Lnrv;
  }



  return;

}

void MCMC::proposalMALA(Array1D<double>& m_t,Array1D<double>& m_cand)
{
  Array1D<double> grads;
  gradlogPosterior_(m_t,grads,NULL);
  cout << "grads= " << grads(0) << " " << grads(1) << endl;
  m_cand=m_t;
  for (int i=0; i < chainDim_; i++) {
    m_cand(i) += epsMALA_*epsMALA_*grads(i)/2.;
    m_cand(i) += epsMALA_*dsfmt_genrand_nrv(&RandomState);
  }


  return;
}


void MCMC::proposalMMALA(Array1D<double>& m_t,Array1D<double>& m_cand)
{

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
  sqrt_mtensorinv=mtensorinv;
  FTN_NAME(dpotrf)(&lu,&chainDim_, sqrt_mtensorinv.GetArrayPointer(),&chainDim_,&chol_info);
  // Catch the error in Cholesky factorization
  if (chol_info != 0 )
      printf("Error in Cholesky factorization, info=%d\n", chol_info);


  for (int i=0; i < chainDim_; i++) {
    m_cand(i) += epsMALA_*epsMALA_*mtggrads(i)/2.;
    for (int j=0; j < i+1; j++) {
      m_cand(i) += epsMALA_*sqrt_mtensorinv(i,j)*dsfmt_genrand_nrv(&RandomState);
    }
  }


  return;
}

void MCMC::appendMAP()
{
  this->fullChain_.PushBack(modeState_);
  return;
}

void MCMC::writeChainTxt(string filename)
{

  // Choose whether write or append
  char* writemode="w";
  if (lastwrite_>=0 || namePrepend_)
    writemode="a";

  // Open the text file
  FILE* f_out;
  if(!(f_out = fopen(filename.c_str(),writemode))){
    printf("writeChain: could not open file '%s'\n",filename.c_str());
    exit(1);
  }

  // Write to the text file
  for(int i=this->lastwrite_+1;i<this->fullChain_.XSize();i++){
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

  return ;

}

void MCMC::writeFullChainTxt(string filename, Array1D<chainstate> fullchain)
{

  // Open the text file in a write mode
  char* writemode="w";

  // Append if the names already prepended
  if (namePrepend_)
    writemode="a";

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

 return ;

}


void MCMC::writeChainBin(string filename)
{

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

  return ;

}

void MCMC::parseBinChain(string filename, Array1D<chainstate>& readchain)
{
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


double MCMC::getMode(Array1D<double>& MAPparams)
{


  //for(int ic=0;ic<this->chainDim_;ic++)
    MAPparams=modeState_.state;

  return modeState_.post;
}


double MCMC::probOldNew(Array1D<double>& a, Array1D<double>& b)
{
 string method=methodinfo_.type;

  double logprob;
  Array1D<double> gradb;

  if(!strcmp(method.c_str(),"mala")){
    gradlogPosterior_(b,gradb,NULL);
    double eps2=this->epsMALA_*this->epsMALA_;
    Array1D<double> bmean(this->chainDim_,0.e0);
    Array1D<double> diagcov(this->chainDim_,0.e0);

    for (int i=0;i<chainDim_;i++){
      bmean(i)=b(i)+eps2*gradb(i)/2.0;
      diagcov(i)=eps2;
    }

    logprob=evallogMVN_diag(a,bmean,diagcov);
  }
  else if(!strcmp(method.c_str(),"mmala")) {
    //logprob=
  }
 else
   return 0.0;

  return logprob;
}

double MCMC::evallogMVN_diag(Array1D<double>& x,Array1D<double>& mu,Array1D<double>& sig2)
{
  double pi=4.0*atan(1.0);

  double value=0.e0;

  // \todo Put sanity checks on dimensions

  for (int i=0;i<this->chainDim_;i++){
    value -= 0.5*log(2.*pi*sig2(i));
    value -= (x(i)-mu(i))*(x(i)-mu(i))/(2.0*sig2(i));
  }
  return value;
}


bool MCMC::inDomain(Array1D<double>& m)
{
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

void MCMC::setLower(double lower, int i)
{
        this->Lower_(i)=lower;
        lower_flag_(i)=1;

    return;
}

void MCMC::setUpper(double upper, int i)
{
    this->Upper_(i)=upper;
    upper_flag_(i)=1;
    return;
}


void MCMC::setChainDim(int chdim)
{
    this->chainDim_=chdim; chaindimInit_=true;
    this->setDefaultDomain();

    return;
}

void MCMC::setDefaultDomain()
{
    this->Lower_.Resize(this->chainDim_,-DBL_MAX);
    this->Upper_.Resize(this->chainDim_,DBL_MAX);
    this->lower_flag_.Resize(this->chainDim_,0);
    this->upper_flag_.Resize(this->chainDim_,0);


    return;
}
