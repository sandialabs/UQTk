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
//#include "tmcmc.h"
#include "gen_defs.h"
#include "lbfgs_routines.h"
#include "tmcmcClass.h"

#define BETA_MAX 0.3

extern "C" {
 void dpotrf_(char *, int *, double *, int*, int *);
}

bool fileExists(std::string fname);
void outProcToFile(const RealVector spls, const int ndim, const int nspl,
                    int nprocs) ;
void outProcToFile(const RealVector spls, const int ndim, const int nspl,
                    std::string fname);
void shuffle_spls(RealVector &spls, RealVector &llik, RealVector &lprior);

void readInitSamples(RealVector &spls, std::string fname);
void parseSetup(dsfmt_t &RandomState, int nspl, int iseed);
void PriorGen(dsfmt_t &RandomState, int nspl, CharVector &distr,
  RealVector &means, RealVector &vars, int iseed);
double pearsonCorrCoef(RealVector X, RealVector Y);
double rescaleSTD(RealVector w, double wmean);

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

double tmcmc(RealVector &spls, RealVector &lprior, RealVector &llik,
          double gm, int nspl,
          int iseed, int nProcs, int ndim, double cv,
          int MFactor, bool basis, int CATSteps, int write_flag) {
  /* TMCMC Algorithm
      Input: spls - On return, contains samples according to posterior
             gm - initial gamma value
             nspl - number of Samples
             iseed - random seed
             nProcs - number of processors
             ndim - dimensionality
             cv - Coefficient of Variance threshold for adapting Beta
             MFactor - Multiplicative factor for chain length to encourage mixing.
      Output: evid - asymptotically unbiased model evidence estimator
  */

  int nscal = ndim*nspl;


  // Roberts and Rosenthal 2011 - Initial gamma if one is not provided.
  if (gm < 0) gm = 2.38 / sqrt(ndim);
  double gm2 = gm * gm;

  std::ofstream gamma_file("gamma.dat");

  // std::cout<<"----------------------------------------------------"<<std::endl;
  // std::cout<<"No. of samples    :"<<nspl<<std::endl;
  // std::cout<<"Dimensionality    :"<<ndim<<std::endl;
  // std::cout<<"BASIS/CATMIPs     :"<<basis<<std::endl;
  // if (basis)
  // std::cout<<"CATSteps          :"<<CATSteps<<std::endl;
  // std::cout<<"MFactor           :"<<MFactor<<std::endl;
  // std::cout<<"----------------------------------------------------"<<std::endl;

  /* initialize random number generator */
  dsfmt_t RandomState;
  dsfmt_init_gen_rand(&RandomState,iseed);

  /* Generate Initial Random Samples */
  /* Check if an initial sample file exists, otherwise
    generate from a d-dim hypercube */
  spls.clear();
  readInitSamples(spls, "tmcmc_prior_samples.dat");
  // Catch the error in reading prior samples
  if (spls.size() != nscal) {
    printf("Number of prior samples read from file not equal to number of samples requested\n");
    exit(1);
  }

  assert(spls.size() == nscal);
  if (write_flag == 1){
    std::cout << "Samples read from prior_samples.dat\n";
  }


  /* Compute initial likelihoods using outside model*/
  outProcToFile(spls, ndim, nspl, nProcs); // Make mcmcstates
  std::string ll_stream="./tmcmc_getLL.sh "+std::to_string(nProcs);
  system(ll_stream.c_str());
  std::ifstream input_file("tmcmc_ll.dat");

  llik.resize(nspl);
  double check;
  std::string line;
  for (int j = 0; j < nspl; j++ ) {
    std::getline(input_file, line);
    check = std::atof(line.c_str());
    if (check == 0 || check < -pow(10, 300)) {
      llik[j] = -pow(10, 300);
    } else {
      llik[j] = check;
    }
  }


  /* Compute initial log priors, if applicable */
  std::string lp_stream = "./tmcmc_getLP.sh "+std::to_string(nProcs);
  system(lp_stream.c_str());
  std::ifstream lp_input("tmcmc_lp.dat");

  lprior.resize(nspl);
  for (int j = 0; j < nspl; j++ ) {
    std::getline(lp_input, line);
    check = std::atof(line.c_str());
    if (check == 0.0 || check < -pow(10, 300)) {
      lprior[j] = -pow(10, 300);
    } else {
      lprior[j] = check;
    }
  }


  /* Output first set of samples to file*/
  outProcToFile(spls,ndim,nspl,std::string("samples.dat.0"));
  outProcToFile(llik,1,   nspl,std::string("loglik.dat.0") );
  outProcToFile(lprior, 1,nspl,std::string("logprior.dat.0"));

  RealVector Sm;
  double accRatio = 1.0;
  int iter = 0;

  double beta = 0.0, dBeta = 0.0, evid = 0.0;
  double max_llik;

  do { // Start algorithm
    iter++;

    /* shuffle samples */
    shuffle_spls(spls,llik,lprior);

    /* compute weights */
    RealVector w(nspl,0.0);
    double wsum, wmean, w2mean, wstd;

    dBeta = std::min(BETA_MAX,1.0-beta);

    /* used to normalize log_likelihood values for weight computations */
    max_llik = *std::max_element(llik.begin(), llik.end());

    /* Adapt delta beta as needed */
    do {
      for (int j=0; j < nspl; j++) w[j] = exp(dBeta*(llik[j]-max_llik));
      wsum   = std::accumulate(w.begin(), w.end(), 0.0);
      wmean  = wsum / w.size();
      w2mean = std::inner_product(w.begin(), w.end(), w.begin(), 0.0)/ w.size();
      wstd   = sqrt(w2mean- pow(wmean, 2));

      if (wstd/wmean > (cv + 1.0) || wstd == 0) dBeta *= 0.9;
      if (wstd/wmean > (cv + 0.5) || wstd == 0) dBeta *= 0.95;
      if (wstd/wmean > (cv + 0.05) || wstd == 0) dBeta *= 0.99;
      if (wstd/wmean > (cv + 0.005) || wstd == 0) dBeta *= 0.999;
      if (wstd/wmean > (cv + 0.0005) || wstd == 0) dBeta *= 0.9999;
      if (wstd/wmean > (cv + 0.00005) || wstd == 0) dBeta *= 0.99999;
      if (wstd/wmean > (cv + 0.000005) || wstd == 0) dBeta *= 0.999999;
      if (wstd/wmean > (cv + 0.0000005) || wstd == 0) dBeta *= 0.9999999;
      if (wstd/wmean > (cv + 0.00000005) || wstd == 0) dBeta *= 0.99999999;

      if (dBeta < 1.0e-10)
        break;

    } while (wstd/wmean > (cv + 0.00000005) || wstd == 0);

    if (write_flag == 1){
      std::cout<<"DBeta: " << dBeta<<" Wmean: "<<wmean;
      std::cout<<" Wstd: "<<wstd<<" Cv: "<<wstd/wmean<<std::endl;
    }

    beta += dBeta;
    evid += log(wmean);
    if (write_flag == 1){
      std::cout<<"Iteration "<<iter<<" Beta= "<<beta;
      std::cout<<" wMean= "<< wmean << " Evid=   " << evid<<std::endl<<std::flush;
    }


    /* Save mean ll for Bayes factor */
    Sm.push_back(wmean);

    /* rescale w and do cumulative sum */
    RealVector wb(nspl);
    std::transform(w.begin(),w.end(),w.begin(),
                   std::bind2nd(std::multiplies<double>(), 1.0/wsum));
    std::partial_sum(&w[0],&w[0]+nspl,&wb[0]);
    wb.insert ( wb.begin() , 0.0 );

    /* Covariance matrix */
    RealVector theta0(ndim);
    for (int i=0; i<ndim; i++ ) {
      theta0[i] = 0.0;
      for (int j=0; j<nspl; j++ ) {
        theta0[i] += spls[j*ndim+i]*w[j];
      }
    }

    RealVector cvmat(ndim*ndim,0.0);
    for (int j=0; j < nspl; j++) {
      for (int i1=0; i1<ndim; i1++) {
        for (int i2=0; i2<ndim; i2++) {
          double outp=w[j] * (spls[j*ndim+i1]-theta0[i1]) *
                            (spls[j*ndim+i2]-theta0[i2]);
          cvmat[i1*ndim+i2] += outp;
          cvmat[i2*ndim+i1] += outp;
  	    }
      }
    }

    /* Control Parameter, Covariance rescaling */
    std::transform(cvmat.begin(), cvmat.end(), cvmat.begin(),
                   std::bind2nd(std::multiplies<double>(),gm2));


    /* Cholesky factorization of the proposal covariance, in-place */
    int chol_info=0;
    char lu='L';
    dpotrf_(&lu, &ndim, &cvmat[0], &ndim, &chol_info);

    /* generate random samples into [0,1] */
    RealVector spl01(nspl);
    for (int j=0; j<nspl; j++)
      spl01[j] = dsfmt_genrand_open_open(&RandomState);

    /* get bin IDs and count */
    IntVector pos(nspl,0);
    for (int j=0; j < nspl; j++) {
      pos[j] = std::lower_bound(wb.begin(),wb.end(),spl01[j])-wb.begin()-1;
    }

    /* Count number of times a sample is "picked" by PRNG */
    IntVector splPos, splCount;
    for (int j=0; j < nspl; j++) {
      int icount=0;
      for (int ispl=0; ispl<nspl; ispl++) {
        if (pos[ispl]==j) icount += MFactor;
      }
      if (icount>0) {
        splPos.push_back(j);
        splCount.push_back(icount);
      }
    }

    /* Initialize samples that were retained, cardinality, and
    likelihood values */
    RealVector splSt, llikSt, lpriorSt, gradSt;
    IntVector splCard;
    int nsplSt = splPos.size();

    /* Resampling Step */
    for (int ispl=0; ispl<nsplSt; ispl++) {
      if (basis) { // Resample according to BASIS and CATMIPs
        int isplCount = splCount[ispl];
        for (size_t i = 0; i < isplCount; ++i) {
          for (int j = 0; j < ndim; ++j) {
            splSt.push_back(spls[splPos[ispl]*ndim + j]);
          }
          splCard.push_back(CATSteps);
          llikSt.push_back(llik[splPos[ispl]]);
          lpriorSt.push_back(lprior[splPos[ispl]]);
        }
      } else {

        for (int i=0; i<ndim; i++ ) {
          splSt.push_back(spls[splPos[ispl]*ndim+i]);
        }
        splCard.push_back(splCount[ispl]);
        llikSt.push_back(llik[splPos[ispl]]);
        lpriorSt.push_back(lprior[splPos[ispl]]);
      }
    }


    RealVector splSave, llikSave, lpriorSave, XiSave, gradSave;
    int nSteps = *std::max_element(splCard.begin(), splCard.end());

    /* Run single steps of the Markov chains at a time, for the chains
        that need several jumps */
    if (basis) nsplSt = nspl; // Post resampling, nspl # of chains
    for (int isbSteps=0; isbSteps < nSteps; isbSteps++ ) {
      RealVector splCand(nsplSt*ndim);
      for (int ispl=0; ispl<nsplSt; ispl++) {
        /* generate candidate */
        RealVector xi(ndim);
        for (size_t i=0; i < ndim; i++)  {
          xi[i] = dsfmt_genrand_nrv(&RandomState);
        }

  	    for (size_t i=0; i < ndim; i++) {
          splCand[ispl*ndim+i] = splSt[ispl*ndim+i];
          double Lnrv=0.0;

          for (size_t j=0; j < (i+1); ++j) {
              Lnrv += cvmat[j*ndim+i] * xi[j];
          }

          splCand[ispl*ndim+i] += Lnrv;
        } /* done generating candidate */
      }


      /* Compute new likelihoods */
      RealVector splsComp;
      int compCount=0;
      for (int ispl=0; ispl<nsplSt; ispl++) {
        for (int i=0; i<ndim; i++ ) {
          splsComp.push_back(splCand[ispl*ndim+i]);
        }
        compCount++;
      }

      outProcToFile(splsComp,ndim,compCount,nProcs);
      std::string ll_stream="./tmcmc_getLL.sh "+ std::to_string(nProcs);
      system(ll_stream.c_str());
      std::ifstream input_file("tmcmc_ll.dat");

      /* Collect loglikelihood for proposals as applicable */
      RealVector llikComp(compCount);
      for (int j = 0; j < compCount; j++ ) {
        std::getline(input_file, line);
        check = std::atof(line.c_str());
        if (check == 0.0 || check < -pow(10, 300)) {
          llikComp[j] = -pow(10, 300);
        } else {
          llikComp[j] = check;
        }
      }

      std::string lp_stream = "./tmcmc_getLP.sh "+std::to_string(nProcs);
      system(lp_stream.c_str());
      std::ifstream lp_input("tmcmc_lp.dat");

      /* Collect logprior for proposals as applicable */
      RealVector lpriorComp(compCount);
      for (int j = 0; j < compCount; j++ ) {
        std::getline(lp_input, line);
        check = std::atof(line.c_str());
        // input_file >> check2;
        if (check == 0.0 || check < -pow(10, 300)) {
          lpriorComp[j] = -pow(10, 300);
        } else {
          lpriorComp[j] = check;
        }
      }

      /* decide who jumps */
      int icomp=0;
      int acceptCount = 0;
      RealVector splNew(nsplSt*ndim), llikNew(nsplSt), lpriorNew(nsplSt);
      RealVector gradNew(nsplSt*ndim);

      for (int ispl=0; ispl<nsplSt; ispl++) {
        double alpha = dsfmt_genrand_urv(&RandomState);
        double AcceptRatio = -1;

        AcceptRatio = beta * (llikComp[icomp] - llikSt[ispl])
                      + (lpriorComp[icomp] - lpriorSt[ispl]);

        if (log(alpha) < AcceptRatio) { // Accept proposal
          for (int i=0; i < ndim; i++) {
            splNew[ispl*ndim+i] = splsComp[icomp*ndim+i];
          }
          lpriorNew[ispl] = lpriorComp[icomp];
          llikNew[ispl] = llikComp[icomp];
          acceptCount++;

        } else { // Reject Proposal
          for (int i=0; i<ndim; i++ ) {
            splNew[ispl*ndim+i] = splSt[ispl*ndim+i] ;
          }
          lpriorNew[ispl] = lpriorSt[ispl];
          llikNew[ispl] = llikSt[ispl];
        }

        icomp++;
      }

      /* Save Samples for Next Iteration */
      if ((!basis && isbSteps % MFactor == 0) ||
          (basis && isbSteps == (nSteps - 1))) {

        for (int ij=0; ij<nsplSt*ndim; ij++) {
          splSave.push_back(splNew[ij]);
        }


        for (int ij=0; ij<nsplSt; ++ij) {
          llikSave.push_back(llikNew[ij]);
          lpriorSave.push_back(lpriorNew[ij]);
        }
      }


      /* Clear Proposals for next iteration */
      splSt.clear();
      gradSt.clear();
      llikSt.clear();
      lpriorSt.clear();

      gamma_file << "====================" << "\n";
      gamma_file << "Itera: " << iter << "\n";
      gamma_file << "Ratio: " << accRatio << "\n";
      gamma_file << "Gamma: " << gm2      << "\n";

      if (basis) {
        splSt = splNew;
        llikSt = llikNew;
        lpriorSt = lpriorNew;
        gradSt = gradNew;
      } else {
        for (int ispl=0; ispl<nsplSt; ispl++) {
          splCard[ispl] -= 1;
          if (splCard[ispl]>0) {
            for (int i=0; i<ndim; i++ ) {
              splSt.push_back(splNew[ispl*ndim+i]);
            }

            llikSt.push_back(llikNew[ispl]);
            lpriorSt.push_back(lpriorNew[ispl]);
  	      }
        }
      }

      /* Reduce length of chains remaining in samples */
      for (int ispl=0; ispl<splCard.size();) {
        if (splCard[ispl]==0)
          splCard.erase(splCard.begin()+ispl) ;
        else
          ispl++;
      }

      accRatio = (double) acceptCount / (double) nsplSt;

      nsplSt = llikSt.size();
      assert(splCard.size()==llikSt.size());
      assert(splSt.size()==llikSt.size()*ndim);
    }

    // Rescaling based on Catanach (Thesis 2017)
    // Assumptions made on optimal acceptance rate
    double G = 2.1;
    gm = gm * exp(G * (accRatio - 0.234));
    gm2 = pow(gm, 2.0);

    /* Set samples for next temperature iteration */
    spls = splSt;
    llik = llikSt;
    lprior = lpriorSt;
    assert(llik.size()==nspl);
    assert(lprior.size()==nspl);
    assert(spls.size()==nscal);

    outProcToFile(splSt,ndim,nspl,
                    std::string("samples.dat.")+std::to_string(iter));
    outProcToFile(llikSt,1,   nspl,
                    std::string("loglik.dat." )+std::to_string(iter));

    outProcToFile(lpriorSt,1, nspl,
                    std::string("logprior.dat.") + std::to_string(iter));

  } while ( ( iter < 1000 ) && (beta<1-1.e-10) );

  if (write_flag == 1){
    std::cout << "TMCMC Algorithm Done" << std::endl;
  }
  return (evid);
} /* done tmcmc */

void outProcToFile(const RealVector spls, const int ndim,
                    const int nspl, int nprocs) {
  // Separate spls vector into nproc pieces, into mcmcstates files

  /* no. of mcmc states per file */
  int nsplP = (int) nspl/nprocs;
  int nAdd  = nspl-nsplP*nprocs;

  /* save samples to files */
  int isplSt ;
  int isplEn = 0;
  for (int ifile=0; ifile < nprocs; ifile++) {

    isplSt  = isplEn;
    isplEn += nsplP;
    if (nAdd>0) {
      isplEn += 1;
      nAdd   -= 1;
    }

    char fname[20];
    sprintf(fname,"%s%d%s","mcmcstates_",ifile+1,".dat");
    FILE *myfile  = fopen(fname,"w") ;
    for (int j = isplSt; j < isplEn; j++) {
      for (int i = 0; i < ndim; i++)
        fprintf(myfile,"%24.18e ",spls[j*ndim+i]);
      fprintf(myfile,"\n");
    }
    fclose(myfile);
  }
  assert(isplEn==nspl);

  return ;

}

void outProcToFile(const RealVector spls, const int ndim, const int
nspl, std::string fname) {
  // Output sample vector into fname file
  // assert(spls.size()==ndim*nspl);

  FILE *myfile  = fopen(fname.c_str(),"w") ;
  for (int j = 0; j < nspl; j++) {
    for (int i = 0; i < ndim; i++)
      fprintf(myfile,"%24.18e ",spls[j*ndim+i]);
    fprintf(myfile,"\n");
  }
  fclose(myfile);

  return ;

}

void shuffle_spls(RealVector &spls, RealVector &llik, RealVector &lprior) {
  // Shuffle the samples randomly
  int nspl = llik.size();
  int ndim = spls.size()/nspl;

  IntVector idx(nspl);
  for (int j=0; j<nspl; j++) idx[j]=j;

  RealVector splsTmp(nspl*ndim), llikTmp(nspl), lpriorTmp(nspl);
  shuffle (idx.begin(), idx.end(), std::default_random_engine());
  for (int j = 0; j < nspl; j++) {
    llikTmp[j] = llik[idx[j]];
    lpriorTmp[j] = lprior[idx[j]];
    for (int i = 0; i < ndim; i++) splsTmp[j*ndim+i] = spls[idx[j]*ndim+i];
  }

  for (int j = 0; j < nspl; j++) {
    llik[j] = llikTmp[j];
    lprior[j] = lpriorTmp[j];
    for (int i = 0; i < ndim; i++) spls[j*ndim+i] = splsTmp[j*ndim+i];
  }

  return ;

}

bool fileExists(std::string fname) {
  // Check if a file exists
  std::ifstream file(fname);
  return file.is_open();
}

void readInitSamples(RealVector &spls, std::string fname) {
  // Read fname file and put it in spls vector
  std::string line;
  std::string token;
  std::ifstream DAT;
  std::stringstream iss;
  DAT.open(fname);
  double temp;
  std::string::size_type sz;
  int i = 0;

  while (std::getline(DAT, line)) {

    while(1){
      try {
        temp = std::stod (line,&sz);
      } catch (const std::exception& ex) {
        break;
      }
      line = line.substr(sz);
      spls.push_back(temp);
    }
  }
}

void parseSetup(dsfmt_t &RandomState, int nspl, int iseed) {
  std::ifstream input_file("setup.dat");
  // Read setup.dat and separate distribution and parameters
  std::string line;
  std::string token;
  CharVector distr;
  RealVector paramA;
  RealVector paramB;
  std::stringstream iss;

  while (std::getline(input_file, line)) {
    int step = 0;
    iss << line;
    while (std::getline(iss, token, ' ')) {
      switch (step) {
        case 0:
          if (token[0] != 'n' && token[0] != 'u') {
            std::cout << "Too many parameters for distribution" << std::endl;
            throw("Improperly formatted setup file");
          } else {
            distr.push_back(token[0]);
            step++;
          }
          break;

        case 1:
          paramA.push_back(std::stod(token));
          step++;
          break;

        case 2:
          paramB.push_back(std::stod(token));
          step++;
          if (distr.back() == 'n' || distr.back() == 'u') {
            step = 0;
          }
          break;

        default:
          break;
      }
    }
    line.clear();
    iss.clear();
  }

  PriorGen(RandomState, nspl, distr, paramA, paramB, iseed);
}

void PriorGen(dsfmt_t &RandomState, int nspl, CharVector &distr,
  RealVector &means, RealVector &vars, int iseed) {
  // Read sample.dat.0, generate starting random variates based on
  // parameters (means, vars / [a,b]) and distr
  std::ofstream DAT;
  DAT.open("samples.dat.0");
  dsfmt_init_gen_rand(&RandomState, iseed);
  DAT << std::setprecision(18);
  for (int j = 0; j < nspl; ++j) {
    for (size_t i = 0; i < distr.size(); ++i) {
      if (i != 0) {
        DAT << " ";
      }
      double rv;

      switch (distr[i]) {
        case 'u':
          rv = dsfmt_genrand_urv_sm(&RandomState, means[i], vars[i]);
          DAT << rv;
          break;

        case 'n':
          rv = dsfmt_genrand_nrv_sm(&RandomState, means[i], vars[i]);
          DAT << rv;
          break;
      }

    }
    DAT << "\n";
  }
  DAT.close();


}

double pearsonCorrCoef(RealVector X, RealVector Y) {
    RealVector Xtmp(X.size()), Ytmp(Y.size());
    assert(X.size() == Y.size());
    double xmean = std::accumulate(X.begin(), X.end(), 0.0) / X.size();
    double ymean = std::accumulate(Y.begin(), Y.end(), 0.0) / X.size();

    std::transform(X.begin(), X.end(), Xtmp.begin(),
                    std::bind2nd(std::minus<double>(), xmean));

    std::transform(Y.begin(), Y.end(), Ytmp.begin(),
                    std::bind2nd(std::minus<double>(), ymean));
    double sx = sqrt((1.0 / (X.size() - 1)) * std::inner_product(Xtmp.begin(),
                    Xtmp.end(), Xtmp.begin(), 0.0));
    double sy = sqrt((1.0 / (Y.size() - 1)) * std::inner_product(Ytmp.begin(),
                    Ytmp.end(), Ytmp.begin(), 0.0));

    double num = (std::inner_product(X.begin(), X.end(), Y.begin(), 0.0)
                    - X.size() * xmean * ymean);
    double denom = (X.size() - 1) * sx * sy;

    return (num / denom);
}

double rescaleSTD(RealVector w, double wmean) {
  RealVector wTmp(w.size());

  std::transform(w.begin(), w.end(), wTmp.begin(),
                  std::bind2nd(std::minus<double>(), wmean));

  double lTerm = std::inner_product(wTmp.begin(), wTmp.end(), wTmp.begin(), 0.0);
  double rTerm = pow(std::accumulate(wTmp.begin(), wTmp.end(), 0.0), 2.0) / w.size();

  double var = (lTerm - rTerm) / (w.size() - 1);

  return (sqrt(var));
}
