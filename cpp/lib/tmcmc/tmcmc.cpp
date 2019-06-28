#include "tmcmc.h"

#define BETA_MAX 0.8

extern "C" {
 void dpotrf_(char *, int *, double *, int*, int *);
}

void outProcToFile(const RealVector spls, const int ndim, const int nspl,
                    const int nprocs) ;
void outProcToFile(const RealVector spls, const int ndim, const int nspl,
                    std::string fname);
void shuffle_spls(RealVector &spls, RealVector &llik);

bool fileExists(std::string fname);
void readInitSamples(RealVector &spls, std::string fname);
void parseSetup(dsfmt_t &RandomState, int nspl, int iseed);
void PriorGen(dsfmt_t &RandomState, int nspl, CharVector &distr,
  RealVector &means, RealVector &vars, int iseed);
double pearsonCorrCoef(RealVector X, RealVector Y);
double rescaleSTD(RealVector w, double wmean);

double tmcmc(RealVector &rngs, double gm, int nspl, RealVector dtvec,
            int iseed, const int nProcs, int ndim = 2, bool spc = false) {
  // Wrapper for systematic tests
  double cv = 1.0;
  double a = 1.0/9.0;
  double b = 8.0/9.0;
  double mala = false;
  double tauScale = 0.0001;
  double betaThres = 0.0;
  int MFactor = 1;
  bool basis = false;
  int CATSteps = 1;
  return tmcmc(rngs, gm, nspl, dtvec, iseed, nProcs, ndim, spc, cv, a, b,
                mala, tauScale, betaThres, MFactor, basis, CATSteps);
}

double tmcmc(RealVector &rngs, double gm, int nspl, RealVector dtvec,
          int iseed, const int nProcs, int ndim, bool spc, double cv,
          double a, double b, bool mala, double tauScale, double betaThres,
          int MFactor, bool basis, int CATSteps) {
  /* TMCMC Algorithm
      Input: rngs - ranges for all samples
             gm - initial gamma value
             nspl - number of Samples
             dtvec - change in Beta (optional)
             iseed - random seed
             nProcs - number of processors
             ndim - number of dimensions
             spc - SpectralClustering (on/off)
             cv - Coefficient of Variance threshold for adapting Beta
             a - gamma multiplier (a + bR)gm
             b - gamma multiplier (a + bR)gm, See Minson, Simon, Beck 2013
             mala - MALA Proposal (on/off)
             tauScale - tau scaling parameter (MALA), l > 0
             betaThres - Turns off MALA when beta is above betaThres
             MFactor - Multiplicative factor for chain length to encourage mixing.
      Output: evid - asymptotically unbiased model evidence estimator
  */
  #ifdef USE_HDF5
  using namespace H5;
  #endif

  /* Initial Tau for MALA */
  if (tauScale < 0) {
    tauScale = 0;
  }
  double tau = pow(tauScale, 2) / (2.0 * pow(ndim, 1.0 / 3.0));

  /* Read ndims from rngs */
  if (rngs.size() != 0) {
    ndim = rngs.size() / 2;
  }
  int nscal = ndim*nspl;

  std::string nstages;
  if (dtvec.size() == 0) {
    nstages = "Adaptive";
  } else {
    nstages = std::to_string(dtvec.size());
  }
  //double dt = 1.0 / Ntemp;
  //double dts[] = {0.02,0.03,0.04,0.05,0.06,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1};
  //RealVector dtvec(dts, dts + sizeof(dts) / sizeof(double) );

  // Roberts and Rosenthal 2011 - Initial gamma if one is not provided.
  if (gm < 0) gm = 2.38 / sqrt(ndim);
  double gm2 = gm * gm;

  std::ofstream gamma_file("gamma.dat");

  std::cout<<"----------------------------------------------------"<<std::endl;
  std::cout<<"No. of samples    :"<<nspl<<std::endl;
  std::cout<<"No. of dimensions :"<<ndim<<std::endl;
  std::cout<<"No. of stages     :"<<nstages<<std::endl;
  std::cout<<"MALA              :"<<mala<<std::endl;
  std::cout<<"BASIS/CATMIPs     :"<<basis<<std::endl;
  if (mala)
  std::cout<<"Tau               :"<<tau<<std::endl;
  if (basis)
  std::cout<<"CATSteps          :"<<CATSteps<<std::endl;
  std::cout<<"MFactor           :"<<MFactor<<std::endl;
  std::cout<<"----------------------------------------------------"<<std::endl;

  /* initialize random number generator */
  dsfmt_t RandomState;
  dsfmt_init_gen_rand(&RandomState,iseed);

  RealVector spls;
  // Track if a setup or sample file was used.
  bool usedSetup = false;
  bool usedSample = false;

  #ifdef USE_HDF5
  /* Open up H5 File for use */
  H5File h5file("RunResults.h5", H5F_ACC_RDWR);
  h5file.openFile("RunResults.h5", H5F_ACC_RDWR);
  hsize_t attrDims[1] = {1};
  hsize_t sampleDims[2] = {nspl, ndim};
  hsize_t loglikDims[2] = {nspl, 1};
  DataSpace attrDataSpace = DataSpace (1, attrDims);
  DataSpace sampleDataSpace = DataSpace (2, sampleDims);
  DataSpace logDataSpace = DataSpace (2, loglikDims);

  /* Generate stage group in file */
  int i;
  std::ifstream TMCMCiter("TMCMCiter.dat");
  TMCMCiter >> i;
  TMCMCiter.close();
  std::string stageName = "Stage" + std::to_string(i);
  Group stage = h5file.openGroup(stageName.c_str());

  /* Extract previous set of samples */
  if (i > 0) {
    int lastSampleIter;
    std::ifstream lastSample("lastSample.dat");
    lastSample >> lastSampleIter;
    lastSample.close();
    std::string oldStageName = "Stage" + std::to_string(i - 1);
    std::string oldIterName = "Iteration_" + std::to_string(lastSampleIter);
    Group oldStage = h5file.openGroup(oldStageName.c_str());
    Group oldIter = oldStage.openGroup(oldIterName.c_str());
    DataSet oldSamples = oldIter.openDataSet("Samples");
    double *sampleArray = new double[nspl*ndim];
    oldSamples.read(sampleArray, PredType::NATIVE_DOUBLE);
    for (int i = 0; i < nspl; ++i) {
      for (int j = 0; j < ndim; ++j) {
        spls.push_back(sampleArray[i*ndim + j]);
      }
    }
    delete[] sampleArray;
    usedSample = true;
    std::cout << "Samples read from old samples in HDF5\n";
  } else if (fileExists("setup.dat")){
    parseSetup(RandomState, nspl, iseed);
    readInitSamples(spls, "samples.dat.0");
    std::cout << "Samples generated from setup.dat, then deleted, (HDF5)\n";
    usedSetup = true;
  } else {
    std::cout << "Setup file needed for HDF5!";
    exit(0);
  }
  #else

  /* Generate Initial Random Samples */
  /* Check if an initial sample file exists, or a setup file
     last resort, generate from a d-dim hypercube */
  if (fileExists("samples.dat.0")) {
    readInitSamples(spls, "samples.dat.0");
    std::cout << "Samples read from samples.dat.0\n";
    usedSample = true;

  } else if (fileExists("setup.dat")) {
    parseSetup(RandomState, nspl, iseed);
    readInitSamples(spls, "samples.dat.0");
    std::cout << "Samples generated from setup.dat\n";
    usedSetup = true;

  } else {
    /* Generate from d-dim hypercube */
    for (int j=0; j<nscal; j++) {
      spls.push_back(dsfmt_genrand_open_open(&RandomState));
    }
    /* Rescale to ranges */
    for (int j=0; j<nspl; j++ ) {
      int ipos=j*ndim;
      for (int i=0; i<ndim; i++ ) {
        spls[ipos+i] = rngs[2*i]+spls[ipos+i]*(rngs[2*i+1]-rngs[2*i]);
      }
    }
    std::cout << "Samples generated from d-dim hypercube\n";
  }
  #endif

  assert(spls.size() == nscal);

  /* Compute initial likelihoods using outside model*/
  outProcToFile(spls, ndim, nspl, nProcs); // Make mcmcstates
  std::string ll_stream="./tmcmc_getLL.sh "+std::to_string(nProcs);
  system(ll_stream.c_str());
  std::ifstream input_file("tmcmc_ll.dat");


  RealVector llik(nspl);
  double check;
  for (int j = 0; j < nspl; j++ ) {
    input_file >> check;
    if (check == 0 || check < -pow(10, 300)) {
      llik[j] = -pow(10, 300);
    } else {
      llik[j] = check;
    }
  }

  RealVector gradLog;
  if (mala) {
    std::ifstream grad_input("gradlog.dat");
    double check4;
    while (grad_input >> check4) {
      // Just in case of empty spaces in file
      if (check4 == 0) {
        continue;
      }
      gradLog.push_back(check4);
    }
  }

  /* Compute initial log priors, if applicable */
  RealVector lprior(nspl);
  if (usedSetup || usedSample) {
    std::string lp_stream = "./tmcmc_getLP.sh " + std::to_string(0);
    system(lp_stream.c_str());
    std::ifstream lp_input("tmcmc_lp.dat");
    for (int i = 0; i < nspl; ++i ) {
      lp_input >> lprior[i];
    }
  }


  #ifdef USE_HDF5
  /* Put Samples into HDF5 */
  std::string iterName = "Iteration_0";
  Group iterGroup = stage.createGroup(iterName.c_str());

  DataSet sampleDataset = iterGroup.createDataSet("Samples",
                            PredType::NATIVE_DOUBLE, sampleDataSpace);
  DataSet loglikDataset = iterGroup.createDataSet("LogLikelihood",
                            PredType::NATIVE_DOUBLE, logDataSpace);
  DataSet logpriorDataset = iterGroup.createDataSet("LogPrior",
                            PredType::NATIVE_DOUBLE, logDataSpace);

  double *splsDatasetData = new double[spls.size()];
  double *llikDatasetData = new double[llik.size()];

  for (size_t i = 0; i < spls.size(); ++i) {
    splsDatasetData[i] = spls[i];
  }
  for (size_t i = 0; i < llik.size(); ++i) {
    llikDatasetData[i] = llik[i];
  }

  sampleDataset.write(splsDatasetData, PredType::NATIVE_DOUBLE);
  loglikDataset.write(llikDatasetData, PredType::NATIVE_DOUBLE);


  if (usedSetup || usedSample) {
    double *lpriorDatasetData = new double[lprior.size()];
    for (size_t i = 0; i < lprior.size(); ++i) {
      lpriorDatasetData[i] = lprior[i];
    }

    logpriorDataset.write(lpriorDatasetData, PredType::NATIVE_DOUBLE);
    delete[] lpriorDatasetData;
  }

  Attribute betaAttr = iterGroup.createAttribute("Beta",
                          PredType::NATIVE_DOUBLE, attrDataSpace);
  double betaAttrData[1] = {0.0};
  betaAttr.write(PredType::NATIVE_DOUBLE, betaAttrData);
  iterGroup.close();
  #endif

  /* Output first set of samples to file*/
  outProcToFile(spls,ndim,nspl,std::string("samples.dat.0"));
  outProcToFile(llik,1,   nspl,std::string("loglik.dat.0") );
  if (usedSetup || usedSample) {
    outProcToFile(lprior, 1,nspl,std::string("logprior.dat.0"));
  }

  RealVector Sm;
  double accRatio = 1.0;
  double pearsonR = 1.0;
  int iter = 0;
  double beta = 0.0, dBeta = 0.0, evid = 0.0;

  do { // Start algorithm
    iter++;

    /* shuffle samples */
    shuffle_spls(spls,llik);

    /* compute weights */
    RealVector w(nspl,0.0);
    double wsum, wmean, w2mean, wstd;
    if (dtvec.size()>=iter)
      dBeta = std::min(dtvec[iter-1],1.0-beta);
    else
      dBeta = std::min(BETA_MAX,1.0-beta);

    /* Adapt delta beta as needed */
    do {
      for (int j=0; j < nspl; j++) w[j] = exp(dBeta*llik[j]);
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


    } while (wstd/wmean > (cv + 0.00000005) || wstd == 0);

    std::cout<<"DBeta: " << dBeta<<" Wmean: "<<wmean;
    std::cout<<" Wstd: "<<wstd<<" Cv: "<<wstd/wmean<<std::endl;

    beta += dBeta;
    evid += log(wmean);
    std::cout<<"Iteration "<<iter<<" Beta= "<<beta;
    std::cout<<" wMean= "<< wmean << " Evid=   " << evid<<std::endl<<std::flush;

    /* Turn off MALA after threshold. Theory is not understood */
    if (beta >= betaThres) {
      mala = false;
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
    if (spc) {
      // Incomplete - need to link individual covariance to clusters.
      std::string cluster_string = "./cluster.sh " + std::to_string(iter - 1);
      system(cluster_string.c_str());
      std::ifstream covfile("Cov.dat");
      for (size_t i = 0; i < ndim * ndim; ++i) covfile >> cvmat[i];
    } else {
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
    }

    /* Control Parameter, Covariance rescaling */
    if (!mala) {
      std::transform(cvmat.begin(), cvmat.end(), cvmat.begin(),
                   std::bind2nd(std::multiplies<double>(),gm2));
    }

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
          if (mala) {
            gradSt.push_back(gradLog[splPos[ispl] * ndim + i]);
          }
        }
        splCard.push_back(splCount[ispl]);
        llikSt.push_back(llik[splPos[ispl]]);
        lpriorSt.push_back(lprior[splPos[ispl]]);
      }
    }

    RealVector splSave, llikSave, lpriorSave, XiSave, gradSave;
    int nSteps = *std::max_element(splCard.begin(), splCard.end());
    std::cout << "Max Steps: " << nSteps << std::endl;

    /* Run single steps of the Markov chains at a time, for the chains
        that need several jumps */
    if (basis) nsplSt = nspl; // Post resampling, nspl # of chains
    for (int isbSteps=0; isbSteps < nSteps; isbSteps++ ) {
      RealVector splCand(nsplSt*ndim);
      BoolVector isInside(nsplSt,true);
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

          if (mala) {
            XiSave.push_back(Lnrv); // Save for use in AcceptRatio
            splCand[ispl * ndim + i] += tau * gradSt[ispl*ndim + i] +
                                        sqrt(2.0 * tau) * Lnrv;
          } else {
            splCand[ispl*ndim+i] += Lnrv;
          }

          // Ensure within bounds
          if ((splCand[ispl*ndim+i] < rngs[2*i])
              || (splCand[ispl*ndim+i] > rngs[2*i+1])) {
            isInside[ispl] = false;
            break;
          }
        } /* done generating candidate */
      }

      /* Compute new likelihoods */
      RealVector splsComp;
      int compCount=0;
      for (int ispl=0; ispl<nsplSt; ispl++) {
        if (isInside[ispl]) {
          for (int i=0; i<ndim; i++ ) {
            splsComp.push_back(splCand[ispl*ndim+i]);
          }
          compCount++;
	      }
      }

      std::cout << "Jump: " << isbSteps;
      std::cout << ", Saving " << compCount << " out of ";
      std::cout << nsplSt << " to file" << std::endl;
      outProcToFile(splsComp,ndim,compCount,nProcs);
      std::string ll_stream="./tmcmc_getLL.sh "+ std::to_string(nProcs);
      system(ll_stream.c_str());
      std::ifstream input_file("tmcmc_ll.dat");

      /* Collect loglikelihood, gradlog, logprior for proposals as applicable */
      RealVector llikComp(compCount);
      double check2;
      for (int j = 0; j < compCount; j++ ) {
        input_file >> check2;
        if (check2 == 0.0 || check2 < -pow(10, 300)) {
          llikComp[j] = -pow(10, 300);
        } else {
          llikComp[j] = check2;
        }
      }

      RealVector gradComp;
      if (mala) {
        std::ifstream gradInput("gradlog.dat");
        double check3;
        while (gradInput >> check3) {
          if (check3 == 0) {
            continue;
          }
          gradComp.push_back(check3);
        }
      }

      RealVector lpriorComp(compCount);
      if (usedSetup || usedSample) {
        std::string lp_stream = "./tmcmc_getLP.sh " + std::to_string(iter - 1);
        system(lp_stream.c_str());
        std::ifstream lp_input("tmcmc_lp.dat");
        double check5;
        for (int i = 0; i < compCount; ++i) {
          lp_input >> check5;
          if (check5 == 0.0 || check5 < -pow(10, 300)) {
            lpriorComp[i] = -pow(10, 300);
          } else {
            lpriorComp[i] = check5;
          }
        }
      }

      /* decide who jumps */
      int icomp=0;
      int acceptCount = 0;
      RealVector splNew(nsplSt*ndim), llikNew(nsplSt), lpriorNew(nsplSt);
      RealVector gradNew(nsplSt*ndim);

      for (int ispl=0; ispl<nsplSt; ispl++) {
        if (isInside[ispl]) { // Require within bounds
          double alpha = dsfmt_genrand_urv(&RandomState);
          double AcceptRatio = -1;

          if (!mala && (usedSetup || usedSample)) { // Standard MH Update
            AcceptRatio = beta * (llikComp[icomp] - llikSt[ispl])
                          + (lpriorComp[icomp] - lpriorSt[ispl]);

          } else if (mala) { // MALA - Asymmetry.
            if (usedSetup || usedSample) {
              AcceptRatio = beta * (llikComp[icomp] - llikSt[ispl])
                          + (lpriorComp[icomp] - lpriorSt[ispl]);
            } else {
              AcceptRatio = beta * (llikComp[icomp] - llikSt[ispl]);
            }

            for (size_t i = 0; i < ndim; ++i) {
              AcceptRatio += 0.5 * pow(XiSave[icomp * ndim + i], 2);
              AcceptRatio += -(1.0 / (4.0 * tau)) *
                             pow((-(tau * gradComp[icomp * ndim + i])
                                  -(tau * gradSt[ispl * ndim + i])
                                  -(sqrt(2.0 * tau) * XiSave[icomp*ndim+i])),2);
            }

          } else { // No Setup File, No MALA
            AcceptRatio = beta * (llikComp[icomp] - llikSt[ispl]);
          }

          if (log(alpha) < AcceptRatio) { // Accept proposal
            for (int i=0; i < ndim; i++) {
              splNew[ispl*ndim+i] = splsComp[icomp*ndim+i];
              if (mala) {
                gradNew[ispl*ndim+i] = gradComp[icomp*ndim+i];
              }
            }
            if (usedSetup || usedSample) {
              lpriorNew[ispl] = lpriorComp[icomp];
            }
            llikNew[ispl] = llikComp[icomp];
            acceptCount++;

	        } else { // Reject Proposal
            for (int i=0; i<ndim; i++ ) {
              splNew[ispl*ndim+i] = splSt[ispl*ndim+i] ;
              if (mala) {
                gradNew[ispl*ndim+i] = gradSt[ispl*ndim+i];
              }
            }
            if (usedSetup || usedSample) {
              lpriorNew[ispl] = lpriorSt[ispl];
            }
            llikNew[ispl] = llikSt[ispl];
	        }

          icomp++;
	      } else { // Out of Bounds - Automatic Reject
          for (int i=0; i<ndim; i++ ) {
            splNew[ispl*ndim+i] = splSt[ispl*ndim+i] ;
            if (mala) {
              gradNew[ispl*ndim+i] = gradSt[ispl*ndim+i];
            }
          }
          llikNew[ispl] = llikSt[ispl];
	      }
      }

      /* Save Samples for Next Iteration */
      if ((!basis && isbSteps % MFactor == 0) ||
          (basis && isbSteps == (nSteps - 1))) {

        for (int ij=0; ij<nsplSt*ndim; ij++) {
          splSave.push_back(splNew[ij]);
        }

        if (mala) {
          gradSave.clear();
          for (int ij=0; ij<nsplSt*ndim; ij++) {
            gradSave.push_back(gradNew[ij]);
          }
        }

        for (int ij=0; ij<nsplSt; ++ij) {
          llikSave.push_back(llikNew[ij]);
        }
        for (int ij = 0; ij < nsplSt; ++ij) {
            lpriorSave.push_back(lpriorNew[ij]);
        }
      }


      /* Clear Proposals for next iteration */
      splSt.clear();
      gradSt.clear();
      llikSt.clear();
      lpriorSt.clear();

      // Rescaling based on Muto & Beck 2008 (See Minson et al 2013)
      // gamma = a + b*R, R = current acceptance rate.
      // double a = (double)1.0/(double)9.0;
      // double b = (double)8.0/(double)9.0;
      // gm2 = pow(accRatio * b + a, 2);


      // Rescaling MALA - No literature.
      if (tau + 0.01 * (accRatio - 0.5) < 0) {
        tau = 0.00000001;
      } else {
        tau = tau + 0.01 * (accRatio - 0.5);
      }

      gamma_file << "====================" << "\n";
      gamma_file << "Itera: " << iter << "\n";
      gamma_file << "Ratio: " << accRatio << "\n";
      gamma_file << "Gamma: " << gm2      << "\n";
      gamma_file << "Tau  : " << tau      << "\n";

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
              if (mala) gradSt.push_back(gradNew[ispl*ndim+i]);
            }

            llikSt.push_back(llikNew[ispl]);

            if (usedSetup || usedSample) {
              lpriorSt.push_back(lpriorNew[ispl]);
            }
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

      pearsonR = pearsonCorrCoef(spls, splSt);
      std::cout << "Correlation: " << pearsonR << std::endl;
      nsplSt = llikSt.size();
      std::cout<<splCard.size()<<" - "<<llikSt.size()<<std::endl;
      assert(splCard.size()==llikSt.size());
      assert(splSt.size()==llikSt.size()*ndim);
    }

    // Rescaling based on Catanach (Thesis 2017)
    // Assumptions made on optimal acceptance rate
    double G = 2.1;
    gm = gm * exp(G * (accRatio - 0.234));
    gm2 = pow(gm, 2.0);

    /* Set samples for next temperature iteration */
    spls = splSave;
    gradLog = gradSave;
    llik = llikSave;
    lprior = lpriorSave;
    assert(llik.size()==nspl);
    assert(spls.size()==nscal);


    #ifdef USE_HDF5
    iterName = "Iteration_" + std::to_string(iter);
    iterGroup = stage.createGroup(iterName.c_str());

    sampleDataset = iterGroup.createDataSet("Samples",
                      PredType::NATIVE_DOUBLE, sampleDataSpace);
    loglikDataset = iterGroup.createDataSet("LogLikelihood",
                      PredType::NATIVE_DOUBLE, logDataSpace);
    logpriorDataset = iterGroup.createDataSet("LogPrior",
                      PredType::NATIVE_DOUBLE, logDataSpace);

    for (size_t i = 0; i < spls.size(); ++i) {
      splsDatasetData[i] = spls[i];
    }
    for (size_t i = 0; i < llik.size(); ++i) {
      llikDatasetData[i] = llik[i];
    }

    sampleDataset.write(splsDatasetData, PredType::NATIVE_DOUBLE);
    loglikDataset.write(llikDatasetData, PredType::NATIVE_DOUBLE);

    if (usedSetup || usedSample) {
      double *lpriorDatasetData = new double[lprior.size()];
      for (size_t i = 0; i < lprior.size(); ++i) {
        lpriorDatasetData[i] = lprior[i];
      }

      logpriorDataset.write(lpriorDatasetData, PredType::NATIVE_DOUBLE);
      delete[] lpriorDatasetData;
    }

    Attribute betaAttr = iterGroup.createAttribute("Beta",
                            PredType::NATIVE_DOUBLE, attrDataSpace);
    double betaAttrData[1] = {beta};
    betaAttr.write(PredType::NATIVE_DOUBLE, betaAttrData);
    iterGroup.close();
    #endif

    outProcToFile(spls,ndim,nspl,
                    std::string("samples.dat.")+std::to_string(iter));
    outProcToFile(llik,1,   nspl,
                    std::string("loglik.dat." )+std::to_string(iter));

    if (usedSetup || usedSample) {
      outProcToFile(lprior,1, nspl,
                    std::string("logprior.dat.") + std::to_string(iter));
    }

  } while ( ( iter < 1000 ) && (beta<1-1.e-10) );

  #ifdef USE_HDF5
  delete[] splsDatasetData;
  delete[] llikDatasetData;

  std::ofstream lastSample("lastSample.dat");
  lastSample.clear();
  lastSample << iter;
  lastSample.close();
  stage.close();
  h5file.close();
  #endif

  // Obsolete?
  // outProcToFile(spls,ndim,nspl,std::string("chain.dat"));
  // outProcToFile(llik,1,   nspl,std::string("loglk.dat"));

  // if (usedSetup || usedSample) {
  //   outProcToFile(lprior, 1, nspl, std::string("logpr.dat"));
  // }

  std::cout << "Algorithm Done" << std::endl;
  return (evid);
} /* done tmcmc */


void outProcToFile(const RealVector spls, const int ndim,
                    const int nspl, const int nprocs) {
  // Separate spls vector into nproc pieces, into mcmcstates files
  assert(spls.size()==ndim*nspl);

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

//void outProcToFile(const RealVector spls, const int ndim,
//                      const int nspl, const int nprocs) {
//
//  assert(spls.size()==ndim*nspl);
//
//  /* no. of mcmc states per file */
//  int nsplP = (int) nspl/nprocs;
//  if ( nsplP*nprocs < nspl) nsplP += 1;
//
//  /* save samples to files */
//  for (int ifile=0; ifile < nprocs; ifile++) {
//
//    int isplSt = ifile * nsplP;
//   int isplEn = isplSt+nsplP ;
//    if ( isplEn > nspl ) isplEn = nspl;
//
//    char fname[20];
//    sprintf(fname,"%s%d%s","mcmcStates_",ifile+1,".dat");
//    FILE *myfile  = fopen(fname,"w") ;
//    for (int j = isplSt; j < isplEn; j++) {
//      for (int i = 0; i < ndim; i++)
//        fprintf(myfile,"%24.18e ",spls[j*ndim+i]);
//      fprintf(myfile,"\n");
//    }
//    fclose(myfile);
//  }
//
//  return ;
//
//}

void outProcToFile(const RealVector spls, const int ndim, const int
nspl, std::string fname) {
  // Output sample vector into fname file
  assert(spls.size()==ndim*nspl);

  FILE *myfile  = fopen(fname.c_str(),"w") ;
  for (int j = 0; j < nspl; j++) {
    for (int i = 0; i < ndim; i++)
      fprintf(myfile,"%24.18e ",spls[j*ndim+i]);
    fprintf(myfile,"\n");
  }
  fclose(myfile);

  return ;

}

void shuffle_spls(RealVector &spls, RealVector &llik) {
  // Shuffle the samples randomly
  int nspl = llik.size();
  int ndim = spls.size()/nspl;

  IntVector idx(nspl);
  for (int j=0; j<nspl; j++) idx[j]=j;

  RealVector splsTmp(nspl*ndim), llikTmp(nspl);
  shuffle (idx.begin(), idx.end(), std::default_random_engine());
  for (int j = 0; j < nspl; j++) {
    llikTmp[j] = llik[idx[j]];
    for (int i = 0; i < ndim; i++) splsTmp[j*ndim+i] = spls[idx[j]*ndim+i];
  }

  for (int j = 0; j < nspl; j++) {
    llik[j] = llikTmp[j];
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

  while (std::getline(DAT, line)) {
    iss << line;
    while (std::getline(iss, token, ' ')) {
      spls.push_back(std::stod(token));
    }
    iss.clear();
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

  #ifdef USE_HDF5
  using namespace H5;
  H5File h5file("RunResults.h5", H5F_ACC_RDWR);
  hsize_t setupDim[2] = {paramA.size(), 1};
  hsize_t paramDim[2] = {paramA.size(), 2};
  DataSpace setupDataSpace = DataSpace (2, setupDim);
  DataSpace paramDataSpace = DataSpace (2, paramDim);

  DataSet setupDataset = h5file.createDataSet("Setup Distributions",
        PredType::NATIVE_CHAR, setupDataSpace);
  DataSet paramDataset = h5file.createDataSet("Setup Parameters",
        PredType::NATIVE_DOUBLE, paramDataSpace);

  char setupArray[paramA.size()];
  double paramArray[paramA.size()][2];
  for (size_t i = 0; i < paramA.size(); ++i) {
    setupArray[i] = distr[i];
    paramArray[i][0] = paramA[i];
    paramArray[i][1] = paramB[i];
  }

  setupDataset.write(setupArray, PredType::NATIVE_CHAR);
  paramDataset.write(paramArray, PredType::NATIVE_DOUBLE);
  setupDataset.close();
  paramDataset.close();
  h5file.close();
  #endif
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
