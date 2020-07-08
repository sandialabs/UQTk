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
/// \file mcmc.h
/// \author K. Sargsyan, C. Safta, B. Debusschere, 2012 -
/// \brief Header file for the Markov chain Monte Carlo class

#ifndef UQTKMCMC_H_SEEN
#define UQTKMCMC_H_SEEN


#include "dsfmt_add.h"

#include <iostream>
#include <string.h>
#include <stdio.h>
#include <sstream>

using namespace std; // needed for python string conversion

class LikelihoodBase{
public:
  virtual double eval(Array1D<double>&){
      return 3.14;
  };
  virtual ~LikelihoodBase(){};
};

class base{
public:
    virtual double fun(double* x, int n){return 0.0;}
};

class main{
public:
    base* b_;

    main(base& b){
        b_ = &b;
    }
};
//*****************************************

/// \class MCMC
/// \brief Markov Chain Monte Carlo base class.
///        Implemented the basic and most general MCMC algorithms
class MCMC{
public:
    // Constructors:
    
    /// \brief Constructor when given a pointer to a logPosterior function and a pointer to additional information (e.g. data)
    MCMC(double (*logposterior)(Array1D<double>&, void *), void *postinfo);
    /// \brief Constructor when given a Likelihood Base class object
    //MCMC(LikelihoodBase& L);
    /// \brief Dummy Constructor, used for TMCMC
    MCMC();
    
    // Destructor:
    ~MCMC(){};
    
    // Set or initialization functions:
    
    /// \brief Set the write flag function given an integer
    void setWriteFlag(int I);
    /// \brief Set the gradient function given a pointer to a logPosterior function, a 1D array of doubles, and a pointer to additional information (e.g. data)
    void setGradient(void (*gradlogPosterior)(Array1D<double>&, Array1D<double>&, void *));
    /// \brief Set the metric tensor function given a pointer to a metric tensor function, a 2D array of doubles, and a pointer to additional information (e.g. data)
    void setMetricTensor(void (*metricTensor)(Array1D<double>&, Array2D<double>&, void *));
    /// \brief Set the function accept and reject functions given a pointer to the accept and reject functions
    void setFcnAccept(void (*fcnAccept)(void *));
    void setFcnReject(void (*fcnReject)(void *));
    /// \brief Set the dimensions of the chain given an integer
    void setChainDim(int chdim);
    /// \brief Initialize proposal covariance matrix given as a 2d-array
    ///        For aMCMC, this matrix is used only before adaptivity starts
    void initChainPropCov(Array2D<double>& propcov);
    /// \brief Initialize proposal covariance matrix given its 1d-array diagonal
    ///        For aMCMC, this matrix is used only before adaptivity starts
    void initChainPropCovDiag(Array1D<double>& sig);
    /// \brief Set output specification struct, type('txt' or 'bin'), filename, frequency of outputs to the file and to screen.
    void setOutputInfo(string outtype, string file,int freq_file, int freq_screen);
    /// \brief Set the indicator to confirm that the names of parameters are prepended in the output file
    void namesPrepended();
    /// \brief Set random generation seed
    void setSeed(int seed);
    /// \brief Set lower bounds
    void setLower(double lower, int i);
    /// \brief Set upper bounds
    void setUpper(double upper, int i);
    /// \brief Set default unbounded domain
    void setDefaultDomain();

    
    // Get functions:
    
    /// \brief Returns proposal covariance matrix
    void getChainPropCov(Array2D<double>& propcov);
    /// \brief Get the name of the chain file
    string getFilename();
    /// \brief Get the value of the write flag as an integer
    int getWriteFlag();
    /// \brief Get samples of the chain with burnin and thining
    void getSamples(int burnin, int every,Array2D<double>& samples);
    /// \brief Get all samples of the chain
    void getSamples(Array2D<double>& samples);
    /// \brief Get gradient function by passing in log posterior pointer
    void getGradient(void (*gradlogPosterior)(Array1D<double>&, Array1D<double>&, void *));
    /// \brief Get metric tensor function
    void getMetricTensor(void (*metricTensor)(Array1D<double>&, Array2D<double>&, void *));
    /// \brief Get the accept and reject functions given a pointer
    void getFcnAccept(void (*fcnAccept)(void *));
    void getFcnReject(void (*fcnReject)(void *));
    /// \brief Get the output file type as a string
    string getOutputType();
    /// \brief Get the frequency of output to file
    int getFileFreq();
    /// \brief Get the frequency of output to the screen
    int getScreenFreq();
    /// \brief Get if the names are prepended
    bool getNamesPrepended();
    /// \brief Get the seed used for random generation
    int getSeed();
    /// \brief Get the lower bounds based on an index i
    double getLower(int i);
    /// \brief Get the upper bounds based on an index i
    double getUpper(int i);
    
    // Chain Functions:
    
    /// \brief Reset the chain state
    void resetChainState();
    /// \brief Reset to a new chain file
    void resetChainFilename(string filename);
    /// \brief An auxiliary function to parse the binary file and produce an array of chain-states
    void parseBinChain(string filename, Array1D<chainstate>& readchain);
    /// \brief Write an array of chain-states to a file
    void writeFullChainTxt(string filename, Array1D<chainstate> fullchain);
    /// \brief Get full chain as an array of chain-states
    void getFullChain(Array1D<chainstate>& readchain);
    /// \brief Append MAP state to the end
    void appendMAP();
    /// \brief Get MAP parameters
    double getMode(Array1D<double>& MAPparams);
    
    // Run functions:
    
    /// \brief The optimization routine
    void runOptim(Array1D<double>& start);
    /// \brief The actual function that generates MCMC
    virtual void runChain(int ncalls, Array1D<double>& chstart);
    /// \brief Start an MCMC chain with trivial initial condition
    void runChain(int ncalls);
    
    // Evaluation Functions:
    
    /// \brief Check to see if a new mode was found during last call to runChain
    bool newModeFound();
    /// \brief Get the chain's acceptance ratio
    void getAcceptRatio(double * accrat);
    /// \brief Get the chain's acceptance ratio as a double
    double getAcceptRatio();
    /// \brief Get the MCMC chain dimensionality
    int GetChainDim() const;
    /// \brief Function to evaluate the log-posterior
    double evalLogPosterior(Array1D<double>& m);
    /// \brief Function to evaluate the gradient of log-posterior
    void evalGradLogPosterior(Array1D<double>& m, Array1D<double>& grads);
    /// \brief Check if a point is in the domain
    bool inDomain(Array1D<double>& m);
    
    // Struct for the chain state
    struct chainstate{
      int step;
      Array1D<double> state;
      double alfa;
      double post;
    };
    
    
private:
    int WRITE_FLAG; // Write Flag
    int FLAG; // Flag
    LikelihoodBase* L_; // Pointer to the Likelihood base passed in through contstructor
    struct outputInfo{
        string outtype;
        string filename;
        int freq_file;
        int freq_screen;
    } outputinfo_;
    dsfmt_t RandomState; // (Was pubic before so we will need to see if it needs to be)
    int chainDim_; // Chain dimensions
    double (*logPosterior_)(Array1D<double>&, void *); // Pointer to log-posterior function
    void (*gradlogPosterior_)(Array1D<double>&, Array1D<double>&, void *); //Pointer to gradient log-posterior function
    void (*metricTensor_)(Array1D<double>&, Array2D<double>&, void *); // Pointer to metric tensorr function
    void (*fcnAccept_)(void *); // Pointer to accept function
    void (*fcnReject_)(void *); // Pointer to reject function
    void *postInfo_; // Void pointer to the posterior info (data)
    Array2D<double> chcov; // Chain proposal distributions (before the adaptivity starts)
    Array2D<double> propLCov_; // The Cholesky factor(square-root) of proposal covariance
    int seed_; // Random seed for MCMC
    
    /// \brief Pure virtual proposal function that will exist in all instances of the derived classes
    void virtual proposal() = 0;
    
    double probOldNew(Array1D<double>& a, Array1D<double>& b); // Evaluate old|new probabilities and new|old probabilities
    double evallogMVN_diag(Array1D<double>& x,Array1D<double>& mu,Array1D<double>& sig2); // Evaluate MVN
    
    chainstate currState_; // The current chain state
    chainstate modeState_; // The current MAP state
    Array1D<chainstate> fullChain_; // Array of chain states
    
    void updateMode(); // Function to update the chain mode
    
    void writeChainTxt(string filename); // Write the full chain as a text
    void writeChainBin(string filename); // Write the full chain as a binary file
    int lastwrite_; // Indicates up to which state 
    bool namesPrepend = false;
    
    bool newMode_ = false // Flag to indicate whether a new mode is found during the last call of runChain, initalized as false
    
    double accRatio_ = -1.0; // Acceptance ratio of the chain, initialized as -1.0
    
    // Flags to indicate whether the corresponding parameters are initialized or not
    bool chaindimInit_ = false;
    bool propcovInit_ = false;
    bool methodInit_ = false;
    bool outputInit_ false;
    
    bool gradflag_ = false; // Flag to check whether gradient information is given or not
    bool tensflag_ = false; // Flag that indicates whether tensor information is given or not
    
    bool fcnAcceptFlag_ = false; // Flag that indicates whether the accept function is given or not
    bool fcnRejectFlag_ = false; // Flag that indicates whether the reject function is given or not

    /// \brief Lower bounds
    Array1D<double> Lower_;
    /// \brief Upper bounds
    Array1D<double> Upper_;

    ///\brief Lower bound existence flags
    Array1D<int> lower_flag_;
    ///\brief Upper bound existence flags
    Array1D<int> upper_flag_;

};

//*****************************************

/// \class Hamiltonian Sampling
/// \brief Hamiltonian Markov Chain Monte Carlo class. Derived from the base class for MCMC
///        Implemented the HMCMC algorithms
class Hamilton:public MCMC{};

//*****************************************

/// \class MALA or Langevian Sampling
/// \brief MALA Markov Chain Monte Carlo class. Derived from the base class for Hamiltonian
///        Implemented the HMCMC algorithms, however L = 1 for the MALA class of algorithms
class MALA:public Hamilton{
public:
    // Initialization and set functions for private variables that are necessary to the MALA algorithms
    
    /// \brief Initialize epsilon for MALA
    void initEpsMALA(double eps_mala);
    
    
    // Get functions:
    /// \brief Get epsilon for MALA
    double getEpsMALA():
    
private:
    double eps_mala; // Epsilon for MALA algorithm
    
    // Flags to see if corresponding values are initalized or not
    bool epsMalaInit_ = false;
    
    double default_eps_mala_ = 0.1; // Default epsilon value for mala
};

//*****************************************

/// \class MMALA or Manifold variant of Langevian Sampling
/// \brief MMALA Markov Chain Monte Carlo class. Derived from the base class for MALA
///        Implemented the MMALA algorithms
class MMALA:public MALA{};

//*****************************************

/// \class Single-Site MCMC
/// \brief Single-Site Markov Chain Monte Carlo class. Derived from the base class for MCMC
///        Implemented the algorithms for single-site (Metropolis-within-Gibbs)
class SS:public MCMC{
    /// \brief Get function for number of sub steps
    int getNSubSteps();
private:
    int nSubSteps_ = this -> GetChainDim();
};

//*****************************************

/// \class Adaptive MCMC
/// \brief Adaptive Markov Chain Monte Carlo class. Derived from the base class for MCMC
///        Implemented the algorithms for aMCMC
class AMCMC:public MCMC{
public:
    // Initialization and set functions for private variables that are necessary to the aMCMC algorithms
    
    /// \brief Initialize adaptivity step parameters for aMCMC
    void initAdaptSteps(int adaptstart,int adaptstep, int adaptend);
    /// \brief Initialize the scaling factor gamma for aMCMC
    void initAMGamma(double gamma);
    /// \brief Initialize the covariance 'nugget' for aMCMC
    void initEpsCov(double eps_cov);
    
    // Get functions for any variables
    
    /// \brief Get function for the size=3 vector (t_start,t_step,t_end) that indicates when the adaptivity starts, how often the proposal covariance is updated and when the adaptivity ends, respectively
    void getAdaptSteps(Array1D<int> adaptstep_);
    /// \brief Get function for the coefficient behind the covariance scaling factor
    double getGamma();
    /// \brief Get function for the offset epsilon for Cholesky to be computationally feasible
    double getEpsCov();
    /// \brief Get function for number of sub steps
    int getNSubSteps();
    
    // Print functions:
    
    /// \brief Function to print the chain information on the screen
    void printChainSetup();

    
private:
    double gamma; // The coefficient behind the covariance scaling factor
    double eps_cov; // The offset epsilon for Cholesky to be computationally feasible
    int nSubSteps_ = 1;
    Array2D<double> curcov; // Covariance of the chain values so far
    Array1D<double> curmean; // Mean of the chain values sampled so far
    Array1D<int> adaptstep; // a size=3 vector (t_start,t_step,t_end) that indicates when the adaptivity starts, how often the proposal covariance is updated and when the adaptivity ends, respectively
    
    // Flags to indicate if corresponding values are initated or not
    bool adaptstepInit_ = false;
    bool gammaInit_ = false;
    bool epscovInit_ = false;
    
    //Default values for Gamma and EPS_Cov
    double default_gamma_ = 0.01;
    double default_eps_cov_ = 1e-8;
};

//*****************************************
/// \class Transitional MCMC
/// \brief Transitional Markov Chain Monte Carlo class. Derived from the base class for MCMC
///        Implemented the algorithms for TMCMC
class TMCMC:public MCMC{
public:
    
    // Initalization and set functions:
    
    /// \brief Set the default values for TMCMC function (Probably don't need this)
    void initDefaults();
    /// \brief Initialize the number of processes for TMCMC
    void initTMCMCNprocs(int tmcmc_nprocs);
    /// \brief Initialize the coefficient behind the covariance scaling
    ///        factor for TMCMC
    void initTMCMCGamma(double tmcmc_gamma);
    /// \brief Initialize the maximum allowed coefficient of variation for
    ///        the weights in TMCMC
    void initTMCMCCv(double tmcmc_cv);
    /// \brief Initialize the multiplicative factor for chain length to
    ///        encourage mixing in TMCMC
    void initTMCMCMFactor(int tmcmc_MFactor);
    /// \brief Initialize the choice to resample according to BASIS and
    ///        CATMIPs in TMCMC
    void initTMCMCBasis(bool tmcmc_basis);
    /// \brief Initialize the CATMIPs resampling parameter for TMCMC
    void initTMCMCCATSteps(int tmcmc_CATSteps);
    
    // Get functions:
    
    /// \brief Get the number of processes for TMCMC
    int getTMCMCNprocs();
    /// \brief Get the coefficient behind the covariance scaling
    ///        factor for TMCMC
    double getTMCMCGamma();
    /// \brief Get the maximum allowed coefficient of variation for
    ///        the weights in TMCMC
    double getTMCMCCv();
    /// \brief Initialize the multiplicative factor for chain length to
    ///        encourage mixing in TMCMC
    int getTMCMCMFactor();
    /// \brief Get the choice to resample according to BASIS and
    ///        CATMIPs in TMCMC
    bool getTMCMCBasis();
    /// \brief Get the CATMIPs resampling parameter for TMCMC
    int getTMCMCCATSteps();
    
private:
    int TMCMCNprocs; // The number of processes to use for parallel likelihood evaluation
    double TMCMCGamma; // The coefficient behind the covariance scaling factor
    double TMCMCCv; // The maximum allowed coefficient of variation for the weights
    int TMCMCMFactor; // The multiplicative factor for chain length to encourage mixing
    bool TMCMCBasis; // Resample according to BASIS and CATMIPs
    int TMCMCCATSteps; // CATMIPs resampling parameter
    std::vector<std::vector<double>> tmcmc_rngs; // The ranges for all samples
    
    // Initalization Flags
    bool tmcmcNprocsInit_ = false;
    bool tmcmcGammaInit_ = false;
    bool tmcmcCvInit_ = false;
    bool tmcmcMFactorInit_ = false;
    bool tmcmcBasisInit_ = false;
    bool tmcmcCATStepsInit_ = false;
    bool tmcmcRngsInit_ = false;
    
    // TMCMC Defaults
    int default_tmcmc_nprocs_ = 4;
    double default_tmcmc_gamma_ = -1.0;
    double default_tmcmc_cv_ = 0.1;
    int default_tmcmc_MFactor_ = 1;
    bool default_tmcmc_basis_ = true;
    int default_tmcmc_CATSteps_ = 1;
    std::vector<std::vector<double>> default_tmcmc_rngs_;
};

