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
#ifndef DFI_H_
#define DFI_H_

#include "Array1D.h"
#include "Array2D.h"
#include "Array3D.h"
#include "mcmc.h"
#include "quad.h"
#include "math.h"
#include "arrayio.h"

#include "dsfmt_add.h"
#include <sys/time.h>
#include "arraytools.h"
#include <string>
#include <iostream>
//for setting stdout precision
#include <iomanip>

#include "sampling.hpp"

//for UQTk KDE functionality
#include "tools.h"

//for PCE surrogate construction
#ifndef PCSET_H_SEEN
#define __wsu
#include "PCSet.h"
#endif //PCSET_H_SEEN

#include <stdlib.h>
#include <sstream>

//data container for surrogate model
class DFIsurr{

	public:

	//model surrogate containers
	
	//PCE dimension
	int PCEdim;

	//surrogate limits
	Array1D<double> surrLo;
	Array1D<double> surrHi;

        //surrogate defined flag
        bool surrDefined;

        //PCset defines and initializes polynomial chaos basis function set
        PCSet* surrModel;

	//number of PCE terms
        int numPCETerms;
	//cotainer for coeefficients
        Array2D<double> PCEcoefficients;
        //container for evaluated PCE basis functions
        Array2D<double> psiPCE;

	// surrogate evaulation function
	void evaluateSurr(Array1D<double> & modelOutput, Array1D<double> & params);
};




//data container
class dataPosteriorInformation{

	public:

	int seed;
	//dimension of the data space
	int dataDim;
	//parameter dimension
	int paramDim;
	//the number of constraints to impose
	int numConstraints;


	//boolean flag indicating that the data likelihood function is running for noise optimzation
	bool errorOpt;
	//boolean flag indicating that the the data chain has achieved burn-in
	bool dataChainBurnedIn;
	//counter for entry in data chain
	int dataChain_count;


	//container for passing noisy data sample to 1D noise optimization chain
	Array1D<double> optErrorParameters;

	//containers for data signal, i.e. y=f(x)
	Array1D<double> trueDatax;
	Array1D<double> trueDatay;
	//container for error signal (e.g. noise + bias)
	Array1D<double> error;
	//container for nominal parameters, which define the nominal true data (trueDatay)
	Array1D<double> nominalParameters;
	//container for nominal error parameters
	Array1D<double> nominalErrorParameters;
	//container for additional parameters needed to run the model (but not inferred)
	Array1D<double> hyperparameters;
	//container for parameters for surrogate data model
	Array1D<double> surrParameters;

	//ABC data likelihood containers
	//container to store labels describing the statistics to enforce as constraints
	vector<string> statLabels;
	//container to store the values of the statistics to enforce as constraints
	vector<double> statValues;
	//container to store the values of the ABC deltas for approximately enforcing statistics as contraints
	vector<double> statDeltas;


	//parameter chain write options
	int paramWriteFlag;
	//burn in chain file
	string burninParamWriteFile;
	//parameter chainfile
	string mainParamWriteFile;

	//length of parameter chains
	int parameterBurnInNumSamples;
	int parameterChainNumSamples;


	//surrogate model
	DFIsurr surrModelObj;
	int numSurr;
	vector<DFIsurr> surrModels;

};

//data container
class parameterPosteriorInformation{

	public:
	Array1D<double> dataChainState;
	//container for passing optimal error parameters to inner likelihood during error optimation
	Array1D<double> optErrorParams;
	Array1D<double> hyperparameters;
	bool errorOpt;
	//the input 'x' possibly required to compute the model function y=f(x)
	Array1D<double> trueDatax;

	// pointer to surrogate model object
	DFIsurr * surrModelObj_;

	vector<DFIsurr> * surrModels_;

};

//============================================
//user defined functions

//run the true data model
void userRunModel(Array1D<double> &modelDataY, Array1D<double> &modelDataX, Array1D<double> & parameters, Array1D<double> &hyperparameters);
//compute parameter (truth model and error model) parameter posterior
double userComputeParamLogPosterior(parameterPosteriorInformation * paramPostInfo, Array1D<double> parameters);
//compute parameter (truth model and error model) likelihood
double userComputeParamLogLikelihood(parameterPosteriorInformation * paramPostInfo, Array1D<double> modelDataOut, Array1D<double> parameters, Array1D<double> hyperparameters);
//compute statistics from inner parameter chain
void userComputeStatistics(Array1D<double> &parameterStatistics, Array1D<MCMC::chainstate> & parameterChain);
//define the target data signal
void userDefineData(dataPosteriorInformation & dataPostInfo);
//define the constraints
void userDefineConstraints(dataPosteriorInformation & dataPostInfo);
//specify the nominal values of the parameters
void userSpecifyNominalParams(dataPosteriorInformation & dataPostInfo);
//=============================================

//data log postrior function to be passed to UQTk MCMC object for data chain
double dataInferenceLogPosterior(Array1D<double>& m, void *info);
//parameter log posterior function to be passed to UQTk MCMC object for parameter chains
double parameterInferenceLogPosterior(Array1D<double>& beta, void *info);

//define inner parameter inference
void parameterInference(dataPosteriorInformation *dataPostInfo, Array1D<double> &m, Array1D<MCMC::chainstate>& parameterChainEntries);

//compute parameter (truth model and error model) parameter posterior
double computeParamLogPosterior(parameterPosteriorInformation * paramPostInfo, Array1D<double> parameters);
//compute parameter (truth model and error model) likelihood
double computeParamLogLikelihood(parameterPosteriorInformation * paramPostInfo, Array1D<double> modelDataOut, Array1D<double> parameters, Array1D<double> hyperparameters);
//compute statistics from inner parameter chain
void computeStatistics(Array1D<double> &parameterStatistics, Array1D<MCMC::chainstate> & parameterChain);



class DFI{

	private:
	//seed
	int seed;


	//metadata container to pass by argument to data MCMC chain
	dataPosteriorInformation dataPostInfo;
	//metadata container to pass by argument to parameter MCMC chain
	parameterPosteriorInformation paramPostInfo;

	//I/O ===============
	//log file name
	stringstream logFileName;
	//log file
	ofstream logFile;
	//===================

	//container for noisy data
	Array1D<double> noisyData;
	//data scale (defines an approximate scale for the noise to be added)
	double dataScale;


	//length of main data chain after burn-in
	int dataChainNumSamples;
	//length of burn-in chainlets
	int dataChainNumSamples_burnin;
	//length of error model parameter optimization chain
	double errorOptChainNumSamples;
	//target data chain acceptance ratio (a burn-in parameter)
	double targetDataChainAcceptanceRatio;
	//data chain acceptance ratio
	double dataChainAcceptanceRatio;
	//data chain mode ratio
	double dataPosteriorMode;

	//intial data chain (Gaussian) proposal covariance
        double dataChainPropCov_init;
	//scale for setting proposal covariance
        double dataChainPropCov_fac;


	//data chain proposal covariance matrix
	Array2D<double> dataChainPropCovMatrix;

	//=====wrappers for user specified functions=========
	//define the target data signal
	void defineData(dataPosteriorInformation & dataPostInfo){
		userDefineData(dataPostInfo);
	};
	//define the constraints
	void defineConstraints(dataPosteriorInformation & dataPostInfo){
		userDefineConstraints(dataPostInfo);
	};
	//specify the nominal values of the parameters
	void specifyNominalParams(dataPosteriorInformation & dataPostInfo){
		userSpecifyNominalParams(dataPostInfo);
	};
	//run the true data model
	void runModel(Array1D<double> &modelDataY, Array1D<double> &modelDataX, Array1D<double> & parameters, Array1D<double> &hyperparameters){

	//if a surrogate is defined, use it 
	if (dataPostInfo.surrModelObj.surrDefined){

		/* map parameters of interest to surrogate parameters */
		//userSurrMap(dataPostInfo.surrParameters, parameters, hyperparameters);

		/* check if parameters are within surrogate bounds */
		bool inBounds=true;
		for (int i=0; i<parameters.XSize();i++){
			if (  (parameters(i) < dataPostInfo.surrModelObj.surrLo(i)) || (parameters(i) >dataPostInfo.surrModelObj.surrHi(i)) ){
				inBounds=false;
				std::cout<<"Out of bounds: param("<<i+1<<"<<) = "<<parameters(i)<<std::endl;
			}
		}

		if (inBounds){
			dataPostInfo.surrModelObj.evaluateSurr(modelDataY, parameters);
		}else{
			userRunModel(modelDataY, modelDataX, parameters, hyperparameters);
		}

	}else{
		//run user specified detailed model
		userRunModel(modelDataY, modelDataX, parameters, hyperparameters);
	}

	};
	//===================================================


	public:
	//constructor
	DFI();
	//constructor
	DFI(string inputfile);
	//destructor
	~DFI();
	//perform data inference
	void dataInference();
	//perform data refitting
	void dataRefit();
	//build KDE densities
	void buildKDE(Array1D<int> KDEdim);
	//generate samples from a pdf
	void genSamples(Array2D<double> &pdf);
	//build surrogate model
	void buildSurrogateModel();
	//load surrogate model
	void loadSurrogateModel();
	//test surrogate model
	void testSurrogateModel();


};

#endif  //DFI_H_
