/* =====================================================================================

                      The UQ Toolkit (UQTk) version 3.1.2
                          Copyright (2022) NTESS
                        https://www.sandia.gov/UQToolkit/
                        https://github.com/sandialabs/UQTk

     Copyright 2022 National Technology & Engineering Solutions of Sandia, LLC (NTESS).
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
#include "dfi.h"
//user defined functions
#include "userFunctions.h"

//constructor
DFI::DFI(){

	//set seed for random number generator
	seed=13;

	std::cout<<"constructor..."<<std::endl;
	std::cout<<endl;

	//chain burn-in indicator
        dataPostInfo.dataChainBurnedIn=false;
        dataPostInfo.dataChain_count=0;

        ///set data dimension
        dataPostInfo.dataDim=10;
        noisyData.Resize(dataPostInfo.dataDim,0.0);
        ///set the seed
        dataPostInfo.seed=seed;

        //length of main data chain after burn-in
        dataChainNumSamples=10000;
        //length of burn-in chainlets
        dataChainNumSamples_burnin=500;
        //length of error model parameter optimization chain
        errorOptChainNumSamples=1000;
        //target data chain acceptance ratio (a burn-in parameter)
        targetDataChainAcceptanceRatio=0.20;
        dataChainAcceptanceRatio=0.0;
        dataPosteriorMode=0.0;

        //intial data chain proposal covariances
        dataChainPropCov_fac = 0.1;

        dataChainPropCovMatrix.Resize(dataPostInfo.dataDim,dataPostInfo.dataDim,0.0);

	//by default surrogate is not defined
	dataPostInfo.surrModelObj.surrDefined=false;


	//set lengths of parameter chains
        dataPostInfo.parameterBurnInNumSamples=10000;
        dataPostInfo.parameterChainNumSamples=50000;

        //write parameter chains to file (1=yes)
        dataPostInfo.paramWriteFlag=0;
        dataPostInfo.burninParamWriteFile="parameterChainBurnin.dat";
        dataPostInfo.mainParamWriteFile="parameterChain.dat";
}


//another constructor with user input file
DFI::DFI(string inputfile){

	std::cout<<"*constructor using user input file definitions* "<<inputfile<<std::endl;
	std::cout<<endl;
	std::cout<<"parsing inputs by keyword..."<<std::endl;


	ifstream DFIinput(inputfile.c_str());
	string line;
	double dvalue;
	int ivalue;
	string svalue;

	int linecount=0;
	while ( getline(DFIinput, line)  ){
		istringstream iss(line);
		//std::cout<<line<<std::endl;
		int readcount=0;
      		//parse the line
       		while (iss >>svalue ){
      			//keyword searches
			if (strcmp(svalue.c_str(),"seed")==0) {
				iss>>seed;
			        std::cout<<"seed= "<<seed<<std::endl;;
        		}
       			if (strcmp(svalue.c_str(),"data_chain_samples")==0) {
        			//length of main data chain after burn-in
				iss>>dataChainNumSamples;
			        std::cout<<"# of data chain samples to run= "<<dataChainNumSamples<<std::endl;;
			}
       			if (strcmp(svalue.c_str(),"data_chain_burnin_samples")==0) {
				//length of burn-in chainlets
				iss>>dataChainNumSamples_burnin;
			        std::cout<<"# of data chain burn-in samples to run= "<<dataChainNumSamples_burnin<<std::endl;;
			}
       			if (strcmp(svalue.c_str(),"parameter_chain_burnin_samples")==0) {
        			//length of main data chain after burn-in
				iss>>dataPostInfo.parameterBurnInNumSamples;
			        std::cout<<"# of parameter chain burn-in samples to run= "<<dataPostInfo.parameterBurnInNumSamples<<std::endl;;
			}
       			if (strcmp(svalue.c_str(),"parameter_chain_samples")==0) {
				//length of burn-in chainlets
				iss>>dataPostInfo.parameterChainNumSamples;
			        std::cout<<"# of parameter chain samples to run= "<<dataPostInfo.parameterChainNumSamples<<std::endl;;
			}
       			if (strcmp(svalue.c_str(),"error_optimization_chain_samples")==0) {
				//length of error model parameter optimization chain
				iss>>errorOptChainNumSamples;
			        std::cout<<"# of error optimization chain samples to run= "<<errorOptChainNumSamples<<std::endl;;
			}
       			if (strcmp(svalue.c_str(),"data_chain_propcov_fac")==0) {
        			//intial data chain proposal covariances scale factor
				iss>>dataChainPropCov_fac;
			        std::cout<<"data chain proposal covariance chain factor= "<<dataChainPropCov_fac<<std::endl;;
			}
       			if (strcmp(svalue.c_str(),"data_space_dimenion")==0) {
				///set data dimension
			        iss>>dataPostInfo.dataDim;
			        std::cout<<"dimensionality of data space= "<<dataPostInfo.dataDim<<std::endl;;
			}
       			if (strcmp(svalue.c_str(),"target_chain_accept")==0) {
				//target data chain acceptance ratio (a burn-in parameter)
			        iss>>targetDataChainAcceptanceRatio;
			        std::cout<<"target data chain acceptance ratio== "<<targetDataChainAcceptanceRatio<<std::endl;;
			}
       			if (strcmp(svalue.c_str(),"TEST")==0) {
				std::cout<<"TEST"<<std::endl;
			}
		}	
	}

	//structure passed to data likelihood function
        dataPostInfo.dataChainBurnedIn=false;
        dataPostInfo.dataChain_count=0;

        ///set data dimension
        //dataPostInfo.dataDim=10;
        noisyData.Resize(dataPostInfo.dataDim,0.0);
        ///set the seed
        dataPostInfo.seed=seed;

        //length of main data chain after burn-in
        //dataChainNumSamples=1000;
        //length of burn-in chainlets
        //dataChainNumSamples_burnin=100;
        //length of error model parameter optimization chain
        //errorOptChainNumSamples=1000;
        //target data chain acceptance ratio (a burn-in parameter)
        //targetDataChainAcceptanceRatio=0.20;
        dataChainAcceptanceRatio=0.0;
        dataPosteriorMode=0.0;

        //intial data chain proposal covariances
        //dataChainPropCov_fac = 0.1;
        //dataChainPropCov_init =pow(1.0e-7,2);

        dataChainPropCovMatrix.Resize(dataPostInfo.dataDim,dataPostInfo.dataDim,0.0);

	//by default surrogate is not defined
	dataPostInfo.surrModelObj.surrDefined=false;

	//set lengths of parameter chains
        //dataPostInfo.parameterBurnInNumSamples=10000;
        //dataPostInfo.parameterChainNumSamples=50000;

        //write parameter chains to file (1=yes)
        dataPostInfo.paramWriteFlag=0;
        dataPostInfo.burninParamWriteFile="parameterChainBurnin.dat";
        dataPostInfo.mainParamWriteFile="parameterChain.dat";

	//debug
	//exit(1);


}


DFI::~DFI(){



}


//perform data inference
void DFI::dataInference(){

	/*=====	setup problem ===================*/
        /*use user specified nominal model parameter values*/
        specifyNominalParams(dataPostInfo);
	/*set total parameter dimension (sum of inferred parameters: true model parameters and error model parameters)*/
	dataPostInfo.paramDim=dataPostInfo.nominalParameters.XSize()+dataPostInfo.nominalErrorParameters.XSize();
        /*use user defined constraints*/
        defineConstraints(dataPostInfo);
        /*total number of constraints*/
        dataPostInfo.numConstraints=dataPostInfo.statValues.size();

        /*use user defined initial noisy data chain state*/
        defineData(dataPostInfo);

        /*generate guess for nominal (smooth) true data*/
	runModel(dataPostInfo.trueDatay, dataPostInfo.trueDatax, dataPostInfo.nominalParameters, dataPostInfo.hyperparameters);

 	std::cout<<"Noisy initial data for data inference:"<<endl;
	/*define a scale for the added noise*/
	dataScale=dataPostInfo.trueDatay(0);
	//dataScale=1.0e-3;
        for (int i=0; i< dataPostInfo.dataDim; i++){
                /*initial data is the nominal data polluted with some added noise*/
                noisyData(i)=dataPostInfo.trueDatay(i) + dataScale*dataPostInfo.error(i);
                std::cout<<"x("<<i+1<<"): "<< dataPostInfo.trueDatax(i)<<", y("<<i+1<<"): "<<noisyData(i)<<endl;
        }
	/*=======================================*/


        //test that data posterior function is deterministic (call multiple times with constant data
	/*
        cout<<dataInferenceLogPosterior(initData, (void*) &dataPostInfo)<<endl;
        cout<<dataInferenceLogPosterior(initData, (void*) &dataPostInfo)<<endl;
        cout<<dataInferenceLogPosterior(initData, (void*) &dataPostInfo)<<endl;
        exit(1);
	*/

        //================ERROR PARAMETER OPTIMIZATION CHAIN ================================
        std::cout<<"============================================="<<std::endl;
        cout<<endl;
        cout<<"    STARTING OPTIMIZATION TO FIND FIRST GUESS FOR ERROR MODEL PARAMETERS"<<endl;
        cout<<endl;
        cout<<"============================================="<<endl;
        
	dataPostInfo.errorOpt=true;
        Array1D<double>optimalErrorParams=dataPostInfo.nominalErrorParameters;
        Array2D<double>errorOptChainPropCovMatrix(optimalErrorParams.XSize(),optimalErrorParams.XSize(),0.0);
        for (int i=0; i< optimalErrorParams.XSize(); i++){
                errorOptChainPropCovMatrix(i,i) = pow(optimalErrorParams(i)*0.1,2.0);
        }

        AMCMC errorOptimizationChain(dataInferenceLogPosterior, &dataPostInfo);
        errorOptimizationChain.setChainDim(optimalErrorParams.XSize());
        errorOptimizationChain.setWriteFlag(1);
        errorOptimizationChain.setSeed(seed);
        errorOptimizationChain.initAdaptSteps(errorOptChainNumSamples/2,100,errorOptChainNumSamples+1);  //( start,frequency,stop)
        errorOptimizationChain.initChainPropCov(errorOptChainPropCovMatrix);
        errorOptimizationChain.setOutputInfo("txt","optimalErrorParameter.chain",1,errorOptChainNumSamples);
        errorOptimizationChain.runChain(errorOptChainNumSamples,optimalErrorParams);
        double optErrorInfPosteriorMode = errorOptimizationChain.getMode(optimalErrorParams);

        cout<<"Optimal noise error parameters: ";
        for (int i=0; i<optimalErrorParams.XSize()-1;i++){
                cout<<optimalErrorParams(i)<<", ";
        }
        cout<<optimalErrorParams(optimalErrorParams.XSize()-1)<<endl;

        //write optimal error model parameters to file
        stringstream optErrParamFilename;
        optErrParamFilename <<"optimalErrorParams.dat";
        ofstream optErrParamsFile;
        optErrParamsFile.open(optErrParamFilename.str().c_str());
        for (int i=0;i<optimalErrorParams.XSize();i++){
                optErrParamsFile<< setprecision(15)<<optimalErrorParams(i)<<endl;
        }
        //

        //set the initial nominal error model parameters
        for (int i=0; i<optimalErrorParams.XSize();i++){
                dataPostInfo.nominalErrorParameters(i)=optimalErrorParams(i);
        }

        //set the initial noisy data state for the data chain using the optimal noise strength 
        for (int i=0; i< dataPostInfo.dataDim; i++){
                noisyData(i)= dataPostInfo.trueDatay(i) + exp(dataPostInfo.nominalErrorParameters(0))*dataPostInfo.error(i);
        }
        //=================================================================



        //================BURN-IN DATA CHAIN ================================
        cout<<"============================================="<<endl;
        cout<<endl;
        cout<<"    STARTING OUTER (DATA) BURN-IN CHAINLETS"<<endl;
        cout<<endl;
        cout<<"============================================="<<endl;

	dataPostInfo.errorOpt=false;

	//define initial MCMC proposal distributioni (iid-Gaussian) covariance matrix
        dataChainPropCov_init=noisyData(0)*dataChainPropCov_fac;
	for (int i=0; i< dataPostInfo.dataDim; i++){
                dataChainPropCovMatrix(i,i) = pow(dataChainPropCov_init,2.0);
        }

	int dataBurninCount=0;
        while (dataChainAcceptanceRatio< targetDataChainAcceptanceRatio){
                dataBurninCount++;
                cout<<"Starting data burn-in chain "<<dataBurninCount<<endl;
                AMCMC dataChain(dataInferenceLogPosterior,(void*) &dataPostInfo);
                dataChain.setChainDim(dataPostInfo.dataDim);
                dataChain.setWriteFlag(1);
                dataChain.setSeed(seed);
                dataChain.initAdaptSteps(dataChainNumSamples_burnin+1,1,dataChainNumSamples_burnin+1);
                dataChain.initChainPropCov(dataChainPropCovMatrix);
                dataChain.setOutputInfo("txt","dataBurnin.chain",1,dataChainNumSamples_burnin);
                dataChain.runChain(dataChainNumSamples_burnin,noisyData);

                dataPosteriorMode = dataChain.getMode(noisyData);

                cout << "Previous data chain (burn-in) acceptance ratio: " << dataChainAcceptanceRatio << endl;
                dataChain.getAcceptRatio(&dataChainAcceptanceRatio);
                cout << "New data chain (burn-in) acceptance ratio: " << dataChainAcceptanceRatio << endl;

                cout << "Previous data chain proposal covariance:  " << dataChainPropCov_init << endl;
                if (dataChainAcceptanceRatio<targetDataChainAcceptanceRatio){

                        if (dataChainPropCov_init>pow(exp(optimalErrorParams(0))/4.0 ,2.0) ){
                                ///redefine initial MCMC proposal distribution covariance matrix to adjust the acceptance ratio
                                dataChainPropCovMatrix.SetValue(0.0);
                                //adjust data proposal covariances
                                dataChainPropCov_init = dataChainPropCov_init*fmax(fmin(dataChainAcceptanceRatio/targetDataChainAcceptanceRatio,1.5),0.65);
                                for (int i=0; i< dataPostInfo.dataDim; i++){
                                        dataChainPropCovMatrix(i,i) = dataChainPropCov_init;
                                }
                        }else{
                                //we dont want to shrink the data proposal too small since then the data chain will not explore the space
				//increase acceptance by relaxing the ABC deltas
                                for (int i=0;i<dataPostInfo.numConstraints;i++){
                                        dataPostInfo.statDeltas[i]=0.9*dataPostInfo.statDeltas[i];
                                }
                        }

                }
                cout << "New data chain proposal covariance:  " << dataChainPropCov_init << endl;

	//end of burn-in loop
        }
	//=================================================================





	//================MAIN DATA CHAIN ================================
        cout<<endl;
        cout<<"============================================="<<endl;
        cout<<endl;
        cout<<"    STARTING OUTER (DATA) MAIN CHAIN"<<endl;
        cout<<endl;
        cout<<"============================================="<<endl;
        //

	dataPostInfo.dataChainBurnedIn=true;
        AMCMC dataChain(dataInferenceLogPosterior,(void*) &dataPostInfo);
        dataChain.setChainDim(dataPostInfo.dataDim);
        dataChain.setWriteFlag(1);
        dataChain.setSeed(seed);
        dataChain.initAdaptSteps(dataChainNumSamples+1,1,dataChainNumSamples+1);
        dataChain.initChainPropCov(dataChainPropCovMatrix);
        //output to file and screen
        dataChain.setOutputInfo("txt","dataMain.chain",1,dataChainNumSamples);
        dataChain.runChain(dataChainNumSamples,noisyData);
        //=================================================================



        //=== output final constraint values ===============
        cout<<"Final constraint ABC delta values:"<<endl;
        cout<<endl;
        for (int i=0;i<dataPostInfo.numConstraints;i++){
                cout<<dataPostInfo.statLabels[i]<<":"<<endl;
                cout<<"StatValue= "<<dataPostInfo.statValues[i]<<endl;
                cout<<"StatDelta= "<<dataPostInfo.statDeltas[i]<<endl;
                cout<<endl;
        }
	//==================================================

}


//perform data refitting
void DFI::dataRefit(){

        //read the number of columns in chain file
        ifstream countFileColumns("dataMain.chain");
        string firstLine;
        getline(countFileColumns,firstLine);
        stringstream parseLine(firstLine);
        string colstring;
        int numColumns=0;
        while (parseLine>>colstring){
                numColumns++;
        }
        cout<<"there are "<<numColumns<<" columns"<<endl;


        ///set data dimension (as read from data chain file)
        dataPostInfo.dataDim=numColumns-3;

        //userSpecifyNominalParams(&dataPostInfo);
        specifyNominalParams(dataPostInfo);
        dataPostInfo.paramDim=dataPostInfo.nominalParameters.XSize()+dataPostInfo.nominalErrorParameters.XSize();

        dataPostInfo.errorOpt=false;

        //write parameter chains to file (1=yes)
        //dataPostInfo.paramWriteFlag=0;

        //define the initial noisy data chain state
        defineData(dataPostInfo);

        //determine number of entries in data chain file
        ifstream countFileLines("dataMain.chain");
        string line;
        int dataChainLength=0;
        while ( getline(countFileLines,line) ){
                dataChainLength++;
 	}
        cout<<"reading data chain file containing "<<dataChainLength<<" entries"<<endl;

        //read data chain from file, (3 extra entries needed for chain entry number, alfa, log likelihood)
        Array2D<double> chain(dataChainLength,dataPostInfo.dataDim+3,0.0);
        read_datafile(chain,"dataMain.chain");


        //preform parameter inference
        Array1D<double> dataSample(dataPostInfo.dataDim,0.0);
        //parameter chain container
        Array1D<MCMC::chainstate> parameterChainEntries;

        stringstream chainoutname;
        stringstream chainsampsfilename;

        Array1D<double> lo(dataPostInfo.paramDim,0.0);
        Array1D<double> hi(dataPostInfo.paramDim,0.0);

        bool firstMaxMin=true;
        for (int i=0;i<dataChainLength;i++){
                for (int j=0;j<dataPostInfo.dataDim;j++){
                        dataSample(j)=chain(i,j+1);
                }
                std::cout<<"Refitting data chain sample "<<i<<" of "<<dataChainLength<<std::endl;
		chainoutname <<"parameterChainBurnin_"<<i<<".refit";
                dataPostInfo.burninParamWriteFile=chainoutname.str();;
                //clear the stringstream
                chainoutname.str("");
                chainoutname <<"parameterChain_"<<i<<".refit";
                dataPostInfo.mainParamWriteFile=chainoutname.str();
                chainoutname.str("");

		//need to pass reference to dataPostInfo, function expects a pointer
		parameterInference(&dataPostInfo, dataSample, parameterChainEntries);

		cout<<endl; //skip a line to make the log file more readable
                //write out chain samples and compute parameter limits
                chainsampsfilename.str("");
                chainsampsfilename <<"parameterChain_"<<i<<".refit";
                ofstream chainsampsfile;
                chainsampsfile.open(chainsampsfilename.str().c_str());
                int paramDim =parameterChainEntries(0).state.XSize();
                //loop through each chain entry
                for (int j=0;j<parameterChainEntries.XSize();j++){

                        for (int k=0;k<paramDim;k++){
                                chainsampsfile<< setprecision(15)<<parameterChainEntries(j).state(k)<<" ";

                                if (firstMaxMin){
                                        lo(k)=parameterChainEntries(j).state(k);
                                        hi(k)=parameterChainEntries(j).state(k);
                                }else{
                                        lo(k)=min(lo(k),parameterChainEntries(j).state(k));
                                        hi(k)=max(hi(k),parameterChainEntries(j).state(k));
                                }
                        }
                        firstMaxMin=false;
                        chainsampsfile<< endl;

                }


        }

	//kde limits filename (parameter ranes over which to build KDE) 
        stringstream limitsfilename;
        limitsfilename <<"kdeLimits.dat";
        //create file for write
        ofstream limitsfile;
        limitsfile.open(limitsfilename.str().c_str());
        //write parameter limits to file
        for (int i=0;i<parameterChainEntries(0).state.XSize()-1;i++){
                limitsfile <<setprecision(15)<< lo(i)<<" "<<hi(i)<<endl;
        }

        return;




}


void DFI::buildSurrogateModel(){


	//===== setup problem ===================
        //use user specified nominal model parameter values
        specifyNominalParams(dataPostInfo);
	//make sure data is defined
	defineData(dataPostInfo);
	//=======================================

	//PCE order 
	int PCEorder=2;
	//PCE dimension equal to truth model paramter dimension
	int PCEdim=dataPostInfo.nominalParameters.XSize();
	dataPostInfo.surrModelObj.PCEdim=PCEdim;

	//surrogate limits
	dataPostInfo.surrModelObj.surrHi.Resize(PCEdim,0.0);
	dataPostInfo.surrModelObj.surrLo.Resize(PCEdim,0.0);

	for (int i=0; i<PCEdim; i++){
		dataPostInfo.surrModelObj.surrLo(i) = 0.8*dataPostInfo.nominalParameters(i);
		dataPostInfo.surrModelObj.surrHi(i) = 1.2*dataPostInfo.nominalParameters(i);
		std::cout<<"lo= "<< dataPostInfo.surrModelObj.surrLo(i) <<", hi= "<< dataPostInfo.surrModelObj.surrHi(i)<<std::endl;
	}

	// write surrogate limits to file 
        stringstream surrLimitsFilename;
        surrLimitsFilename <<"limits_1.surr";
        //create file for write
        ofstream surrLimitsFile;
        surrLimitsFile.open(surrLimitsFilename.str().c_str());
        //write parameter limits to file
        for (int i=0;i<PCEdim;i++){
                surrLimitsFile <<"param"<<i+1<<" "<<setprecision(15)<< dataPostInfo.surrModelObj.surrLo(i)<<" "<<dataPostInfo.surrModelObj.surrHi(i)<<endl;
        }
	surrLimitsFile.close();

	// Legendre-uniform PCE
	dataPostInfo.surrModelObj.surrModel = new PCSet("NISPnoq",PCEorder,PCEdim,"LU");
	// Gauss-Hermite PCE
	//dataPostInfo.surrModelObj.surrModel = new PCSet("NISPnoq",PCEorder,PCEdim,"GH");
	dataPostInfo.surrModelObj.numPCETerms = dataPostInfo.surrModelObj.surrModel ->GetNumberPCTerms();

	std::cout<<"PCE dimension: "<<PCEdim<<std::endl;
	std::cout<<"PCE order: "<<PCEorder<<std::endl;
	std::cout<<"Number of PC terms: "<<dataPostInfo.surrModelObj.numPCETerms<<std::endl;

	int numQuadPts=dataPostInfo.surrModelObj.numPCETerms;
	//set the quadrature rule: grid-type, full/sparse,
	dataPostInfo.surrModelObj.surrModel->SetQuadRule("CC","full",numQuadPts);
	//container for quadrature points
	Array2D<double> quadPoints;
	dataPostInfo.surrModelObj.surrModel->GetQuadPoints(quadPoints);
	numQuadPts=quadPoints.XSize();

	//====== write out quadrature point ===================
	std::cout<<"Write out quadrature points"<<std::endl;

	for (int i=0;i<quadPoints.XSize();i++){
		for (int j=0;j<quadPoints.YSize();j++){
			std::cout<<"Quadpoint ("<<i<<","<<j<<") : "<<quadPoints(i,j)<<std::endl;
		}
	}
	//=====================================================



	//====== evaluate model at quadrature points ==================
	std::cout<<"Evaluate model at quadrature points"<<std::endl;

	Array2D<double> modelDataQuad(quadPoints.XSize(),dataPostInfo.dataDim,0.0);
	Array1D<double> modelDataOut(dataPostInfo.dataDim,0.0);
	Array1D<double> quadPoint(dataPostInfo.nominalParameters.XSize(),0.0);
	//for (int i=0; i<numQuadPts; i++){
	std::cout<<"Num. Quad Points= "<<numQuadPts<<std::endl;
	std::cout<<"Quad points (Xsize) "<<quadPoints.XSize()<<std::endl;
	std::cout<<"Quad points (Ysize) "<<quadPoints.YSize()<<std::endl;

	for (int i=0; i<quadPoints.XSize(); i++){
		std::cout<<"Quad point "<<i+1<<" (";


		for (int j=0; j<quadPoints.YSize();j++){
			quadPoint(j) = ((dataPostInfo.surrModelObj.surrLo(j)+dataPostInfo.surrModelObj.surrHi(j))/2.0) + quadPoints(i,j)*((dataPostInfo.surrModelObj.surrHi(j)-dataPostInfo.surrModelObj.surrLo(j))/2.0);
		}

		//run the model
		for (int j=0; j<quadPoints.YSize()-1;j++){
			std::cout<<quadPoint(j)<<", ";
		}
		std::cout<<quadPoint(quadPoints.YSize()-1)<<")"<<std::endl;

		runModel(modelDataOut, dataPostInfo.trueDatax, quadPoint, dataPostInfo.hyperparameters);

		//store model data evaluated at each quadrature point (i)
		for (int j=0; j<dataPostInfo.dataDim; j++){
			modelDataQuad(i,j) = modelDataOut(j);
		}

	}
	//==============================================================

	//====== write out data at quadrature points ===================
	std::cout<<"Write out model data at quadrature points:"<<std::endl;
	//loop over quadrature points
	for (int i=0; i<quadPoints.XSize(); i++){

		std::cout<<"(";
		for (int j=0; j<quadPoints.YSize()-1; j++){
			quadPoint(j) = ((dataPostInfo.surrModelObj.surrLo(j)+dataPostInfo.surrModelObj.surrHi(j))/2.0) + quadPoints(i,j)*((dataPostInfo.surrModelObj.surrHi(j)-dataPostInfo.surrModelObj.surrLo(j))/2.0);
			std::cout<<quadPoint(j)<<", ";
		}
		quadPoint(quadPoints.YSize()-1) = ((dataPostInfo.surrModelObj.surrLo(quadPoints.YSize()-1)+dataPostInfo.surrModelObj.surrHi(quadPoints.YSize()-1))/2.0) + quadPoints(i,quadPoints.YSize()-1)*((dataPostInfo.surrModelObj.surrHi(quadPoints.YSize()-1)-dataPostInfo.surrModelObj.surrLo(quadPoints.YSize()-1))/2.0);
		std::cout<<quadPoint(quadPoints.YSize()-1)<<") : ";

		//loop over data points
		for (int j=0; j<dataPostInfo.dataDim; j++){
			std::cout<<modelDataQuad(i,j)<<" ";
		}
		std::cout<<std::endl;
	}
	std::cout<<"Finished writing out model data at quadrature points"<<std::endl;
	//==============================================================

	Array1D<double> coeffVector(dataPostInfo.surrModelObj.numPCETerms,0.0); //vector of PC mode coefficients
	Array1D<double> dataAtQuadPt(numQuadPts,0.0); //vector of data values at each quad point

	dataPostInfo.surrModelObj.PCEcoefficients.Resize(dataPostInfo.dataDim,dataPostInfo.surrModelObj.numPCETerms,0.0);

	//===== determine PCE mode coefficient by Galerkin projection =======
	std::cout<<"Compute PCE coefficients..."<<std::endl;
	std::cout<<"Number of Quad points = "<<numQuadPts<<std::endl;
	//loop over each data point - time points (1 coefficient for each time point....)
	for (int i=0; i<dataPostInfo.dataDim; i++){

		//for (int j=0; j<dataPostInfo.surrModelObj.numPCETerms; j++){
		for (int j=0; j<numQuadPts; j++){
			 dataAtQuadPt(j) = modelDataQuad(j,i);
		}

		// Evaluate PC mode coefficients by Galerkin Projection
		dataPostInfo.surrModelObj.surrModel->GalerkProjection(dataAtQuadPt,coeffVector);

		for (int j=0; j<dataPostInfo.surrModelObj.numPCETerms; j++){
			dataPostInfo.surrModelObj.PCEcoefficients(i,j) =coeffVector(j);
                }

        }
	//====================================================================


	//=== write coefficients to file =====================================
	std::cout<<"Write PCE coefficients to file..."<<std::endl;
        stringstream PCEcoeffsFilename;
        PCEcoeffsFilename <<"PCEcoeffs_1.dat";
        //create file for write
        ofstream PCEcoeffsFile;
        PCEcoeffsFile.open(PCEcoeffsFilename.str().c_str());
        //write parameter limits to file
        for (int i=0;i<dataPostInfo.surrModelObj.PCEcoefficients.XSize();i++){
        	for (int j=0;j<dataPostInfo.surrModelObj.PCEcoefficients.YSize();j++){

                	PCEcoeffsFile <<setprecision(15)<< dataPostInfo.surrModelObj.PCEcoefficients(i,j)<<" ";
		}
		PCEcoeffsFile <<std::endl;
        }
	//====================================================================

	return;

}





void DFI::loadSurrogateModel(){

	//===== setup problem ===================
        //use user specified nominal model parameter values
        specifyNominalParams(dataPostInfo);
	//make sure data is defined
	defineData(dataPostInfo);
	//=======================================

	dataPostInfo.surrModels.resize(dataPostInfo.numSurr);

	for (int k=0;k<dataPostInfo.numSurr;k++){

	//==== check if surrogate coefficients file exists =========
	stringstream pcecoeffsFilename;
        pcecoeffsFilename <<"PCEcoeffs_"<<k+1<<".dat";
        ifstream infile(pcecoeffsFilename.str().c_str());
        std::cout<<std::endl;
	if (infile.good()){
                std::cout<<"PCE surrogate coefficients file ("<<pcecoeffsFilename.str()<<") found"<<std::endl;
        }else{
                std::cout<<"Error: PCE surrogate coefficients file ("<<pcecoeffsFilename.str()<<") not found"<<std::endl;
		return;
        }
	//==========================================================

	// read surrogate limits from file 
	stringstream surrLimitsFilename;
        surrLimitsFilename <<"limits_"<<k+1<<".surr";
	ifstream surrLimitsInput(surrLimitsFilename.str().c_str());
        string line;
        double dvalue;
        int ivalue;
        string svalue;
        vector<double> hi;
	vector<double> lo;
	// PCE dimension equals the number of lines in limits.surr 
	int PCEdim=0;
        while ( getline(surrLimitsInput, line)  ){
                istringstream iss(line);
        	std::cout<<line<<std::endl;
		int readcount=0;
		//parse the line
		while (iss >>svalue )   {
			std::cout<<svalue<<endl;
			iss>>dvalue;
			lo.push_back(dvalue);
			iss>>dvalue;
			hi.push_back(dvalue);
			PCEdim++;
		}
	}

	int PCEorder=2;

	dataPostInfo.surrModelObj.PCEdim=PCEdim;

	dataPostInfo.surrModelObj.surrModel = new PCSet("NISPnoq",PCEorder,PCEdim,"LU");
	dataPostInfo.surrModelObj.numPCETerms = dataPostInfo.surrModelObj.surrModel ->GetNumberPCTerms();

	// Find number of data points in surrogate coefficients file 
	int surrDim=0;
	while ( getline(infile, line)  ){
		surrDim++;
	}

	dataPostInfo.surrModelObj.PCEcoefficients.Resize(surrDim,dataPostInfo.surrModelObj.numPCETerms,0.0);
	read_datafile(dataPostInfo.surrModelObj.PCEcoefficients,pcecoeffsFilename.str().c_str());

	// store surrogate limits 
	dataPostInfo.surrModelObj.surrHi.Resize(PCEdim,0.0);
	dataPostInfo.surrModelObj.surrLo.Resize(PCEdim,0.0);
	for (int i=0; i<PCEdim; i++){
		dataPostInfo.surrModelObj.surrLo(i) = lo[i];
		dataPostInfo.surrModelObj.surrHi(i) = hi[i];
		std::cout<<"lo= "<< dataPostInfo.surrModelObj.surrLo(i) <<", hi= "<< dataPostInfo.surrModelObj.surrHi(i)<<std::endl;
	}

	//=== evaluate model at nominal parameter values=======
	Array1D<double> sampleParam(dataPostInfo.nominalParameters.XSize());
	for (int i=0;i<dataPostInfo.nominalParameters.XSize();i++){
		sampleParam(i)=dataPostInfo.nominalParameters(i);
	}
	Array1D<double> modelDataOut(dataPostInfo.dataDim,0.0);
	//run model
	userRunModel(modelDataOut, dataPostInfo.trueDatax, sampleParam, dataPostInfo.hyperparameters);

	std::cout<<"Detailed model output for nominal parameter values:"<<std::endl;
	for (int i=0; i<dataPostInfo.dataDim; i++){
		std::cout<<modelDataOut(i)<<" ";
	}
	std::cout<<std::endl;
	//=====================================================

	Array1D<double> surrModelDataOut(surrDim,0.0);
	Array2D<double> point(1,PCEdim);

	//map parameter to interval
	for (int i=0;i<PCEdim;i++){
		point(0,i) = (dataPostInfo.nominalParameters(i) - (dataPostInfo.surrModelObj.surrHi(i) +  dataPostInfo.surrModelObj.surrLo(i))/2.0 )/( (dataPostInfo.surrModelObj.surrHi(i) -  dataPostInfo.surrModelObj.surrLo(i))/2.0)  ;
	}

	std::cout<<"Surrogate model output for nominal parameter values:"<<std::endl;

	// sum over each data point (in "PCEcoeffs.dat" rows index data points, columns index PCE coefficients) 
	for (int i=0; i<surrDim; i++){

 		dataPostInfo.surrModelObj.surrModel->EvalBasisAtCustPts(point,dataPostInfo.surrModelObj.psiPCE);
		//sum over PCE terms
		for (int j=0;j<dataPostInfo.surrModelObj.numPCETerms;j++){
			surrModelDataOut(i)+=dataPostInfo.surrModelObj.PCEcoefficients(i,j)*dataPostInfo.surrModelObj.psiPCE(0,j);
		}
		std::cout<<surrModelDataOut(i)<<" ";

	}
	std::cout<<std::endl;
	//=====================================================

	//=== evaluate using member function==================

	std::cout<<"Surrogate model function:"<<std::endl;
	Array1D<double> temp_(dataPostInfo.nominalParameters.XSize(),0.0);
	for (int i=0;i<dataPostInfo.nominalParameters.XSize();i++){
		temp_(i)=dataPostInfo.nominalParameters(i);
	}

		dataPostInfo.surrModelObj.evaluateSurr(surrModelDataOut, temp_);

		for (int i=0; i<dataPostInfo.dataDim; i++){

			std::cout<<surrModelDataOut(i)<<" ";

		}
		std::cout<<std::endl;
	//====================================================


	//surrogate model is defined
	dataPostInfo.surrModelObj.surrDefined=true;

	}

	return;
}



void DFI::testSurrogateModel(){


	//check if surrogate is defined
	if (dataPostInfo.surrModelObj.surrDefined){

	double lo = dataPostInfo.surrModelObj.surrLo(0);
	double hi = dataPostInfo.surrModelObj.surrHi(0);

	int numDim=1;
	int numSamples=100;

	//run

	Array1D<double> modelDataOut(dataPostInfo.dataDim,0.0);
	Array1D<double> surrModelDataOut(dataPostInfo.dataDim,0.0);

	Array1D<double> sampleParam(dataPostInfo.nominalParameters.XSize(),0.0 );
	srand(seed);
	double surrError=0.0;
	for (int i=0; i<numSamples; i++){

		//generate a random sample for each parameter
		double rand_num = (double) rand() /RAND_MAX;

		std::cout<<std::endl;
		std::cout<<"Parameter sample "<<i+1<<": ";
		for (int i=0; i<sampleParam.XSize();i++){
			sampleParam(i)= dataPostInfo.surrModelObj.surrLo(i) +  rand_num*(dataPostInfo.surrModelObj.surrHi(i) -  dataPostInfo.surrModelObj.surrLo(i));
			std::cout<<sampleParam(i)<<" ";
		}
		std::cout<<std::endl;

        	//run model
        	userRunModel(modelDataOut, dataPostInfo.trueDatax, sampleParam, dataPostInfo.hyperparameters);
		std::cout<<"Detailed model"<<std::endl;
		for (int j=0; j<dataPostInfo.dataDim; j++){
			std::cout<<modelDataOut(j)<<" ";
		}
		std::cout<<std::endl;
		//evaluate surrogate
		dataPostInfo.surrModelObj.evaluateSurr(surrModelDataOut, sampleParam);
		std::cout<<"Surrogate model"<<std::endl;
		for (int j=0; j<dataPostInfo.dataDim; j++){
			std::cout<<surrModelDataOut(j)<<" ";
		}
		std::cout<<std::endl;

		for (int j=0; j<dataPostInfo.dataDim; j++){
			surrError+= pow(  (modelDataOut(j) - surrModelDataOut(j))/modelDataOut(j) ,2.0 );
		}

	}
	surrError = pow(surrError/(numSamples*dataPostInfo.dataDim) ,0.5);
	std::cout<<"RMS relative error= "<<surrError<<std::endl;

	}else{

		std::cout<<"Error: surrogate model is not defined"<<std::endl;

	}


	return;

}


//build KDE densities from samples
void DFI::buildKDE(Array1D<int> KDEdim){

	//first count how many parameter chain (refit) files are available
	stringstream refitChainFileName;
	bool refitFileExist=true;
	int numFiles=0;
	int numColumns=0;
	int numRows=0;
	//number of dimensions for KDE density
        int numDim=KDEdim.XSize();;

	//containers for KDE limits
	Array1D<double> lo(1,0.0);
        Array1D<double> hi(1,0.0);

	//flag for checkiung if we have already read one set of max/min from any file
        bool firstMaxMin=true;
	while(refitFileExist){
		refitChainFileName.str("");
		refitChainFileName<<"parameterChain_"<<numFiles<<".refit";
	        ifstream infile(refitChainFileName.str().c_str());
		if (infile.good()){
			numFiles++;
			std::cout<<refitChainFileName.str()<<" found, computing parameter bounds"<<std::endl;

			//find the number of parameters (number of columns) and refit chain length (number of rows) from the first discovered file
			if (numFiles==1){

				//determine number of entries in data chain file
				ifstream countFileLines(refitChainFileName.str().c_str());
				string line;
				while ( getline(countFileLines,line) ){
					numRows++;
				}

				//read the number of columns in chain file
				ifstream countFileColumns(refitChainFileName.str().c_str());
				string firstLine;
				getline(countFileColumns,firstLine);
				stringstream parseLine(firstLine);
				string colstring;
				while (parseLine>>colstring){
					numColumns++;
				}

				//estimate KDE limits from this first file
        			lo.Resize(numColumns,0.0);
        			hi.Resize(numColumns,0.0);

			}

			Array2D<double> readChainTemp(numRows,numColumns,0.0);
			read_datafile(readChainTemp,refitChainFileName.str().c_str());
			for (int j=0;j<numRows;j++){
				for (int k=0;k<numColumns;k++){
					if (firstMaxMin){
						lo(k)=readChainTemp(j,k);
						hi(k)=readChainTemp(j,k);
					}else{
						//we have read at least one set of values already, so compare with the current max/min
						lo(k)=min(lo(k),readChainTemp(j,k));
						hi(k)=max(hi(k),readChainTemp(j,k));
                                		}
                        		}
				//we have read at least one set of values already
                        	firstMaxMin=false;
			}




		}
		else{
			//break out
			refitFileExist=false;
		}
	}
	std::cout<<std::endl;
	std::cout<<numFiles<<" parameter chain (refit) files found, with "<<numColumns<<" parameters and "<<numRows<<" chain entries"<<std::endl;
	std::cout<<"The estimated KDE limits for the "<<numColumns<<" parameters are:"<<std::endl;
	for (int i=0;i<numColumns; i++){
		std::cout<<"min(param"<<i+1<<")="<<lo(i)<<", max(param"<<i+1<<")="<<hi(i)<<std::endl;
	}
	std::cout<<std::endl;


	//======= build KDE posteriors ====================
	int n_kde=64;
	//read in each file
	//container for reading samples from file
        Array2D<double> refitChain(numRows,numColumns,0.0);
	//container for samples from select dimensions for constructing marginal KDE
        Array2D<double> refitChainMarg(numRows,numDim,0.0);
	//
	//Array2D<double> grid(n_kde,numColumns);
	Array2D<double> grid(n_kde,numDim);
	Array2D<double> points;
	//define grid
	for (int j=0;j<numDim;j++){
		for (int i=0;i<n_kde;i++){
			grid(i,j)=lo(KDEdim(j)-1) + i*( hi(KDEdim(j)-1)-lo(KDEdim(j)-1) )/(n_kde-1);
		}
	}


	//generate KDE points array
	generate_multigrid(points,grid);

	int totpts = (int) points.XSize();
	Array1D<double> pooleddensKDE(totpts,0.0);
	stringstream KDEFileName;
	for (int i=0; i<numFiles; i++){
		refitChainFileName.str("");
		refitChainFileName<<"parameterChain_"<<i<<".refit";
		std::cout<<"reading "<<refitChainFileName.str()<<" and computing KDE density"<<std::endl;
		read_datafile(refitChain,refitChainFileName.str().c_str());

		//just copy what is needed
		for (int j=0; j<numRows; j++){
			for (int k=0; k<numDim; k++){
				refitChainMarg(j,k)=refitChain(j,KDEdim(k)-1);
			}
		}

		//compute KDE
		Array1D<double> densKDE(totpts,0.0);
		getPdf_cl(refitChainMarg, points, densKDE, 1, 1.0);


		for (int j=0;j<totpts;j++){
			pooleddensKDE(j)=pooleddensKDE(j)+densKDE(j);
		}

		//===write pooled KDE density to file==========
        	KDEFileName.str("");
		KDEFileName <<"KDE_"<<i<<".dat";
	        //create file for write
	        ofstream KDEFile;
        	KDEFile.open(KDEFileName.str().c_str());
	        //write parameter limits to file
	        for (int j=0;j<totpts;j++){
			for (int k=0;k<numDim;k++){
                		KDEFile <<setprecision(15)<< points(j,k)<<" ";
			}
			KDEFile<<setprecision(15)<<densKDE(j)<<std::endl;
        	}
		//==============================================

	}

	//average
	for (int i=0;i<totpts;i++){
		pooleddensKDE(i)=pooleddensKDE(i)/( (double) numFiles);
	}

	//===write pooled KDE density to file==========
        stringstream pooledKDEFileName;
        pooledKDEFileName <<"pooledKDE.dat";
        //create file for write
        ofstream pooledKDEFile;
        pooledKDEFile.open(pooledKDEFileName.str().c_str());
        //write parameter limits to file
        for (int i=0;i<totpts;i++){
		//for (int j=0;j<numColumns;j++){
		for (int j=0;j<numDim;j++){
                	pooledKDEFile <<setprecision(15)<< points(i,j)<<" ";
		}
		pooledKDEFile<<setprecision(15)<<pooleddensKDE(i)<<std::endl;
        }
	//=============================================
	//generate samples using rejection sampling


	//=================================================

	return;
}

void DFI::genSamples(Array2D<double> & pdf){

	std::cout<<"generating samples from the pdf and writing to file"<<std::endl;

	return;
}

/*
//compute model evidence
void DFI::computeEvidence(){

}
*/

void DFIsurr::evaluateSurr(Array1D<double> & modelOutput, Array1D<double> & params){

	/* DEBUG
	std::cout<<"PCE coeffs: ==================="<<std::endl;
        for (int i=0; i<PCEcoefficients.XSize(); i++){
        	for (int j=0; j<PCEcoefficients.YSize(); j++){

                	std::cout<<PCEcoefficients(i,j)<<" ";
		}
                std::cout<<std::endl;
	}
	std::cout<<"PCE coeffs: ==================="<<std::endl;
	*/

	Array2D<double>point (1,PCEdim,0.0);
	for (int i=0; i<PCEdim; i++){
		point (0,i) = (params(i) - (surrHi(i) + surrLo(i))/2.0 )/( (surrHi(i) - surrLo(i)/2.0) )    ;
	}
	surrModel->EvalBasisAtCustPts(point,psiPCE);
	/* DEBUG
	std::cout<<"psiPCE: ==================="<<std::endl;
        for (int i=0; i<psiPCE.XSize(); i++){
        	for (int j=0; j<psiPCE.YSize(); j++){

                	std::cout<<psiPCE(i,j)<<" ";
		}
                std::cout<<std::endl;
	}
	std::cout<<"psiPCE: ==================="<<std::endl;
	*/
	//loop over each data point (rows)
        for (int i=0; i<PCEcoefficients.XSize(); i++){
		//sum over PCE terms
		modelOutput(i)=0.0;
		for (int j=0;j<numPCETerms;j++){
                        modelOutput(i)+=PCEcoefficients(i,j)*psiPCE(0,j);
                }

        }

	return;
}

double dataInferenceLogPosterior(Array1D<double>& m, void *info){

        cout<<endl;
        cout<<"Computing data (log) likelihood"<<endl;

        //cast void pointer to appropriate type
        dataPosteriorInformation * dataPostInfo = (dataPosteriorInformation *) info;

        //parameter chain container
        Array1D<MCMC::chainstate> parameterChainEntries;

        //for a given data set (m), compute the parameter chain (parameterChainEntries)
        parameterInference(dataPostInfo, m, parameterChainEntries);

        //if data chain has burned in, write parameter samples to file
        if (dataPostInfo->dataChainBurnedIn){

                stringstream paramChainFilePath;
                paramChainFilePath<<"./paramChains/param_"<<dataPostInfo->dataChain_count<<".chain";
                ofstream paramChainFile;
                paramChainFile.open(paramChainFilePath.str().c_str());
                for (int i=0;i<parameterChainEntries.XSize();i++){
                        for (int j=0;j<parameterChainEntries(0).state.XSize();j++){
                                paramChainFile<<setprecision(15)<<parameterChainEntries(i).state(j)<<" ";
                        }
                        paramChainFile<< setprecision(15)<<parameterChainEntries(i).alfa<<" "<<parameterChainEntries(i).post <<endl;
                }

                stringstream dataProposalFilePath;
                dataProposalFilePath<<"./paramChains/data_"<<dataPostInfo->dataChain_count<<".sample";
                ofstream dataProposalFile;
                dataProposalFile.open(dataProposalFilePath.str().c_str());
                for (int j=0;j<m.XSize();j++){
                        dataProposalFile<<setprecision(15)<<m(j)<<" ";
                }
                dataProposalFile<<endl;

         	dataPostInfo->dataChain_count++;
	}

	 //compute statistics using the parameter chain
        Array1D<double> parameterStatistics(dataPostInfo->numConstraints,0.0);
	cout<<"	Calculating parameter statistics from parameter MCMC chain"<<endl;
	userComputeStatistics( parameterStatistics, parameterChainEntries);
        for(int j=0; j<parameterStatistics.XSize(); j++){
		cout<<"	stat "<<j+1<<" ("<<dataPostInfo->statLabels[j]<<", Delta="<<dataPostInfo->statDeltas[j]<<"): "<<parameterStatistics(j)<<", target: "<<dataPostInfo->statValues[j]<<endl;
        }

        //compute approximate Bayesian likelhood (ABC)
        double dataInferenceLogLikelihood=0.0;
        Array1D<double> deltas(dataPostInfo->numConstraints,0.0);

        //compute ABC kernels
        for (int i=0; i<dataPostInfo->numConstraints; i++){
                deltas(i)=dataPostInfo->statDeltas[i]*pow( (dataPostInfo->statValues[i] - parameterStatistics(i))/dataPostInfo->statValues[i] , 2.0);
        }

        for(int i=0; i<dataPostInfo->numConstraints; i++){
                dataInferenceLogLikelihood -= deltas(i);
        }

        for(int i=0; i<dataPostInfo->numConstraints-1; i++){
		cout<<"	d"<<i+1<<": "<<deltas(i)<<", ";
        }
        cout<<"d"<<dataPostInfo->numConstraints<<": "<<deltas(dataPostInfo->numConstraints-1)<<endl;

	cout<<"	-(";
        for(int i=0; i<dataPostInfo->numConstraints-1; i++){
                cout <<"d"<<i+1<<" + ";
        }
        cout<<"d"<<dataPostInfo->numConstraints<<"): "<< dataInferenceLogLikelihood<<endl;
	cout<<"	***************************"<<endl;

        //product of likelihood and prior
        double dataInfLogPost=dataInferenceLogLikelihood * 1.0;

        cout <<"        Exiting with data log posterior function = "<<dataInfLogPost<<endl;

        return dataInfLogPost;

}


void parameterInference(dataPosteriorInformation *dataPostInfo, Array1D<double> &m, Array1D<MCMC::chainstate> & parameterChainEntries){
	cout<<"	****parameter burn-in******"<<endl;

        parameterPosteriorInformation paramPostInfo;
	//declare pointer
	paramPostInfo.surrModelObj_=&dataPostInfo->surrModelObj;
	paramPostInfo.surrModels_=&dataPostInfo->surrModels;

        //retrieve total parameter dimension dataPostInfo container
        int parameterDimension=dataPostInfo->nominalParameters.XSize() + dataPostInfo->nominalErrorParameters.XSize();

        //if this is a noise strength optimzation case, the state 'm' is the noise strength sample
        if (dataPostInfo->errorOpt){
                //running noise optimization, error parameters are not inferred this parameter chain
                parameterDimension=parameterDimension -dataPostInfo->nominalErrorParameters.XSize();
                //set the optimal noise strength (log)
                paramPostInfo.optErrorParams=m;
                paramPostInfo.errorOpt=true;
                //resize the data container that will be passes to the parameter likelihood function
                paramPostInfo.dataChainState.Resize(dataPostInfo->dataDim,0.0);
                //construct the noisy signal with noise scaled by the proposed noise strength
                for (int i=0;i<dataPostInfo->dataDim; i++){
                        paramPostInfo.dataChainState(i)=dataPostInfo->trueDatay(i) + exp(paramPostInfo.optErrorParams(0))*dataPostInfo->error(i);
                }
                paramPostInfo.trueDatax=dataPostInfo->trueDatax;
                paramPostInfo.hyperparameters=dataPostInfo->hyperparameters;
        }else{
                //if this is a regular data likelihood case, the state 'm' is the noisy data sample
                paramPostInfo.errorOpt=false;
                //load the data chain proposed state
                paramPostInfo.dataChainState=m;
                paramPostInfo.hyperparameters=dataPostInfo->hyperparameters;
                //find the input 'x' to the model function y=f(x)
                paramPostInfo.trueDatax=dataPostInfo->trueDatax;
                paramPostInfo.hyperparameters=dataPostInfo->hyperparameters;
        }

	Array1D<double> initParameter(parameterDimension,0.0);
        for (int i=0;i<dataPostInfo->nominalParameters.XSize();i++){
                initParameter(i)=dataPostInfo->nominalParameters(i);
        }

        if (!dataPostInfo->errorOpt){
                for (int i=0;i<dataPostInfo->nominalErrorParameters.XSize();i++){
                        initParameter(i+dataPostInfo->nominalParameters.XSize())=dataPostInfo->nominalErrorParameters(i);
                }
        }


	cout<<"	Initial parameter state: ";
        for (int i=0;i<parameterDimension-1;i++){
                cout<<i+1<<": "<<initParameter(i)<<", ";
        }
        cout<<parameterDimension<<": "<<initParameter(parameterDimension-1)<<endl;

        double parameterChainPropCov_init = 0.01;
        double parameterPosteriorMode=0.0;
        double parameterPosteriorPreviousMode=0.0;
        double parameterPosteriorRelativeChange=1.0;
        double parameterChainAcceptanceRatio=0.0;

        //define initial MCMC proposal distribution covariance matrix
        Array2D<double> parameterChainPropCovMatrix(parameterDimension,parameterDimension,0.0);
        for (int i=0; i< parameterDimension; i++){
                parameterChainPropCovMatrix(i,i) = pow(parameterChainPropCov_init*initParameter(i),2.0);
        }

        //================BURN-IN PARAMETER CHAIN ================================
        //counter for number of burn-in chains run
	int burnInChainCount=0;
        while (parameterPosteriorRelativeChange > 0.01){
                burnInChainCount++;
		cout<<"	Parameter burn-in chainlet: "<<burnInChainCount<<endl;
                AMCMC parameterChain(parameterInferenceLogPosterior,(void*) &paramPostInfo);
                parameterChain.setChainDim(parameterDimension);
                parameterChain.setSeed(dataPostInfo->seed);
                parameterChain.initAdaptSteps(dataPostInfo->parameterBurnInNumSamples/2,500,dataPostInfo->parameterBurnInNumSamples+1);  //( start,frequency,stop)
                parameterChain.initChainPropCov(parameterChainPropCovMatrix);
                parameterChain.setWriteFlag(dataPostInfo->paramWriteFlag);
                parameterChain.setOutputInfo("txt",dataPostInfo->burninParamWriteFile,1,dataPostInfo->parameterBurnInNumSamples);

                parameterChain.runChain(dataPostInfo->parameterBurnInNumSamples,initParameter);
                //MCMC acceptance ratio
                parameterChain.getAcceptRatio(&parameterChainAcceptanceRatio);
                //chain mode
                parameterPosteriorMode = parameterChain.getMode(initParameter);

                cout <<"        Parameter burn-in chain Acceptance Ratio : " << parameterChainAcceptanceRatio << endl;
                cout <<"        Parameter burn-in chain mode is: "<<parameterPosteriorMode<<endl;
                cout <<"        MAP: ";
                for (int i=0; i<parameterDimension; i++){
                        cout <<initParameter(i)<<" ";
                }
                cout <<endl;

                parameterPosteriorRelativeChange = fabs(parameterPosteriorPreviousMode - parameterPosteriorMode)/parameterPosteriorMode;
                cout <<"        Relative change in parameter burn-in chain mode = "<<parameterPosteriorRelativeChange<<endl;
                parameterPosteriorPreviousMode = parameterPosteriorMode;
        }
	cout<<"	Burn-in chain parameter posterior mode relative change < 0.01, proceeding to main parameter chain"<<endl;
	cout<<"	****parameter main*********"<<endl;
        //============================================================================


        //================MAIN PARAMETER CHAIN ================================
        AMCMC parameterChain(parameterInferenceLogPosterior,(void*) &paramPostInfo);
        parameterChain.setChainDim(parameterDimension);
        parameterChain.setSeed(dataPostInfo->seed);
        parameterChain.initAdaptSteps(dataPostInfo->parameterChainNumSamples/2,500,dataPostInfo->parameterChainNumSamples+1);  //( start,frequency,stop)
        parameterChain.initChainPropCov(parameterChainPropCovMatrix);
        parameterChain.setWriteFlag(dataPostInfo->paramWriteFlag);
        parameterChain.setOutputInfo("txt",dataPostInfo->mainParamWriteFile,1,dataPostInfo->parameterBurnInNumSamples);
        parameterChain.runChain(dataPostInfo->parameterChainNumSamples,initParameter);
        //MCMC acceptance ratio
        parameterChain.getAcceptRatio(&parameterChainAcceptanceRatio);
        //chain mode
        parameterPosteriorMode = parameterChain.getMode(initParameter);
        cout <<"        Parameter chain Acceptance Ratio : " << parameterChainAcceptanceRatio << endl;
        cout <<"        Chain mode is: "<<parameterPosteriorMode<<endl;
        cout <<"        MAP: ";
        for (int i=0; i<parameterDimension; i++){
                cout <<initParameter(i)<<" ";
        }
        cout <<endl;
        //get the parameter chain
        parameterChain.getFullChain(parameterChainEntries);
	//============================================================================



        //===== output select main parameter chain samples ==========
        cout <<"        Select samples from the parameter (inner) main chain"<<endl;
        cout <<"        0: ";
        for (int j=0; j<parameterDimension; j++){
                cout<<parameterChainEntries(0).state(j)<<" ";
        }
        cout<<endl;
        cout <<"        1: ";
        for (int j=0; j<parameterDimension; j++){
                cout<<parameterChainEntries(1).state(j)<<" ";
        }
        cout<<endl;
        for (int i=dataPostInfo->parameterChainNumSamples/10; i<=dataPostInfo->parameterChainNumSamples; i=i + dataPostInfo->parameterChainNumSamples/10){
		cout<<"	"<<i<<": ";
                for (int j=0; j<parameterDimension; j++){
                        cout<<parameterChainEntries(i).state(j)<<" ";
                }
                cout<<endl;
        }
	//============================================================

        return;
}


double parameterInferenceLogPosterior(Array1D<double>& parameters, void *info){
        //this is a wrapper function that calls the user defined function

	//cast void pointer to appropriate type
        parameterPosteriorInformation * paramPostInfo = (parameterPosteriorInformation *) info;


        double parameterLogPosterior=0.0;
        parameterLogPosterior=computeParamLogPosterior( paramPostInfo, parameters);

        return parameterLogPosterior;
}


//compute parameter (truth model and error model) parameter posterior
double computeParamLogPosterior(parameterPosteriorInformation * paramPostInfo, Array1D<double> parameters){
        return userComputeParamLogPosterior(paramPostInfo, parameters);
};
//compute parameter (truth model and error model) likelihood
double computeParamLogLikelihood(parameterPosteriorInformation * paramPostInfo, Array1D<double> modelDataOut, Array1D<double> parameters, Array1D<double> hyperparameters){
        return userComputeParamLogLikelihood(paramPostInfo, modelDataOut, parameters, hyperparameters);
};
//compute statistics from inner parameter chain
void computeStatistics(Array1D<double> &parameterStatistics, Array1D<MCMC::chainstate> & parameterChain){
        return userComputeStatistics(parameterStatistics, parameterChain);
};
