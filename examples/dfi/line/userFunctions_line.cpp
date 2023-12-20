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

     Questions? Contact the UQTk Developers at https://github.com/sandialabs/UQTk/discussions
     Sandia National Laboratories, Livermore, CA, USA
===================================================================================== */
#include "dfi.h"

#ifndef USER_FUNCS_H_
#define USER_FUNCS_H_

void userSpecifyNominalParams(dataPosteriorInformation & dataPostInfo){
	//=========================================================
	//===== USER DEFINED ======================================
	//=========================================================

        //define the initial nominal model parameters
        dataPostInfo.nominalParameters.Resize(2,0.0);
        //rate constant
        dataPostInfo.nominalParameters(0)=1.0;
        dataPostInfo.nominalParameters(1)=0.2;

	//define the initial error model parameters
        dataPostInfo.nominalErrorParameters.Resize(1,0.0);

	//if optimal estimate exists, use it
	ifstream infile("optimalErrorParams.dat");
        if (infile.good()){
                infile>>dataPostInfo.nominalErrorParameters(0);
        }else{
                dataPostInfo.nominalErrorParameters(0)= -2.0;
        }

	// total number of surrogate models used (e.g. >1 if using a multi-surrogate tiling of parameter space)
        dataPostInfo.numSurr=1;

	//=========================================================
	return;
}

void userDefineConstraints(dataPosteriorInformation & dataPostInfo){
	//=========================================================
	//===== USER DEFINED ======================================
	//=========================================================

	//place information (label, value, ABC delta) for each constraint here
        //1
        dataPostInfo.statLabels.push_back("mean of line slope");
        dataPostInfo.statValues.push_back(dataPostInfo.nominalParameters(0) );
        dataPostInfo.statDeltas.push_back(1000.0);
        //2
        dataPostInfo.statLabels.push_back("std of line slope");
        dataPostInfo.statValues.push_back(0.1);
        dataPostInfo.statDeltas.push_back(100.0);
        //3
        dataPostInfo.statLabels.push_back("mean of intercept");
        dataPostInfo.statValues.push_back(dataPostInfo.nominalParameters(1));
        dataPostInfo.statDeltas.push_back(1000.0);
        //4
        dataPostInfo.statLabels.push_back("std of intercept");
        dataPostInfo.statValues.push_back(0.1);
        dataPostInfo.statDeltas.push_back(100.0);



	//=========================================================
	return;
}

void userDefineData(dataPosteriorInformation &dataPostInfo){
	//=========================================================
	//===== USER DEFINED ======================================
	//=========================================================

        /*resize containters to data dimension*/
        dataPostInfo.trueDatax.Resize(dataPostInfo.dataDim,0.0);
        dataPostInfo.trueDatay.Resize(dataPostInfo.dataDim,0.0);
        //extra containers if needed ...
	dataPostInfo.error.Resize(dataPostInfo.dataDim,0.0);

        //generate initial guess of noise signal
        dsfmt_gv_init_gen_rand(dataPostInfo.seed);
        for (int i=0; i< dataPostInfo.dataDim; i++){
                dataPostInfo.error(i)=dsfmt_gv_genrand_nrv(); //normally distributed random variable;
                //dataPostInfo.noise(i) = dsfmt_gv_genrand_urv(); //uniformly distributed random variable
        }

        //define the data input 'x', i.e. y=f(x)
        //====================
        cout <<"Determine data input space 'x', i.e. y=f(x)"<<endl;
        cout <<"Data dimension =  "<<dataPostInfo.dataDim<<endl;
        double x_min=0.0;
        double x_max=1.0;
        for (int i=0; i<dataPostInfo.dataDim; i++){
              dataPostInfo.trueDatax(i)=x_min + (x_max-x_min)*i/(dataPostInfo.dataDim-1.0);
              cout<<"x("<<i+1<<"): "<< dataPostInfo.trueDatax(i)<<endl;
      }
	// write the input space 'x' to file *
        stringstream inputspaceFilename;
        inputspaceFilename <<"modelinput_x.dat";
        ofstream inputspaceFile;
        inputspaceFile.open(inputspaceFilename.str().c_str());
        for (int i=0;i<dataPostInfo.dataDim;i++){
                inputspaceFile<< setprecision(15)<<dataPostInfo.trueDatax(i)<<" ";;
        }
        inputspaceFile.close();

        //====================
        //define the data output 'y', i.e. y=f(x)
        cout<<"Model output data (no error) using nominal parameter values:"<<endl;
        userRunModel(dataPostInfo.trueDatay, dataPostInfo.trueDatax, dataPostInfo.nominalParameters, dataPostInfo.hyperparameters);
        for (int i=0; i<dataPostInfo.dataDim; i++){
                cout<<"x("<<i+1<<"): "<< dataPostInfo.trueDatax(i)<<", y("<<i+1<<"): "<<dataPostInfo.trueDatay(i)<<endl;
        }


	//=========================================================
        return;
}

void userRunModel(Array1D<double> &modelDatay, Array1D<double> & modelDatax, Array1D<double> &parameters, Array1D<double> &hyperparameters){
	//=========================================================
	//===== USER DEFINED ======================================
	//=========================================================

        for (int i=0; i<modelDatay.XSize(); i++){
                //y=f(x;p) where 'p' are some parameters
                //line
                modelDatay(i) = parameters(0)*modelDatax(i) + parameters(1);
        }
        //==================================================================================

        return;
}



double userComputeParamLogPosterior(parameterPosteriorInformation * paramPostInfo, Array1D<double> parameters){
	double parameterLogPosterior=0.0;
	//=========================================================
	//===== USER DEFINED ======================================
	//=========================================================

        int numDataPoints= paramPostInfo->dataChainState.XSize();
        Array1D<double> modelDataOut(numDataPoints,0.0);

	// if a surrogate is defined use it 
        if (paramPostInfo->surrModelObj_->surrDefined){
                // check if parameters are within surrogate bounds
                bool inBounds=true;
                for (int i=0; i<parameters.XSize();i++){
                        if (  (parameters(i) < paramPostInfo->surrModelObj_->surrLo(i)) || (parameters(i) >paramPostInfo->surrModelObj_->surrHi(i)) ){
                                inBounds=false;
                        }
                }

                if (inBounds){
                        paramPostInfo->surrModelObj_->evaluateSurr(modelDataOut, parameters);
                }else{
                        userRunModel(modelDataOut, paramPostInfo->trueDatax, parameters, paramPostInfo->hyperparameters);
                }

                paramPostInfo->surrModelObj_->evaluateSurr(modelDataOut, parameters);
        }else{
                //run the full model
                userRunModel(modelDataOut, paramPostInfo->trueDatax, parameters, paramPostInfo->hyperparameters);
        }

	parameterLogPosterior=userComputeParamLogLikelihood(paramPostInfo,modelDataOut,parameters,paramPostInfo->hyperparameters);
        //==============================================================

        return parameterLogPosterior;
}







double userComputeParamLogLikelihood(parameterPosteriorInformation * paramPostInfo, Array1D<double> modelDataOut, Array1D<double> parameters, Array1D<double> hyperparameters){
	double paramLogLikelihood=0.0;
	//=========================================================
	//===== USER DEFINED ======================================
	//=========================================================

        //container for noise strength variable
        double sig=0.0;

        if (paramPostInfo->errorOpt){
                sig=exp(paramPostInfo->optErrorParams(0));
        }else{
                sig=exp( parameters(parameters.XSize()-1) );
        }

        //compute likelihood function
        double modelError=0.0;
        double sumError=0.0;
        for (int i=0; i < modelDataOut.XSize(); i++){
                //compute error pointwise
                modelError  = paramPostInfo->dataChainState(i)  - modelDataOut(i);  // "DATA BEING TESTED" - "DATA GENERATED FROM INFERENCE"
                //likelihood is the product of pointwise liklihoods (sum of logs)
                sumError += log(sig) + pow(modelError,2.0)/( 2.0*(sig*sig) ); // Error measure: Gaussian measurement noise model
        }


        paramLogLikelihood= -sumError;
	//=========================================================
	return paramLogLikelihood;
}

void userComputeStatistics(Array1D<double> &parameterStatistics, Array1D<MCMC::chainstate> & parameterChain){
	//=========================================================
	//===== USER DEFINED ======================================
	//=========================================================

        //calculate statistics of the parameter chain

        //===========USER DEFINED: MEAN AND STD OF 2 PARAMETERS ======================================
        //mean of parameters 1&2
        for (int i=0;i<parameterChain.XSize(); i++){
                        parameterStatistics(0)+=parameterChain(i).state(0);
                        parameterStatistics(2)+=parameterChain(i).state(1);
        }
                parameterStatistics(0)=parameterStatistics(0)/(1.0*parameterChain.XSize());
                parameterStatistics(2)=parameterStatistics(2)/(1.0*parameterChain.XSize());
        //standard deviation of parameters 1&2
        for (int i=0;i<parameterChain.XSize(); i++){
                        parameterStatistics(1)+=pow(parameterStatistics(0) - parameterChain(i).state(0) ,2.0 );
                        parameterStatistics(3)+=pow(parameterStatistics(2) - parameterChain(i).state(1) ,2.0 );
        }
                parameterStatistics(1)=pow(parameterStatistics(1)/(1.0*parameterChain.XSize()) ,0.5);
                parameterStatistics(3)=pow(parameterStatistics(3)/(1.0*parameterChain.XSize()) ,0.5);
        //============================================================================


        return;
}
#endif  //USER_FUNCS_H_
