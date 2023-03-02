/* =====================================================================================

                      The UQ Toolkit (UQTk) version 3.1.3
                          Copyright (2023) NTESS
                        https://www.sandia.gov/UQToolkit/
                        https://github.com/sandialabs/UQTk

     Copyright 2023 National Technology & Engineering Solutions of Sandia, LLC (NTESS).
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
#include <iostream>
#include "math.h"
#include "Array1D.h"
#include "Array2D.h"
#include "mcmc.h"
#include "amcmc.h"
#include "quad.h"
#include "dsfmt_add.h"
#include "arrayio.h"
#include "arraytools.h"
#include "assert.h"

using namespace std;

/*************************************************
Define LogPosterior functions
*************************************************/
class LogPosterior: public LogPosteriorBase{
public:
	LogPosterior(){};
	~LogPosterior(){};
	double eval(Array1D<double>&);
};
class LogPosterior_in: public LogPosteriorBase{
public:
	LogPosterior_in(){};
	~LogPosterior_in(){};
	double eval(Array1D<double>&);
};
class LogPosterior2: public LogPosteriorBase{
public:
	LogPosterior2(){};
	~LogPosterior2(){};
	double eval(Array1D<double>&);
};

// Rosnebrock function
double LogPosterior::eval(Array1D<double>& x){
	LogPosterior_in L;
	int dim = 2;
	int nCalls = 100;
	Array1D<double> g(dim,.1);
	AMCMC mchaintemp(L);
	mchaintemp.setChainDim(dim);
	mchaintemp.initChainPropCovDiag(g);
	mchaintemp.setWriteFlag(0);
	mchaintemp.setSeed(10);
	mchaintemp.runChain(nCalls,x);
    Array2D<double> samples;
	mchaintemp.getSamples(samples);
	samples = Trans(samples);
	// printarray(samples);

	double lnpost = -(1-x(0))*(1-x(0)) - 100*(x(1) - x(0)*x(0))*(x(1) - x(0)*x(0));
  return lnpost;
}
// Rosnebrock function
double LogPosterior_in::eval(Array1D<double>& x){
	double lnpost = -.5*(x(0)*x(0)/.01 + x(1)*x(1)/.64);
  return lnpost;
}
// Rosnebrock function
double LogPosterior2::eval(Array1D<double>& x){
	double lnpost = -(1-x(0))*(1-x(0)) - 100*(x(1) - x(0)*x(0))*(x(1) - x(0)*x(0));
  return lnpost;
}

/*************************************************
Begin main code
*************************************************/
int main(int argc, char ** argv){

	/*************************************************
	Initial start for MCMC chain
	and set LogPosterior function
	*************************************************/
	int dim = 2;
	int nCalls = 100;
	Array1D<double> x(dim,0);

	/*************************************************
	Initiate and Run MCMC chain
	*************************************************/
	LogPosterior L;
	Array1D<double> g(dim,.1);
	AMCMC mchain(L);
	mchain.setChainDim(dim);
	mchain.initChainPropCovDiag(g);
	mchain.setOutputInfo("txt","chain.txt",nCalls,nCalls);
	mchain.setWriteFlag(0);
	mchain.runChain(nCalls,x);
    Array2D<double> samples;
	mchain.getSamples(samples);
	samples = Trans(samples);
	// printarray(samples);

	/*************************************************
	Initiate and Run MCMC chain
	*************************************************/
	LogPosterior2 L2;
	Array1D<double> g2(dim,.1);
	AMCMC mchain2(L2);
	mchain2.setChainDim(dim);
	mchain2.initChainPropCovDiag(g2);
	mchain2.setOutputInfo("txt","chain.txt",nCalls,nCalls);
	mchain2.setWriteFlag(0);
	mchain2.runChain(nCalls,x);
    Array2D<double> samples2;
	mchain2.getSamples(samples2);
	samples2 = Trans(samples2);
	// printarray(samples2);

	// test to make sure two MCMC objects do not
	// affect each others random number generation
	for (int i = 0; i < nCalls; i++){
		for (int n = 0; n < dim; n++){
			assert(samples(i,n) - samples2(i,n) == 0) ;
		}
	}



	return 0;

}
