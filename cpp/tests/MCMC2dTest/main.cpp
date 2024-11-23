/* =====================================================================================

                      The UQ Toolkit (UQTk) version 3.1.5
                          Copyright (2024) NTESS
                        https://www.sandia.gov/UQToolkit/
                        https://github.com/sandialabs/UQTk

     Copyright 2024 National Technology & Engineering Solutions of Sandia, LLC (NTESS).
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
#include <iostream>
#include "math.h"
#include "Array1D.h"
#include "Array2D.h"
#include "mcmc.h"
#include "tmcmc.h"
#include "amcmc.h"
#include "ss.h"
#include "quad.h"
#include "dsfmt_add.h"
#include "arrayio.h"
#include "arraytools.h"
#include "assert.h"

using namespace std;

/*************************************************
Define LogPosterior function
*************************************************/
class LogPosterior: public LogPosteriorBase{
public:
	LogPosterior(){};
	~LogPosterior(){};
	double eval(Array1D<double>&);
};

// // Rosnebrock function
// double LogPosterior::eval(Array1D<double>& x){
//   double lnpost = -(1-x(0))*(1-x(0)) - 100*(x(1) - x(0)*x(0))*(x(1) - x(0)*x(0));
//   return lnpost;
// }

// Simple 2d Gaussian with zero mean and (.1,.8) variance
double LogPosterior::eval(Array1D<double>& x){
	double lnpost = -.5*(x(0)*x(0)/.01 + x(1)*x(1)/.64);
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
	int nCalls = 500000;
	Array1D<double> x(dim,0);

	LogPosterior L;
	// cout << "L.eval(x) = " << L.eval(x) << endl;

	/*************************************************
	Initiate and Run MCMC chain
	*************************************************/
	Array1D<double> g(dim,.1);

	AMCMC mchain(L);
	mchain.setChainDim(dim);
	//mchain.initMethod("am");
	mchain.initChainPropCovDiag(g);
	mchain.setSeed(13);
	// mchain.printChainSetup();
	// mchain.setOutputInfo("txt","chain.txt",nCalls,nCalls);
	mchain.setWriteFlag(0);
	mchain.runChain(nCalls,x);

    // Get chain states
	Array1D<MCMC::chainstate> chainstates;
	mchain.getFullChain(chainstates);

	// get mean from chainstates
	double mean_x1 = 0;
	double mean_x2 = 0;
	int nBurn = 3000;
	for (int i = nBurn; i < nCalls; i++){
		mean_x1 += chainstates(i).state(0);
		mean_x2 += chainstates(i).state(1);
	}
	mean_x1 *= 1./(nCalls-nBurn+1);
	mean_x2 *= 1./(nCalls-nBurn+1);
	cout << mean_x1 << endl;
	cout << mean_x2 << endl;

	// get variance
	double var_x1 = 0;
	double var_x2 = 0;
	for (int i = nBurn; i < nCalls; i=i+1){
		var_x1 += pow(chainstates(i).state(0) - mean_x1,2);
		var_x2 += pow(chainstates(i).state(1) - mean_x2,2);
	}
	var_x1 *= 1./(nCalls-nBurn);
	var_x2 *= 1./(nCalls-nBurn);
	cout << var_x1 << endl;
	cout << var_x2 << endl;

	// check variance
	assert(fabs((sqrt(var_x1) - .1)) < .01);
	assert(fabs((sqrt(var_x2) - .8)) < .01);

    /*************************************************
    Initiate and Run MCMC-SS chain
    *************************************************/
    Array1D<double> g1(dim,.1);

    SS mchain1(L);
    mchain1.setChainDim(dim);
    mchain1.initChainPropCovDiag(g1);
    mchain1.setSeed(13);
    // mchain.printChainSetup();
    // mchain.setOutputInfo("txt","chain.txt",nCalls,nCalls);
    mchain1.setWriteFlag(0);
    mchain1.runChain(nCalls,x);

    // Get chain states
    Array1D<MCMC::chainstate> chainstates1;
    mchain1.getFullChain(chainstates1);

    // get mean from chainstates
    double mean_x3 = 0;
    double mean_x4 = 0;
    for (int i = nBurn; i < nCalls; i++){
        mean_x3 += chainstates1(i).state(0);
        mean_x4 += chainstates1(i).state(1);
    }
    mean_x3 *= 1./(nCalls-nBurn+1);
    mean_x4 *= 1./(nCalls-nBurn+1);
    cout << mean_x3 << endl;
    cout << mean_x4 << endl;

    // get variance
    double var_x3 = 0;
    double var_x4 = 0;
    for (int i = nBurn; i < nCalls; i=i+1){
        var_x3 += pow(chainstates1(i).state(0) - mean_x3,2);
        var_x4 += pow(chainstates1(i).state(1) - mean_x4,2);
    }
    var_x3 *= 1./(nCalls-nBurn);
    var_x4 *= 1./(nCalls-nBurn);
    cout << var_x3 << endl;
    cout << var_x4 << endl;

    // check variance
    assert(fabs((sqrt(var_x3) - .1)) < .01);
    assert(fabs((sqrt(var_x4) - .8)) < .01);

    /*************************************************
     TMCMC 2d Test
     ***********************************************/

	/*************************************************
    Dimensionality and number of samples requested
    *************************************************/
    dim = 2;
    nCalls = 1000;

    /*************************************************
    Initiate and Run TMCMC
    *************************************************/
    TMCMC mchain2;
    mchain2.setChainDim(dim);
    mchain2.setSeed(1);
    mchain2.setWriteFlag(1);
    mchain2.initTMCMCCv(0.5);
    mchain2.initTMCMCNprocs(4);
    mchain2.runChain(nCalls);

    // Get chain states
    Array1D<MCMC::chainstate> chainstates2;
    mchain2.getFullChain(chainstates2);

    // get mean from chainstates
    double mean_x5 = 0;
    double mean_x6 = 0;
    for (int i = 0; i < nCalls; i++){
        mean_x5 += chainstates2(i).state(0);
        mean_x6 += chainstates2(i).state(1);
    }
    mean_x5 *= 1./nCalls;
    mean_x6 *= 1./nCalls;
    cout << mean_x5 << endl;
    cout << mean_x6 << endl;

    // check mean
    assert(fabs(mean_x5) < .05);
    assert(fabs(mean_x6) < .05);

    // get variance
    double var_x5 = 0;
    double var_x6 = 0;
    for (int i = 0; i < nCalls; i=i+1){
        var_x5 += pow(chainstates2(i).state(0) - mean_x5,2);
        var_x6 += pow(chainstates2(i).state(1) - mean_x6,2);
    }
    var_x5 *= 1./(nCalls-1);
    var_x6 *= 1./(nCalls-1);
    cout << var_x5 << endl;
    cout << var_x6 << endl;

    // check variance
    assert(fabs(var_x5 - 1./(1./1. + 1./0.01)) < .05);
    assert(fabs(var_x6 - 1./(1./1. + 1./0.64)) < .05);

    return 0;


}
