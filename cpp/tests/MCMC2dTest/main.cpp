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
#include <iostream>
#include "math.h"
#include "Array1D.h"
#include "Array2D.h"
#include "mcmc.h"
#include "quad.h"
#include "dsfmt_add.h"
#include "arrayio.h"
#include "arraytools.h"
#include "assert.h"

using namespace std;

/*************************************************
Define Likelihood function
*************************************************/
class Likelihood: public LikelihoodBase{
public:
	Likelihood(){};
	~Likelihood(){};
	double eval(Array1D<double>&);
};

// // Rosnebrock function
// double Likelihood::eval(Array1D<double>& x){
//   double lnpost = -(1-x(0))*(1-x(0)) - 100*(x(1) - x(0)*x(0))*(x(1) - x(0)*x(0));
//   return lnpost;
// }

// Simple 2d Gaussian with zero mean and (.1,.8) variance
double Likelihood::eval(Array1D<double>& x){
	double lnpost = -.5*(x(0)*x(0)/.01 + x(1)*x(1)/.64);
	return lnpost;
}

/*************************************************
Begin main code
*************************************************/
int main(int argc, char ** argv){

/*************************************************
	Initial start for MCMC chain
	and set Likelihood function
	*************************************************/
	int dim = 2;
	int nCalls = 500000;
	Array1D<double> x(dim,0);

	Likelihood L;
	// cout << "L.eval(x) = " << L.eval(x) << endl;

	/*************************************************
	Initiate and Run MCMC chain
	*************************************************/
	Array1D<double> g(dim,.1);

	MCMC mchain(L);
	mchain.setChainDim(dim);
	mchain.initMethod("am");
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

	return 0;

}
