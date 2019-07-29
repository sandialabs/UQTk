/* =====================================================================================
                     The UQ Toolkit (UQTk) version @UQTKVERSION@
                     Copyright (@UQTKYEAR@) Sandia Corporation
                     http://www.sandia.gov/UQToolkit/

     Copyright (@UQTKYEAR@) Sandia Corporation. Under the terms of Contract DE-AC04-94AL85000
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
Define Likelihood functions
*************************************************/
class Likelihood: public LikelihoodBase{
public:
	Likelihood(){};
	~Likelihood(){}; 
	double eval(Array1D<double>&);
};
class Likelihood_in: public LikelihoodBase{
public:
	Likelihood_in(){};
	~Likelihood_in(){}; 
	double eval(Array1D<double>&);
};
class Likelihood2: public LikelihoodBase{
public:
	Likelihood2(){};
	~Likelihood2(){}; 
	double eval(Array1D<double>&);
};

// Rosnebrock function
double Likelihood::eval(Array1D<double>& x){
	Likelihood_in L; 
	int dim = 2; 
	int nCalls = 100; 
	Array1D<double> g(dim,.1);
	MCMC mchaintemp(L); 
	mchaintemp.setChainDim(dim);
	mchaintemp.initMethod("am");
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
double Likelihood_in::eval(Array1D<double>& x){
	double lnpost = -.5*(x(0)*x(0)/.01 + x(1)*x(1)/.64);
  return lnpost; 
}
// Rosnebrock function
double Likelihood2::eval(Array1D<double>& x){
	double lnpost = -(1-x(0))*(1-x(0)) - 100*(x(1) - x(0)*x(0))*(x(1) - x(0)*x(0));
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
	int nCalls = 100;
	Array1D<double> x(dim,0);

	/*************************************************
	Initiate and Run MCMC chain
	*************************************************/
	Likelihood L; 
	Array1D<double> g(dim,.1);
	MCMC mchain(L); 
	mchain.setChainDim(dim);
	mchain.initMethod("am");
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
	Likelihood2 L2; 
	Array1D<double> g2(dim,.1);
	MCMC mchain2(L2); 
	mchain2.setChainDim(dim);
	mchain2.initMethod("am");
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
