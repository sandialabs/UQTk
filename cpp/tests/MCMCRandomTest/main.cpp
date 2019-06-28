/* =====================================================================================
                     The UQ Toolkit (UQTk) version 3.0.4
                     Copyright (2017) Sandia Corporation
                     http://www.sandia.gov/UQToolkit/

     Copyright (2017) Sandia Corporation. Under the terms of Contract DE-AC04-94AL85000
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
Define Likelihood function
*************************************************/
class Likelihood: public LikelihoodBase{
public:
	Likelihood(){};
	~Likelihood(){}; 
	double eval(Array1D<double>&);
};

// Rosnebrock function
double Likelihood::eval(Array1D<double>& x){
	 double lnpost = -(1-x(0))*(1-x(0)) - 100*(x(1) - x(0)*x(0))*(x(1) - x(0)*x(0));
	// double lnpost = -.5*(x(0)*x(0)/.01 + x(1)*x(1)/.64);
	// // 	'''
	// // sample from exp(-.5*(x**2/.1**2 - y**2/.8**2))
	// // '''
	// // y1 = x[0]
	// // y2 = x[1]
	// // return -.5*(y1**2/.1**2 + y2**2/.8**2)
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

	Likelihood L; 
	cout << "L.eval(x) = " << L.eval(x) << endl;

	/*************************************************
	Initiate and Run MCMC chain
	*************************************************/
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
	Initiate and Run a second MCMC chain
	*************************************************/
	MCMC mchain2(L); 
	mchain2.setChainDim(dim);
	mchain2.initMethod("am");
	mchain2.initChainPropCovDiag(g);
	// mchain2.setSeed(130);
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
