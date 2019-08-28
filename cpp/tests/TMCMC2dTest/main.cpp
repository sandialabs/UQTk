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
Begin main code
*************************************************/
int main(int argc, char ** argv){

/*************************************************
	Dimensionality and number of samples requested
	*************************************************/
	int dim = 2; 
	int nCalls = 1000;

	/*************************************************
	Initiate and Run MCMC chain
	*************************************************/
	Array1D<double> g(dim,.1);

	MCMC mchain;
	mchain.setChainDim(dim);
	mchain.initMethod("tmcmc");
	mchain.initChainPropCovDiag(g);
	mchain.setSeed(13);
	mchain.setWriteFlag(1); 
	mchain.runChain(nCalls);

	// Get chain states
	Array1D<MCMC::chainstate> chainstates; 
	mchain.getFullChain(chainstates);

	// get mean from chainstates
	double mean_x1 = 0;
	double mean_x2 = 0;
	for (int i = 0; i < nCalls; i++){
		mean_x1 += chainstates(i).state(0);
		mean_x2 += chainstates(i).state(1);
	}
	mean_x1 *= 1./nCalls;
	mean_x2 *= 1./nCalls;
	cout << mean_x1 << endl;
	cout << mean_x2 << endl;

	// get variance
	double var_x1 = 0;
	double var_x2 = 0;
	for (int i = 0; i < nCalls; i=i+1){
		var_x1 += pow(chainstates(i).state(0) - mean_x1,2);
		var_x2 += pow(chainstates(i).state(1) - mean_x2,2);
	}
	var_x1 *= 1./(nCalls-1);
	var_x2 *= 1./(nCalls-1);
	cout << var_x1 << endl;
	cout << var_x2 << endl;

	// check variance
	assert(fabs(var_x1 - 1./(1./1. + 1./0.01)) < .05);
	assert(fabs(var_x2 - 1./(1./1. + 1./0.64)) < .05);

	return 0; 

}
