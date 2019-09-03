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

using namespace std; 

/*************************************************
Begin main code
*************************************************/
int main(int argc, char ** argv){

/*************************************************
	Dimensionality and number of samples requested
	*************************************************/
	int dim = 3; 
	int nCalls = 5000;

	/*************************************************
	Initiate and Run TMCMC
	*************************************************/
	Array1D<double> g(dim,.1);

	MCMC mchain;
	mchain.setChainDim(dim);
	mchain.initMethod("tmcmc");
	mchain.setSeed(1);
	mchain.setWriteFlag(1); 
    mchain.setOutputInfo("txt","tmcmc_chain.dat",1,1);
    mchain.initTMCMCNprocs(4);
    mchain.initTMCMCCv(0.1);
	mchain.runChain(nCalls);

	return 0;
}
