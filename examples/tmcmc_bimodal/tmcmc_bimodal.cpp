/* =====================================================================================

                      The UQ Toolkit (UQTk) version 3.1.0
                          Copyright (2020) NTESS
                        https://www.sandia.gov/UQToolkit/
                        https://github.com/sandialabs/UQTk

     Copyright 2020 National Technology & Engineering Solutions of Sandia, LLC (NTESS).
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
