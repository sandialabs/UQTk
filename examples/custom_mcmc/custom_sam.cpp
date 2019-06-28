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
/*! \file custom_sam.cpp
*/

#include <unistd.h>
#include <sstream>
#include <map>
#include <math.h>
#include <iostream>
#include <string>

#include <getopt.h>


#include "func.h"
#include "post.h"
#include "mrv.h"
#include "inference.h"

#include "mcmc.h"
#include "tools.h"
#include "arrayio.h"
#include "arraytools.h"



using namespace std;

// Example log-Posterior function
double logPosterior(Array1D<double>& m, void* postinfo);



//  Main program: MCMC sampling
int main (int argc, char *argv[]) 
{
    // Any auxiliary information for posterior can be passed via void*
    double center=0.4;
    void* postinfo=(void*) &center;

    // Chain seed
    int seed=111;
    // Chain dimensionality
    int chdim=2;
    // Gamma factor for AMCMC
    double mcmcgamma=0.1;
    // Number of MCMC steps
    int nmcmc=10000;
    // Number of burnin steps
    int nburn=1000;
    // For good statistics, pick every nstep state
    int nstep=10;

    // Chain start 
    Array1D<double> chstart(chdim,5.0);
    // Standard deviation per dimension for initial non-adaptive part
    Array1D<double> chsig(chdim,0.5);

    // Instantiate chain object
    MCMC mchain(logPosterior,(void *) postinfo);
    cout << "Starting MCMC..................." << endl;
    // Set the seed
    mchain.setSeed(seed);
    // Set chain dimensionality
    mchain.setChainDim(chdim);
   
    // initialize other MCMC parameters
    mchain.initAMGamma(mcmcgamma);
    mchain.initAdaptSteps(nmcmc/10,10,nmcmc);
    mchain.initChainPropCovDiag(chsig);
        
    // print chain setup
    mchain.printChainSetup();
    mchain.runChain(nmcmc, chstart);

    // Map values to output
    Array1D<double> mapparam;
    // Clean (thinned and burn-in removed) chain output
    Array2D<double> pchain;

    // get MAP value
    mchain.getMode(mapparam);
    // thin the chain and remove burn-in
    mchain.getSamples(nburn, nstep, pchain);
    // append MAP values, too
    pchain.insertCol(mapparam,pchain.YSize());
    // transpose for convenient reading/plotting
    pchain=Trans(pchain);

    // Write clean chain
    write_datafile(pchain,"pchain.dat");
    // Write MAP parameter values
    write_datafile_1d(mapparam,"mapparam.dat");


  return 0;
}

// Example log Posterior
double logPosterior(Array1D<double>& m, void* postinfo)
{
    // Recast the auxiliary posterior information
    double* center=(double*) postinfo;

    // Gaussian log-posterior
    double logPost=-0.5*pow(m(0)-(*center),2.0)-0.5*pow(m(1)-(*center),2.0);

    return logPost;
}
