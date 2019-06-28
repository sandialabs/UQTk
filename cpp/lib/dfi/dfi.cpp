/* =====================================================================================
                     The UQ Toolkit (UQTk) version 3.0.4
                     Copyright (2017) Sandia Corporation
                     http://www.sandia.gov/UQToolkit/

     Copyright (2013) Sandia Corporation. Under the terms of Contract DE-AC04-94AL85000
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
#include <math.h>
#include "error_handlers.h"
#include "Array1D.h"
#include "Array2D.h"
#include "quad.h"
#include "tools.h"
#include "arraytools.h"
#include "mcmc.h"
#include "gen_defs.h"
#include "lbfgs_routines.h"
#include "deplapack.h"
#include "depblas.h"

#include "dfi.h"

/*********************************************
Define the inner dfi mcmc chain
*********************************************/
DFIInner::DFIInner(DFISetupBase& d){
    // Constructor sets up the logposterior 
    // from the DFISetupBase class
    d_ = &d; 

    this->nCalls = 30000; 
    this->nBurn = 10000; 

}
double DFIInner::eval(Array1D<double>& beta){
    // evaluates logposterior as determined by DFISetupBase
    return d_->f(beta,z_,sigma_);
}
void DFIInner::getSamples(){ 
    // ndim must be set before running this routine
    // must also set initial start of chain beta0
    MCMC inner_chain(*this);
    inner_chain.setWriteFlag(0); 
    inner_chain.setChainDim(ndim);
    inner_chain.initMethod("am");
    inner_chain.initAdaptSteps(1000,1000,10000);
    // inner_chain.initChainPropCovDiag(gammas_); //optional

    samples_.Clear();
    samples_.Resize(nCalls);
    inner_chain.setOutputInfo("bin","chain.bin",nCalls,nCalls);
    inner_chain.resetChainState(); 
    inner_chain.runChain(nCalls,beta0_);
    inner_chain.getFullChain(samples_);

}
double DFIInner::S(){
    // call user defined function of the inner samples
    this->getSamples();
    return d_->S(samples_);
}

/*************************************************
Main DFI Outer class
*************************************************/
DFI::DFI(DFIInner& dfi_inner){
    dfi_inner_ = &dfi_inner;

    nBeta = dfi_inner_->ndim;
    // sdim = dfi_inner_->sigma_.Length();
    // zdim = dfi_inner_->z_.Length();
}
double DFI::eval(Array1D<double>& zs){
    // get dimensions
    int ndim = zdim + sdim;

    // unravel z and s
    z_.Resize(zdim,0.0); 
    for (int i = 0; i < zdim; i++){
        z_(i) = zs(i); 
    } 
    int scount = 0;
    sigma_.Resize(sdim,0.0);
    for (int i = zdim; i < ndim; i++){
        sigma_(scount) = zs(i);
        scount += 1;  
    }
    dfi_inner_->z_ = z_;
    dfi_inner_->sigma_ = sigma_;  
    
    // get inner samples then evaluate S(z)
    dfi_inner_->getSamples(); 

    // flag to return -inf if sigma is <= 0
    int sigma_flag = 0;
    for (int k = 0; k < sdim; k++){
        if (sigma_(k) <= 0){
            sigma_flag = 1; 
        }
    }

    if (sigma_flag == 1){
        return -1e12; 
    }
    else{
        return dfi_inner_->S();
    }
}
void DFI::runChain(int nCalls, Array1D<double> gammas, Array1D<double> start, int seed, int node){
    
    MCMC outer_chain(*this);
    outer_chain.setWriteFlag(1);
    outer_chain.setSeed(seed);
    outer_chain.setChainDim(start.Length());
    outer_chain.initMethod("ss");
    outer_chain.initChainPropCovDiag(gammas);

    char fn[100];
    snprintf(fn,sizeof fn,"outerchain%03i.txt",node);
    outer_chain.setOutputInfo("txt",fn,5,1);

    outer_chain.runChain(nCalls,start);
}

void DFI::getMLE(Array1D<double>& xstart){

    MCMC outer_chain(*this);
    outer_chain.setChainDim(xstart.Length());
    outer_chain.runOptim(xstart);
    // xstart(0) = 10000; 

  return;
}
//********************************************************************
//********************************************************************

