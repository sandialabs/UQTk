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
/// \file inference.cpp
/// \author K. Sargsyan  2016 -
/// \brief Model inference tools

#include <unistd.h>
#include <sstream>
#include <map>
#include <iostream>
#include <string>
#include <math.h>

#include "func.h"
#include "post.h"
#include "mrv.h"
#include "inference.h"

#include "mcmc.h"
#include "tools.h"
#include "arrayio.h"
#include "arraytools.h"

using namespace std;

void infer_model(Array1D< Array2D<double> (*)(Array2D<double>&, Array2D<double>&, Array2D<double>&, void *) > forwardFuncs, void* funcInfo,
	string likType,
	string priorType, double priora, double priorb,
	Array2D<double>& xdata,Array2D<double>& ydata, Array2D<double>& xgrid,
	int dataNoiseInference, Array1D<double>& datanoise_array,
	int pdim, int order,Array1D<int>& rndInd,Array2D<double>& fixindnom, string pdfType, string pcType,
	int seed, int nmcmc, double mcmcgamma, bool optimflag, Array1D<double>& chstart, Array1D<double>& chsig,
	double likParam, int likParam_int,
	Array2D<double>& pgrid,Array2D<double>& pchain, int nburn, int nstep,
	Array1D<double>& mapparam, Array1D<double>& datavar_map,
	Array1D<double>& pmean_map, Array1D<double>& pvar_map,
	Array1D<double>& fmean_map, Array1D<double>& fvar_map,
	Array1D<double>& postave_datavar,
	Array1D<double>& p_postave_mean, Array1D<double>& p_postave_var, Array1D<double>& p_postvar_mean,
   	Array2D<double>& f_postsam_mean, Array1D<double>& f_postave_mean, Array1D<double>& f_postave_var, Array1D<double>& f_postvar_mean,
   	Array2D<double>& paramPCcfs){

	// Some input argument checks and fixes
	if (likType=="classical" or likType=="koh" or likType=="eov"){
		if (pdfType!="full"){
			cout << "Classical likelihood or koh or eov likelihood types require full PDF type." << endl;
			pdfType="full";
		}
		if (order!=0){
			cout << "Classical or koh or eov likelihood types require zeroth order of PDF PC." << endl;
			order=0;
		}
	}

	// Read x-grid size
    int nxgr = xgrid.XSize();
    assert(xgrid.YSize()==xdata.YSize());

    // Match the likelihood type with the appropriate derived class
    map<string, Post*> lik_dict;
    lik_dict["full"]=new Lik_Full(likParam,likParam_int);
    lik_dict["marg"]=new Lik_Marg(likParam,likParam_int);
	lik_dict["mvn"]=new Lik_MVN(likParam);
	lik_dict["gausmarg"]=new Lik_GausMarg();
	lik_dict["abc"]=new Lik_ABC(likParam);
	lik_dict["abcm"]=new Lik_ABCm(likParam);
	lik_dict["koh"]=new Lik_Koh(likParam);
	lik_dict["classical"]=new Lik_Classical();
    lik_dict["eov"]=new Lik_Eov();

	if (lik_dict.count(likType)==0){
		cout << "Likelihood type " << likType << " is not found. Exiting." << endl;
		exit(1);
	}
	Post* mypost=lik_dict[likType];

	// Parse the data noise information
	if (dataNoiseInference==0)
		mypost->setDataNoise(datanoise_array);
	else if (dataNoiseInference==1)
		mypost->inferDataNoise();
	else if (dataNoiseInference==2)
		mypost->inferLogDataNoise();

	// Set data
	mypost->setData(xdata,ydata);
	// Set model
	mypost->setModel(forwardFuncs, fixindnom, funcInfo);
	// Set the randomized model parameter input
	mypost->setModelRVinput(pdim, order,rndInd, pdfType, pcType);
	// Set prior
	mypost->setPrior(priorType,priora,priorb);

	// Get the chain dimensionality
	int chdim=mypost->getChainDim();
	cout << "MCMC dimensionality : " << chdim << endl;
	if (chstart.Length()==0)
		chstart.Resize(chdim,0.1);
    if (chsig.Length()==0)
		chsig.Resize(chdim,0.01);

    assert (chdim==chstart.Length() and chdim==chsig.Length());

    // Overwrite to chain initial state
    if (dataNoiseInference>0)
		chstart(chdim-1)=datanoise_array(0);

    // Initialize chain
    MCMC mchain(LogPosterior,(void *) mypost);
    mchain.setSeed(seed);
    mchain.setChainDim(chdim);

    // If given parameter grid, compute exact posterior, up to a normalization
    int npgr=pgrid.XSize();
    if (npgr>0){

	    cout << "Compute posterior at a grid..........." << endl;
	    assert(pgrid.YSize()==chdim);
	    double sum=0.0;
	    Array1D<double> pdens(npgr,0.e0);
	    Array1D<double> pdens_log(npgr,0.e0);
	    for(int i=0;i<npgr;i++){
	    	if ((i+1)%100==0)
		        cout << "Grid point " << i+1 << " / " << npgr << endl;
	        Array1D<double> pgrid1p;
	        getRow(pgrid,i,pgrid1p);
	        pdens(i)=exp(mchain.evalLogPosterior(pgrid1p));
	        pdens_log(i)=mchain.evalLogPosterior(pgrid1p);
	        double normalize=1.0;
	        sum+=pdens(i)*normalize;
	    }
	   	write_datafile_1d(pdens,"pdens_unnorm.dat");

	    for(int i=0;i<npgr;i++)
	        pdens(i)/=sum;
	    write_datafile_1d(pdens,"pdens.dat");
	    write_datafile_1d(pdens_log,"pdens_log.dat");
	}

	// if chain array is empty, do MCMC and populate it
	// otherwise, postprocess only
	int npost=pchain.XSize();
	if (npost==0){


	    // initialize other MCMC parameters
	    mchain.initAMGamma(mcmcgamma);
	    if (nmcmc<0){
	    	nmcmc=-nmcmc;
	    	mchain.initAdaptSteps(nmcmc+1,0,nmcmc+1);
	    }
	    else
		    mchain.initAdaptSteps(nmcmc/10,10,nmcmc);

	    mchain.initChainPropCovDiag(chsig);

	    // print chain setup
	    mchain.printChainSetup();
	    // run initial optimization
	    if (optimflag)
		    mchain.runOptim(chstart);
	    // run MCMC
	    if (nmcmc>0){
	    	cout << "Starting MCMC..................." << endl;
	        mchain.runChain(nmcmc, chstart);
	    }

	    if (optimflag || nmcmc>0){
		    // get MAP value
		    mchain.getMode(mapparam);
		    // thin the chain and remove burn-in
	        mchain.getSamples(nburn, nstep, pchain);
	        // append MAP values, too
	        pchain.insertCol(mapparam,pchain.YSize());
	        // transpose for convenient reading/plotting
	        pchain=Trans(pchain);
	    }
    }
    else{
    	assert(pchain.YSize()==chdim);
    	nmcmc=0;
    	optimflag=false;
    	cout << "Posterior samples provided. Forcing nmcmc=0 and optimflag=false." << endl; // move this up?
    	// if pchain is provided, its last row is read as MAP
		getRow(pchain, pchain.XSize()-1, mapparam);
    }

    // Switch the model, if needed
    int pred_mode=1;
    funcInfo=(void*) &pred_mode;
    mypost->setModel(forwardFuncs, fixindnom, funcInfo);


    // Compute function prediction moments at xgrid and parameter moments at MAP values
    // (skip the covariances)
    int npc=0;
    Array2D<double> mapparamPCcf;
	Array1D<double> mapparamPCcf_flat;
    //if (optimflag || nmcmc>0)
    //{
	    Array2D<double> fcov_map,pcov_map;//dummy
	    mypost->momForwardFcn(mapparam, xgrid, fmean_map, fvar_map, false, fcov_map);
		mypost->momParam(mapparam, pmean_map, pvar_map, false, pcov_map);

		// Compute parameter PC coefficients at MAP values
	    mapparamPCcf=mypost->getParamPCcf(mapparam);
	    flatten(mapparamPCcf,mapparamPCcf_flat);
	    npc=mapparamPCcf_flat.Length();
	    // Compute MAP data variance
	    Array1D<double> tmp=mypost->dataSigma(mapparam(mapparam.Length()-1));
		datavar_map=dotmult(tmp,tmp);
	//}

	cout << "===========================================================" << endl;

	npost=pchain.XSize();
    // Repeat postprocessing for each posterior sample
    cout << endl << "Going through " << npost << " posterior samples" << endl;
   	paramPCcfs.Resize(npc,0);
   	postave_datavar.Resize(xdata.XSize(),0.0);
	for(int ipost=0;ipost<npost;ipost++){
		if (ipost % max(1,int(npost/10)) == 0)
			cout << ipost << "/" << npost << endl;

	    Array1D<double> fmean_mcmc,fvar_mcmc, pmean_mcmc,pvar_mcmc;
	    Array2D<double> fcov_mcmc,pcov_mcmc;
	    // Get current parameter column
		Array1D<double> mcmcparam;
        getRow(pchain,ipost,mcmcparam);

 		// Compute function prediction moments at xgrid and parameter moments at current chain values
    	// (skip the covariances)
        mypost->momForwardFcn(mcmcparam, xgrid, fmean_mcmc, fvar_mcmc, false, fcov_mcmc);
 		mypost->momParam(mcmcparam, pmean_mcmc, pvar_mcmc, false, pcov_mcmc);

		// Compute parameter PC coefficients at current chain values
        Array2D<double> mcmcparamPCcf=mypost->getParamPCcf(mcmcparam);

    	Array1D<double> mcmcparamPCcf_flat;
    	flatten(mcmcparamPCcf,mcmcparamPCcf_flat);
    	paramPCcfs.insertCol(mcmcparamPCcf_flat,ipost);
    	// Compute posterior average of data variance
    	// note: if there is no data variance inference, then this is extra work
    	Array1D<double> tmp=mypost->dataSigma(mcmcparam(mcmcparam.Length()-1));
		Array1D<double> tmp2=dotmult(tmp,tmp);
		scaleinplace(tmp2,1./npost);
		addinplace(postave_datavar,tmp2);

		// Compute posterior average of parameter mean and variance, and posterior variance of parameter mean
 		for (int ip=0;ip<pdim;ip++){
 			p_postave_mean(ip)+=pmean_mcmc(ip)/npost;
            p_postave_var(ip)+=pvar_mcmc(ip)/npost;
            p_postvar_mean(ip)+=pow(pmean_mcmc(ip),2.)/npost;
 		}
		// Compute posterior average of function prediction mean and variance, and posterior variance of function prediction mean
        for (int ig=0;ig<nxgr;ig++){
            f_postave_mean(ig)+=fmean_mcmc(ig)/npost;
            f_postave_var(ig)+=fvar_mcmc(ig)/npost;
            f_postvar_mean(ig)+=pow(fmean_mcmc(ig),2.)/npost;
	    }
	    f_postsam_mean.insertCol(fmean_mcmc,ipost);
	}

	// Subtract mean-squared to get to variance
 	for (int ip=0;ip<pdim;ip++)
 		p_postvar_mean(ip)-=pow(p_postave_mean(ip),2.);

	// Subtract mean-squared to get to variance
	for (int ig=0;ig<nxgr;ig++)
	    f_postvar_mean(ig)-=pow(f_postave_mean(ig),2.);

	// Append the MAP column of parameter PC coefficients
	paramPCcfs.insertCol(mapparamPCcf_flat,npost);

	// Some output information
	if (dataNoiseInference>0){
		cout << "Max a posteriori of data variance " << datavar_map(0) << endl;
		cout << "Posterior average of data variance " << postave_datavar(0) << endl;
	}

	return;
}

// Log-posterior function given chain state and auxiliary parameters
double LogPosterior(Array1D<double>& m, void* mypost_void)
{
    Post* mypost = (Post*) mypost_void;

	double logprior=mypost->evalLogPrior(m);
	double loglik=0.0;
	if (logprior>-1.e80){
		loglik=mypost->evalLogLik(m);

	}


	double logpost=loglik+logprior;

	return logpost;
}

