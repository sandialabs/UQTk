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

     Questions? Contact the UQTk Developers at https://github.com/sandialabs/UQTk/discussions
     Sandia National Laboratories, Livermore, CA, USA
===================================================================================== */
/// \file post.cpp
/// \author K. Sargsyan  2016 -
/// \brief Posterior computation class

#include <math.h>
#include <cfloat>
#include <assert.h>

#include "func.h"
#include "mrv.h"
#include "post.h"

#include "tools.h"

#include "PCSet.h"
#include "error_handlers.h"
#include "arrayio.h"
#include "arraytools.h"

#define SIG_MAX 10


// Constructor
Post::Post()
{
    this->extraInferredParams_=0;
}


// Set the x- and y-data
void Post::setData(Array2D<double>& xdata,Array2D<double>& ydata)
{
    int ndata = ydata.XSize();
    int ns = ydata.YSize();

    Array1D<Array1D<double>> ydata2(ndata);
    for(int i=0;i<ndata;i++){
        getRow(ydata, i, ydata2(i));
    }

    this->setData(xdata, ydata2);

	return;
}

// Set the x- and y-data in cases when y-data has different size for different x
void Post::setData(Array2D<double>& xdata,Array1D<Array1D<double> >& ydata)
{
    this->xData_=xdata;
    this->yData_=ydata;

    this->nData_=this->xData_.XSize();
    this->xDim_=this->xData_.YSize();
    assert(this->nData_==this->yData_.XSize());

    this->nEachs_.Resize(this->nData_,0);
    this->yDatam_.Resize(this->nData_,0.e0);
    for(int i=0;i<this->nData_;i++){
        this->nEachs_(i) = this->yData_(i).XSize();
        for(int j=0;j<this->nEachs_(i);j++){
            this->yDatam_(i)+=this->yData_(i)(j);
        }
        this->yDatam_(i)/=this->nEachs_(i);
    }


    return;
}

// Set the magnitude of data noise
void Post::setDataNoise(Array1D<double>& sigma)
{
	this->inferDataNoise_=false;
	this->dataNoiseLogFlag_=false;

	this->dataNoiseSig_=sigma;

    this->extraInferredParams_+=0;

	return;
}

// Indicate inference of data noise stdev
void Post::inferDataNoise()
{

	this->inferDataNoise_=true;
	this->dataNoiseLogFlag_=false;

   this->extraInferredParams_+=1;

	return;
}

// Indicate inference of log of data noise stdev
void Post::inferLogDataNoise()
{
	this->inferDataNoise_=true;
	this->dataNoiseLogFlag_=true;

	this->extraInferredParams_+=1;
	return;
}

// Get data noise, whether inferred or fixed
Array1D<double> Post::dataSigma(double m_last)
{
	Array1D<double> sig;
	if (inferDataNoise_){
		if (dataNoiseLogFlag_){
			sig.Resize(this->nData_,exp(m_last));
		}
		else {
			sig.Resize(this->nData_,m_last);
		}
	}
	else
		sig=this->dataNoiseSig_;



	return sig;
}

// Set a pointer to the forward model f(p,x)
void Post::setModel(Array1D< Array2D<double> (*)(Array2D<double>&, Array2D<double>&, Array2D<double>&, void *) > forwardFuncs, Array2D<double>& fixindnom, void* funcInfo)
{
	this->forwardFcns_=forwardFuncs;
    this->ncat_=forwardFuncs.Length();
    this->fixIndNom_=fixindnom;
    this->funcinfo_=funcInfo;
	return;
}

// Set model input parameters' randomization scheme
void Post::setModelRVinput(int pdim,int order,Array1D<int>& rndInd,string pdfType,string rvpcType)
{

	this->pDim_=pdim;
	this->rndInd_=rndInd;
	this->pdfType_=pdfType;
    this->rvpcType_=rvpcType;

	this->Mrv_ = new Mrv(this->pDim_,this->pdfType_, this->rndInd_,order,this->rvpcType_);
	this->chDim_=this->Mrv_->Parametrize()+this->extraInferredParams_+this->ncat_-1;
	this->Mrv_->getBounds(this->lower_,this->upper_);

    for (int i=0;i<this->ncat_-1;i++){
        this->lower_.PushBack(0.0);
        this->upper_.PushBack(1.0);
    }

	if (inferDataNoise_ and dataNoiseLogFlag_){
		this->lower_.PushBack(-SIG_MAX);
		this->upper_.PushBack(SIG_MAX);
	}
	if (inferDataNoise_ and ~dataNoiseLogFlag_){
		this->lower_.PushBack(0.e0);
		this->upper_.PushBack(SIG_MAX);
	}


  	return;

}

// Get the dimensionailty of the posterior function
int Post::getChainDim()
{
	return chDim_;
}

// Set the prior type and its parameters
void Post::setPrior(string priorType, double priora, double priorb)
{
	this->priorType_=priorType;
	this->priora_=priora;
	this->priorb_=priorb;

	return;
}

// Evaluate log-prior
double Post::evalLogPrior(Array1D<double>& m)
{
    double pi=4.*atan(1.);

	double logPrior=0.0;
	if (this->priorType_=="uniform"){

        Array1D<int> pctermid;
        this->Mrv_->getPCTermId(pctermid);

        for (int ic=0;ic<pctermid.Length();ic++){
            if (pctermid(ic)==0)
                if (m(ic)<this->priora_ or m(ic)>this->priorb_)
                    return -1.e80;
        }


		logPrior=0.0;
	}

    else if (this->priorType_=="uniform_LUpci"){

        Array1D<int> pctermid;
        this->Mrv_->getPCTermId(pctermid);

        for (int ic=0;ic<pctermid.Length();ic++){
            if (pctermid(ic)==0){
                if (ic<pctermid.Length()-1 and pctermid(ic+1)!=0){
                    if (m(ic)-fabs(m(ic+1))<this->priora_ or m(ic)+fabs(m(ic+1))>this->priorb_)
                        return -1.e80;
                }
                else{
                    if (m(ic)<this->priora_ or m(ic)>this->priorb_)
                        return -1.e80;
                }
            }
        }



        logPrior=0.0;
    }

	else if (this->priorType_=="normal"){

   	    Array1D<int> pctermid;
        this->Mrv_->getPCTermId(pctermid);

        for (int ic=0;ic<pctermid.Length();ic++){
	        double ps=this->priorb_;
	        double mism=m(ic)-this->priora_;
	        logPrior-=( 0.5*mism*mism/(ps*ps) + 0.5*log(2.*pi) + log(ps) );
		}
    }

    else if (this->priorType_=="inverse"){

        Array1D<int> pctermid;
        this->Mrv_->getPCTermId(pctermid);

        for (int ic=0;ic<pctermid.Length();ic++){
            if (pctermid(ic)==0)
                if (m(ic)<this->priora_ or m(ic)>this->priorb_)
                    return -1.e80;
            if (pctermid(ic)!=0)
                logPrior-=log(fabs(m(ic)));
        }
    }

    else if (this->priorType_=="wishart"){

        double nu = this->pDim_ + this->priora_;
        double theta = this->priorb_;



        bool covFlag = true;
        Array1D<double> parMean,parVar;
        Array2D<double> parCov;
        this->momParam(m, parMean, parVar, covFlag, parCov);

        //Array2D<double> paramPCcf=getParamPCcf(m);//npc x pdim
        //this->Mrv_->computeMoments(paramPCcf, fcnMean,fcnVar,covFlag, fcnCov);    //funccf is npc x nx


        logPrior = 0.5 * this->pDim_ * nu + 0.5 * (nu - pDim_ - 1.0) * logdeterm(parCov) - theta * trace(parCov);

        //double a = 0.5*(nu-this->pDim_+1.0);
        //double b = theta;
        //double logprior2=a*log(b) + (a-1.)*log(m(1)) - b*m(1) - lgamma(a)
    }

    else if (this->priorType_=="jeffreys"){




        bool covFlag = true;
        Array1D<double> parMean,parVar;
        Array2D<double> parCov;
        this->momParam(m, parMean, parVar, covFlag, parCov);


        logPrior = - 0.5 * (this->pDim_ +1.) * logdeterm(parCov);

    }
    else {
        cout << "Prior type " << this->priorType_ << " unrecognized. Exiting." << endl;
        exit(1);
    }

    return logPrior;
}

// Extract parameter PC coefficients from a posterior input
Array2D<double> Post::getParamPCcf(Array1D<double>& m)
{
    Array1D<double> modelRVparams;
    for (int i=0;i<m.Length()-this->extraInferredParams_-this->ncat_+1;i++)
        modelRVparams.PushBack(m(i));

    Array2D<double> paramPCcf=this->Mrv_->getMultiPCcf(modelRVparams);

    return paramPCcf;
}

// Sample model parameters given posterior input
Array2D<double> Post::samParam(Array1D<double>& m, int ns)
{
    Array2D<double> paramPCcf=this->getParamPCcf(m);

	Array2D<double> paramSamples=this->Mrv_->mcParam(paramPCcf, ns);

	return paramSamples;
}

// Get moments of parameters given posterior input
void Post::momParam(Array1D<double>& m, Array1D<double>& parMean, Array1D<double>& parVar, bool covFlag, Array2D<double>& parCov)
{

    Array2D<double> paramPCcf=this->getParamPCcf(m);

//	Array2D<double> paramSamples=this->Mrv_->samParam(paramPCcf, ns);
    this->Mrv_->computeMoments(paramPCcf, parMean,parVar,covFlag, parCov);

	return;
}

// Sample forward function at a given grid for given posterior input
Array2D<double> Post::samForwardFcn(Array2D<double> (*forwardFunc)(Array2D<double>&, Array2D<double>&, Array2D<double>&, void*),Array1D<double>& m, Array2D<double>& xgrid, int ns)
{

    Array2D<double> paramPCcf=this->getParamPCcf(m);

    Array2D<double> funcSam=this->Mrv_->propMC(forwardFunc, this->fixIndNom_,this->funcinfo_,paramPCcf, xgrid,ns);


	return funcSam;
}

// Get moments of forward function at a given grid for given posterior input
void Post::momForwardFcn(Array2D<double> (*forwardFunc)(Array2D<double>&, Array2D<double>&, Array2D<double>&, void*),Array1D<double>& m, Array2D<double>& xgrid, Array1D<double>& fcnMean, Array1D<double>& fcnVar, bool covFlag, Array2D<double>& fcnCov)
{

    Array2D<double> paramPCcf=this->getParamPCcf(m);

	Array2D<double> funcCf=this->Mrv_->propNISP(forwardFunc, this->fixIndNom_,this->funcinfo_,paramPCcf, xgrid);

    this->Mrv_->computeMoments(funcCf, fcnMean,fcnVar,covFlag, fcnCov);
	return;
}

// Get moments of composite forward function at a given grid for given posterior input
void Post::momForwardFcn(Array1D<double>& m, Array2D<double>& xgrid, Array1D<double>& fcnMean, Array1D<double>& fcnVar, bool covFlag, Array2D<double>& fcnCov)
{
    Array2D<double> paramPCcf=this->getParamPCcf(m);
    int nxgrid=xgrid.XSize();
    fcnMean.Resize(nxgrid,0.e0);
    fcnVar.Resize(nxgrid,0.e0);
    fcnCov.Resize(nxgrid,nxgrid,0.e0);

    Array1D<double> weights(this->ncat_,1.e0);
    for (int i=0;i<this->ncat_-1;i++){
        weights(i)=m(m.Length()-extraInferredParams_-this->ncat_+1+i);
        weights(this->ncat_-1)-=weights(i);
    }

    Array1D< Array1D<double> > fcnMeans(0);
    Array1D< Array1D<double> > fcnVars(0);
    Array1D< Array2D<double> > fcnCovs(0);

    for (int i=0;i<this->ncat_;i++){
        Array2D<double> funcCf=this->Mrv_->propNISP(this->forwardFcns_(i), this->fixIndNom_,this->funcinfo_,paramPCcf, xgrid);
        Array1D<double> fcnMean_i,fcnVar_i;
        Array2D<double> fcnCov_i;
        this->Mrv_->computeMoments(funcCf, fcnMean_i,fcnVar_i,covFlag, fcnCov_i);
        for (int j=0;j<nxgrid;j++)
            fcnMean(j)+=weights(i)*fcnMean_i(j);

        fcnMeans.PushBack(fcnMean_i);
        fcnVars.PushBack(fcnVar_i);
        fcnCovs.PushBack(fcnCov_i);
    }

    for (int j=0;j<nxgrid;j++){
        for (int i=0;i<this->ncat_;i++){
            fcnVar(j)+=weights(i)*fcnVars(i)(j);
            fcnVar(j)+=weights(i)*fcnMeans(i)(j)*fcnMeans(i)(j);
        }
        fcnVar(j)-=fcnMean(j)*fcnMean(j);
    }

    if (covFlag){
        printf("Post:momForwardFcn(): Covariance computation not implemented for categorical variables. Exiting.");
        exit(1);
    }



    return;
}
/****************************************************************************************/
/****************************************************************************************/
/****************************************************************************************/

// Evaluate full log-likelihood
double Lik_Full::evalLogLik(Array1D<double>& m)
{

    for (int ic=0;ic<m.Length();ic++){
    	if (m(ic)<this->lower_(ic) or m(ic)>this->upper_(ic)){
    		return -1.e+80;
    	}
    }

    double pi=4.*atan(1.);

    int neach = this->nEachs_(0);
	Array2D<double> ydatat(neach, this->nData_);
    for (int ix=0;ix<this->nData_;ix++){
        assert(this->nEachs_(ix)==neach);
        for (int ie=0;ie<neach; ie++){
            ydatat(ie, ix) = this->yData_(ix)(ie);
        }
    }

    Array1D<double> dataSig=this->dataSigma(m(m.Length()-1));

    // TODO assumes one function
    Array2D<double> funcSam=samForwardFcn(this->forwardFcns_(0),m, this->xData_, this->nsam_);

    Array2D<double> dataNoiseSam(this->nsam_,this->nData_);
    generate_normal(dataNoiseSam, time(NULL));
    Array2D<double> tmp;
    Array2D<double> diagsig=diag(dataSig);
    prodAlphaMatMat(dataNoiseSam,diagsig,1.0,tmp);
    addinplace(funcSam,tmp);

    Array1D<double> weight(this->nsam_,1.e0);
    Array1D<double> dens(neach,0.e0);
    Array1D<double> bdw(this->nData_,this->bdw_);
    if (this->bdw_<=0)
        get_opt_KDEbdwth(funcSam,bdw);
    //cout << "BD " << bdw(0) << " " << bdw(1) << endl;
    getPdf_figtree(funcSam,ydatat,bdw,dens, weight);


    double logLik=0.;
    for (int ie=0;ie<neach;ie++){
        //cout << dens(ie) << endl;
        //ldens2(ie)=-(data(ie,0)-mu)*(data(ie,0)-mu)/(2.*std*std)-log(std*sqrt(2*pi)) ;
        if (dens(ie)==0)
            return -1.e80;
        else
            logLik +=log(dens(ie));
    }


	return logLik;
}

/****************************************************************************************/
/****************************************************************************************/
/****************************************************************************************/

// Evaluate marginal log-likelihood
double Lik_Marg::evalLogLik(Array1D<double>& m)
{
    for (int ic=0;ic<m.Length();ic++){
    	if (m(ic)<this->lower_(ic) or m(ic)>this->upper_(ic)){
    		return -1.e+80;
    	}
    }

    double pi=4.*atan(1.);

    Array1D<double> dataSig=this->dataSigma(m(m.Length()-1));
    // TODO Assumes one function
    Array2D<double> funcSam=samForwardFcn(this->forwardFcns_(0),m, this->xData_, this->nsam_);



    Array2D<double> dataNoiseSam(this->nsam_,this->nData_);

    generate_normal(dataNoiseSam, time(NULL));
    Array2D<double> tmp;
    Array2D<double> diagsig=diag(dataSig);
    prodAlphaMatMat(dataNoiseSam,diagsig,1.0,tmp);
    addinplace(funcSam,tmp);

    double logLik=0.;

    for (int ix=0;ix<this->nData_;ix++){
    	Array2D<double> funcSam_ix(this->nsam_,1);
    	Array2D<double> ydatat_ix(this->nEachs_(ix),1);
    	for(int is=0;is<this->nsam_;is++)
            funcSam_ix(is,0)=funcSam(is,ix);
        for(int ie=0;ie<this->nEachs_(ix);ie++)
            ydatat_ix(ie,0)=this->yData_(ix)(ie);
	    Array1D<double> weight(this->nsam_,1.e0);
	    Array1D<double> dens(this->nEachs_(ix),0.e0);
	    Array1D<double> bdw(1,this->bdw_);
        if (this->bdw_<=0)
    	    get_opt_KDEbdwth(funcSam_ix,bdw);
	    getPdf_figtree(funcSam_ix,ydatat_ix,bdw,dens, weight);

        for(int ie=0;ie<this->nEachs_(ix);ie++){
			if (dens(ie)==0)
            	return -1.e80;
        	else
            	logLik +=log(dens(ie));
        }


    }


	return logLik;
}

/****************************************************************************************/
/****************************************************************************************/
/****************************************************************************************/

// Evaluate mvn log-likelihood
double Lik_MVN::evalLogLik(Array1D<double>& m)
{
    for (int ic=0;ic<m.Length();ic++){
    	if (m(ic)<this->lower_(ic) or m(ic)>this->upper_(ic)){
    		return -1.e+80;
    	}
    }

    double pi=4.*atan(1.);

    Array1D<double> dataSig=this->dataSigma(m(m.Length()-1));

 	Array1D<double> fcnMean,fcnVar;
    Array2D<double> fcnCov;
 	momForwardFcn(m, this->xData_, fcnMean, fcnVar, true, fcnCov);


    Array2D<double> fcnCovN=fcnCov;

    for (int i=0;i<this->nData_;i++)
    	fcnCovN(i,i)+=(pow(dataSig(i),2.)+pow(this->nugget_,2.));

    double logLik=evalLogMVN(this->yDatam_,fcnMean,fcnCovN);


	return logLik;
}

/****************************************************************************************/
/****************************************************************************************/
/****************************************************************************************/

// Evaluate gaussian-marginal log-likelihood
double Lik_GausMarg::evalLogLik(Array1D<double>& m)
{
    for (int ic=0;ic<m.Length();ic++){
    	if (m(ic)<this->lower_(ic) or m(ic)>this->upper_(ic)){
    		return -1.e+80;
    	}
    }

    double pi=4.*atan(1.);

    Array1D<double> dataSig=this->dataSigma(m(m.Length()-1));

    Array1D<double> fcnMean,fcnVar;
    Array2D<double> fcnCov;
     // TODO assumes one function
 	momForwardFcn(m, this->xData_, fcnMean, fcnVar, false, fcnCov);

    Array1D<double> dataVar=dotmult(dataSig,dataSig);

    double logLik=0.;
    for (int ix=0;ix<this->nData_;ix++){
        double err=fabs(this->yDatam_(ix)-fcnMean(ix));
        logLik-= ( 0.5*err*err/(fcnVar(ix)+dataVar(ix)) + 0.5*log(2.*pi) + 0.5*log(fcnVar(ix)+dataVar(ix)) );
    }

	return logLik;
}


/****************************************************************************************/
/****************************************************************************************/
/****************************************************************************************/

// Evaluate ABC log-likelihood
double Lik_ABC::evalLogLik(Array1D<double>& m)
{
    for (int ic=0;ic<m.Length();ic++){
    	if (m(ic)<this->lower_(ic) or m(ic)>this->upper_(ic)){
    		return -1.e+80;
    	}
    }

    double pi=4.*atan(1.);

    Array1D<double> dataSig=this->dataSigma(m(m.Length()-1));

    Array1D<double> fcnMean,fcnVar;
    Array2D<double> fcnCov;
 	momForwardFcn(m, this->xData_, fcnMean, fcnVar, false, fcnCov);


    Array1D<double> dataVar=dotmult(dataSig,dataSig);


    ////////
    double alpha=1.;
    double norm=0.0;
    double logLik=0.;
    for (int ix=0;ix<this->nData_;ix++){
        double err=fabs(this->yDatam_(ix)-fcnMean(ix));
        norm+=pow(err,2.);
        norm+=pow(alpha*err-sqrt(fcnVar(ix)+dataVar(ix)),2.);
    }
	logLik=-(0.5/(this->abceps_*this->abceps_))*(norm)-0.5*log(2.*pi)-log(this->abceps_);

	return logLik;
}

/****************************************************************************************/
/****************************************************************************************/
/****************************************************************************************/

// Evaluate ABC-mean log-likelihood
double Lik_ABCm::evalLogLik(Array1D<double>& m)
{
    for (int ic=0;ic<m.Length();ic++){
    	if (m(ic)<this->lower_(ic) or m(ic)>this->upper_(ic)){
    		return -1.e+80;
    	}
    }

    double pi=4.*atan(1.);

    //Array1D<double> dataSig=this->dataSigma(m(m.Length()-1));

    Array1D<double> fcnMean,fcnVar;
    Array2D<double> fcnCov;
 	momForwardFcn(m, this->xData_, fcnMean, fcnVar, false, fcnCov);

    double norm=0.0;
    double logLik=0.;
    for (int ix=0;ix<this->nData_;ix++){
        double err=fabs(this->yDatam_(ix)-fcnMean(ix));
        norm+=pow(err,2.);
    }
	logLik=-(0.5/(this->abceps_*this->abceps_))*(norm)-0.5*log(2.*pi)-log(this->abceps_);

	return logLik;
}

/****************************************************************************************/
/****************************************************************************************/
/****************************************************************************************/

// Evaluate Kennedy-O'Hagan log-likelihood
double Lik_Koh::evalLogLik(Array1D<double>& m)
{
    for (int ic=0;ic<m.Length();ic++){
        if (m(ic)<this->lower_(ic) or m(ic)>this->upper_(ic)){
            return -1.e+80;
        }
    }

    double pi=4.*atan(1.);

    Array1D<double> dataSig=this->dataSigma(m(m.Length()-1));
    double modelSig=m(m.Length()-2);

    Array1D<double> fcnMean,fcnVar;
    Array2D<double> fcnCov;
    momForwardFcn(m, this->xData_, fcnMean, fcnVar, false, fcnCov);


    Array2D<double> fcnCovKoh(this->nData_,this->nData_);

    for (int i=0;i<this->nData_;i++){
        fcnCovKoh(i,i)=pow(modelSig,2.0);
        for (int j=0;j<this->nData_;j++){
            double norm=0.0;
            for (int k=0;k<this->xDim_;k++)
                norm+=pow(xData_(i,k)-xData_(j,k),2.0);
            norm=sqrt(norm);
            fcnCovKoh(i,j)=pow(modelSig,2.0)*exp(-0.5*pow(norm/this->corLength_,2.0));
            fcnCovKoh(j,i)=fcnCovKoh(i,j);
        }
    }


    for (int i=0;i<this->nData_;i++)
        fcnCovKoh(i,i)+=(pow(dataSig(i),2.));

    double logLik=evalLogMVN(this->yDatam_,fcnMean,fcnCovKoh);


    return logLik;
}

/****************************************************************************************/
/****************************************************************************************/
/****************************************************************************************/

// Evaluate classical log-likelihood
double Lik_Classical::evalLogLik(Array1D<double>& m)
{
    for (int ic=0;ic<m.Length();ic++){
    	if (m(ic)<this->lower_(ic) or m(ic)>this->upper_(ic)){
    		return -1.e+80;
    	}
    }

    double pi=4.*atan(1.);
    Array1D<double> dataSig=this->dataSigma(m(m.Length()-1));

    Array1D<double> fcnMean,fcnVar;
    Array2D<double> fcnCov;
 	momForwardFcn(m, this->xData_, fcnMean, fcnVar, false, fcnCov);
 	//or     fcnMean=samForwardFcn(this->forwardFcns_(0),m, this->xData_, 1);// assuming one function


    double logLik=0.;
    for (int ix=0;ix<this->nData_;ix++){
        for (int ie=0;ie<this->nEachs_(ix);ie++){
            double err=fabs(this->yData_(ix)(ie)-fcnMean(ix));
            logLik-= ( 0.5*err*err/(dataSig(ix)*dataSig(ix)) + 0.5*log(2.*pi) + log(dataSig(ix)) );
        }
    }


	return logLik;
}

/****************************************************************************************/
/****************************************************************************************/
/****************************************************************************************/
// Evaluate Error-in-variable log-likelihood
double Lik_Eov::evalLogLik(Array1D<double>& m)
{
    for (int ic=0;ic<m.Length();ic++){
        if (m(ic)<this->lower_(ic) or m(ic)>this->upper_(ic)){
            return -1.e+80;
        }
    }

    Array2D <double> xdata_std(this->nData_,this->xDim_);
    read_datafile(xdata_std,"xdata_std.txt");
    int nord=5;

    double pi=4.*atan(1.);
    Array1D<double> dataSig=this->dataSigma(m(m.Length()-1));

    PCSet PCModel("NISP",nord,this->xDim_,"HG",0.0,1.0);
    int npc=PCModel.GetNumberPCTerms();
    /// Get the default quadrature points
    Array2D<double> qdpts;
    Array1D<double> wghts;
    PCModel.GetQuadPointsWeights(qdpts,wghts);
    int totquad=PCModel.GetNQuadPoints();
    PCModel.SetVerbosity(0);

    Array1D<double> fcnMean(this->nData_),fcnVar(this->nData_);
    for (int ix=0;ix<this->nData_;ix++){
        Array2D<double> func_input(totquad,this->xDim_);
        for (int ixd=0;ixd<this->xDim_;ixd++){
            Array1D<double> cfs_in(npc,0.0);
            cfs_in(0)=this->xData_(ix,ixd);
            cfs_in(ixd+1)=xdata_std(ix,ixd);
            Array1D<double> xch;
            PCModel.EvalPCAtCustPoints(xch, qdpts,cfs_in);
                //cout << cfs_in(0) << " "<< cfs_in(1) << " "<< cfs_in(2) << " "<< cfs_in(3) << endl;
            for (int j=0;j<totquad;j++){

                //xch has length totquad
                func_input(j,ixd)=xch(j);

            }

        }


        Array2D<double> fixindnom(0,0,0.e0);
        int mdim=m.Length();// TODO assumes no noise inference
        Array2D<double> m2d(1,mdim,0.e0);
        for(int im=0;im<mdim;im++)
            m2d(0,im)=m(im);
        void* funcinfo=NULL;
        Array2D<double> func_output=this->forwardFcns_(0)(m2d, func_input, fixindnom, funcinfo);
        //func_output is (1,totquad)

        Array1D<double> func_output_1d(totquad,0.0);
        for (int iq=0;iq<totquad; iq++) func_output_1d(iq)=func_output(0,iq);
        Array1D<double> fk;
        PCModel.GalerkProjection(func_output_1d,fk);
        fcnMean(ix)=PCModel.ComputeMean(fk);

        Array1D<double> varfrac_dummy;
        fcnVar(ix)=PCModel.ComputeVarFrac(fk,varfrac_dummy);
    }



    double logLik=0.;
    for (int ix=0;ix<this->nData_;ix++){
        for (int ie=0;ie<this->nEachs_(ix);ie++){
            double err=fabs(this->yData_(ix)(ie)-fcnMean(ix));

            logLik-= ( 0.5*err*err/(dataSig(ix)*dataSig(ix)+fcnVar(ix)) + 0.5*log(2.*pi) + 0.5*log(dataSig(ix)*dataSig(ix)+fcnVar(ix)) );
        }
    }


    return logLik;
}


