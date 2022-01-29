/* =====================================================================================

                      The UQ Toolkit (UQTk) version 3.1.2
                          Copyright (2022) NTESS
                        https://www.sandia.gov/UQToolkit/
                        https://github.com/sandialabs/UQTk

     Copyright 2022 National Technology & Engineering Solutions of Sandia, LLC (NTESS).
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
/// \file mrv.cpp
/// \author K. Sargsyan  2016 -
/// \brief Multivariate random variable class

#include <math.h>
#include <cfloat>
#include <assert.h>

#include "mrv.h"
#include "PCSet.h"
#include "error_handlers.h"
#include "arrayio.h"
#include "arraytools.h"

#define SIG_MAX 10

// Constructor
Mrv::Mrv(int ndim,string pdfType, Array1D<int> rndInd, int order,string pctype)
{
    // Set the appropriate variables
	this->nDim_=ndim;
	this->pdfType_=pdfType;
	this->rndInd_=rndInd;
	this->rDim_=this->rndInd_.Length();
    this->order_=order;
    this->pcType_=pctype;
    this->pcModel_ = new PCSet("NISP",this->order_,this->rDim_,this->pcType_,0.,1.);
    this->nPC_=this->pcModel_->GetNumberPCTerms();

    // Sanity check and reorder the indices of randomized parameters
    if (this->rDim_>0){
        shell_sort(this->rndInd_);
        assert(this->rndInd_(this->rDim_-1)<this->nDim_);
    }

	return;
}


// Parameterize the multivariate PC RV
// Does the bookkeeping, i.e. maps each parameter alpha to a physical parameter and PC order
int Mrv::Parametrize()
{
    // Counters
    int i=0;
    int pccounter=0;

    // Clear the result container
    this->paramId_.Clear();
    this->pctermId_.Clear();

    // Loop through all randomized parameters and fill in paramId and pctermId
    for (int cur_rparam=0;cur_rparam<this->rDim_;cur_rparam++){
        int rParamInd=this->rndInd_(cur_rparam);
        while(i<rParamInd){
            this->paramId_.PushBack(i);
            this->pctermId_.PushBack(0);
            i++;
        }
        pccounter++;

        // The PC PDF type dictates how the parameterization is done
        if(this->pdfType_=="pct"){
            for (int j=0;j<pccounter+1;j++){
                this->paramId_.PushBack(i);
                this->pctermId_.PushBack(j);
            }
        }
        else if(this->pdfType_=="pci"){
            for (int j=0;j<pccounter+1;j+=pccounter){
                this->paramId_.PushBack(i);
                this->pctermId_.PushBack(j);
            }
        }
        else if(this->pdfType_=="full"){
            for (int j=0;j<this->nPC_;j++){
                this->paramId_.PushBack(i);
                this->pctermId_.PushBack(j);
            }
        }
        i++;
    }

    // Fill the rest
    while (i<this->nDim_){
        this->paramId_.PushBack(i);
        this->pctermId_.PushBack(0);
        i++;
    }

    this->pDim_=this->paramId_.Length();

    printf("\n Mapping of R.V. parameters: \n");
    for (int i=0;i<this->pDim_;i++){
        printf("R.V. parameter %d is alpha_(%d,%d), i.e. Model parameter %d, PC term %d.\n",i,this->paramId_(i),this->pctermId_(i),this->paramId_(i),this->pctermId_(i));
    }

	return this->pDim_;
}

// Get the bounds according to the invariance logic,
// i.e. the last randomized parameter has positive coefficients in higher-order terms
void Mrv::getBounds(Array1D<double>& lower, Array1D<double>& upper)
{
	lower.Resize(this->pDim_,-DBL_MAX);
	upper.Resize(this->pDim_,DBL_MAX);
	if(this->pdfType_=="pct"){

        for (int ic=0;ic<this->pDim_;ic++){
        	if(this->pctermId_(ic)!=0){
        		upper(ic)=SIG_MAX;
	        	if(this->paramId_(ic)==this->rndInd_(this->rDim_-1))
	        		lower(ic)=0.e0;
    	    }
        }

    }
    else if(this->pdfType_=="pci"){

        for (int ic=0;ic<this->pDim_;ic++){
            if(this->pctermId_(ic)!=0){
                upper(ic)=SIG_MAX;
                lower(ic)=0.e0;
            }
        }

    }
    else if(this->pdfType_=="full"){
        // No bounds for 'full' representation
    }

	return;
}

// Get matrix of PC coefficients given parameterization
// The matrix is convenient to work with and has dimensions npc-by-ndim
Array2D<double> Mrv::getMultiPCcf(Array1D<double>& rvParams)
{
	assert(rvParams.Length()==this->pDim_);
	Array2D<double> multiPCcf(this->nPC_,this->nDim_,0.e0);

    for (int ic=0;ic<this->pDim_;ic++)
        multiPCcf(this->pctermId_(ic),this->paramId_(ic))=rvParams(ic);

    return multiPCcf;
}

// Evaluate multivariate PC given germ samples and coefficient matrix
Array2D<double> Mrv::evalMultiPC(Array2D<double>& xiSam, Array2D<double>& multiPCcf)
{
	int nsam=xiSam.XSize();
	assert(xiSam.YSize()==this->rDim_);

    Array2D<double> multiPCSam(nsam,0);
    for (int j=0;j<this->nDim_;j++){
        Array1D<double> pdfpccf;
        getCol(multiPCcf,j,pdfpccf);

        Array1D<double> samples_dim;
        this->pcModel_->EvalPCAtCustPoints(samples_dim,xiSam, pdfpccf);
        multiPCSam.insertCol(samples_dim,j);
    }

	return multiPCSam;
}

// Random sample of all parameters given a coefficient matrix
Array2D<double> Mrv::mcParam(Array2D<double>& multiPCcf, int nsam)
{
	Array2D<double> xiSam(nsam,this->rDim_,0.e0);
	this->pcModel_->DrawSampleVar(xiSam);
	Array2D<double> multiPCSam=this->evalMultiPC(xiSam,multiPCcf);

	return multiPCSam;
}

// Get quadrature samples of all parameters given coefficient matrix
Array2D<double> Mrv::quadParam(Array2D<double>& multiPCcf)
{
	Array2D<double> xiSam;
    Array1D<double> weights;

    this->pcModel_->GetQuadPointsWeights(xiSam, weights);
    Array2D<double> multiPCQuadSam=this->evalMultiPC(xiSam,multiPCcf);

	return multiPCQuadSam;
}

// Propagate the multivariate RV with given coefficeints through a given function at given values x
Array2D<double> Mrv::propNISP(Array2D<double> (*forwardFcn)(Array2D<double>&, Array2D<double>&, Array2D<double>&, void*), Array2D<double>& fixindnom,void* funcinfo, Array2D<double>& multiPCcf, Array2D<double>& x)
{
	Array2D<double> multiPCQuadSam=this->quadParam(multiPCcf);
    int ns=multiPCQuadSam.XSize();

	Array2D<double> funcQuad=forwardFcn(multiPCQuadSam,x, fixindnom, funcinfo);
	int nx=funcQuad.YSize();

	Array2D<double> funcCf(nPC_,0);

    for (int ix=0;ix<nx;ix++){
        Array1D<double> funcQuad_ix;
        getCol(funcQuad,ix,funcQuad_ix);
        Array1D<double> fcf;
        this->pcModel_->GalerkProjection(funcQuad_ix,fcf);
        funcCf.insertCol(fcf,ix);
    }

	return funcCf;
}

// Sample values of a given function given input coefficeint matrix
Array2D<double> Mrv::propMC(Array2D<double> (*forwardFcn)(Array2D<double>&, Array2D<double>&, Array2D<double>&, void*), Array2D<double>& fixindnom,void* funcinfo,Array2D<double>& multiPCcf, Array2D<double>& x,int nsam)
{
	Array2D<double> multiPCSam=this->mcParam(multiPCcf, nsam);
    int ns=multiPCSam.XSize();


	Array2D<double> funcSam=forwardFcn(multiPCSam,x, fixindnom, funcinfo);

	return funcSam;
}

// Compute moments given coefficent matrix
void Mrv::computeMoments(Array2D<double>& funcCf, Array1D<double>& fcnMean,Array1D<double>& fcnVar,bool covFlag, Array2D<double>& fcnCov)
{
	int nx=funcCf.YSize();
    fcnMean.Resize(nx);
    fcnVar.Resize(nx);
    if (covFlag)
    	fcnCov.Resize(nx,nx,0.e0);

    for (int ix=0;ix<nx;ix++){
        Array1D<double> fcf;
        getCol(funcCf,ix,fcf);
        fcnMean(ix)=this->pcModel_->ComputeMean(fcf);

        Array1D<double> varfrac;
        double var=this->pcModel_->ComputeVarFrac(fcf,varfrac);
        fcnVar(ix)=var;
        if (covFlag){
        	fcnCov(ix,ix)=var;

	        for (int jx=0;jx<ix;jx++){
	            Array1D<double> fcf_;
	            getCol(funcCf,jx,fcf_);
	            Array1D<double> varfrac_;
	            double var_=this->pcModel_->ComputeVarFrac(fcf_,varfrac_);
	            double cov=0.0;
	            for (int k=1;k<varfrac_.XSize();k++)
	                cov+=sqrt(var*varfrac(k)*var_*varfrac_(k));

	            fcnCov(ix,jx)=cov;//*0.99;
	            fcnCov(jx,ix)=cov;//*0.99;
	        }
	    }


    }

	return;
}





