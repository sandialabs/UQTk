/* =====================================================================================

                      The UQ Toolkit (UQTk) version 3.1.4
                          Copyright (2023) NTESS
                        https://www.sandia.gov/UQToolkit/
                        https://github.com/sandialabs/UQTk

     Copyright 2023 National Technology & Engineering Solutions of Sandia, LLC (NTESS).
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
/// \file lreg.cpp
/// \author K. Sargsyan  2015 -
/// \brief Linear regression class


#include <math.h>
#include <cfloat>
#include <iostream>

#include "lreg.h"
#include "gen_defs.h"
#include "error_handlers.h"

#include "arraytools.h"
#include "arrayio.h"
#include "PCSet.h"

#include "bcs.h"

#include "ftndefs.h"
#include "depblas.h"
#include "deplapack.h"


#include "tools.h"
#include "lbfgs_routines.h"

#include <assert.h>
//#define DEBUG

/// Leave-one-out error computation done as a global function for optimization
/// \todo Find a more elegant way to do this within the class
double loo(int ndim, double* m, void* classpointer);


// Constructor for Radial Basis Function basis class
RBFreg::RBFreg(Array2D<double>& centers, Array1D<double>& widths)
{
    centers_=centers;
    widths_=widths;
    ndim_=centers.YSize();
    nbas_=centers.XSize();
    CHECKEQ(ndim_,widths.XSize());

    return;
}

// Constructor for Polynomial Chaos basis class, given order and dim
PCreg::PCreg(string strpar,int order, int dim)
{
    Array2D<int> mindex;
    computeMultiIndex(dim,order,mindex);

    this->SetMindex(mindex);
    pctype_=strpar;

    nbas_=this->mindex_.XSize();
    ndim_=this->mindex_.YSize();

    return;
}

// Constructor for Polynomial Chaos basis class, given multiindex
PCreg::PCreg(string strpar,Array2D<int>& mindex)
{
    this->SetMindex(mindex);
    pctype_=strpar;

    nbas_=this->mindex_.XSize();
    ndim_=this->mindex_.YSize();


    return;
}

// Constructor of monomial basis class, given order and dim
PLreg::PLreg(int order, int dim)
{
    Array2D<int> mindex;
    computeMultiIndex(dim,order,mindex);

    this->SetMindex(mindex);

    nbas_=this->mindex_.XSize();
    ndim_=this->mindex_.YSize();

    return;
}

// Constructor of monomial basis class, given multiindex
PLreg::PLreg(Array2D<int>& mindex)
{
    this->SetMindex(mindex);

    nbas_=this->mindex_.XSize();
    ndim_=this->mindex_.YSize();

    return;
}

/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////

// Initializing linear regression, parent class
void Lreg::InitRegr()
{
    weights_.Resize(this->nbas_,0.e0);
    dataSetFlag_=false;
    regMode_="m";
}

// void Lreg::testing()
// {
//     Array1D<double> mm;
//     mm.Resize(1,0.4);
//     cout << loo(1, mm.GetArrayPointer(), this) << " CCCC " << endl;
//     mm.Resize(1,0.2);
//     cout << loo(1, mm.GetArrayPointer(), this) << " CCCC " << endl;

// }

// Set the x- and y(1D)- data
void Lreg::SetupData(Array2D<double>& xdata, Array1D<double>& ydata)
{
    xdata_=xdata;
    ydata_=ydata;
    npt_=xdata.XSize();
    CHECKEQ(npt_,ydata.XSize());
    CHECKEQ(ndim_,xdata.YSize());

    dataSetFlag_=true;

    return;
}

// Set the x- and y(2D)- data
void Lreg::SetupData(Array2D<double>& xdata, Array2D<double>& ydata)
{

    int nx=xdata.XSize();
    CHECKEQ(ndim_,xdata.YSize());
    CHECKEQ(nx,ydata.XSize());
    int numeach=ydata.YSize();


    npt_=nx*numeach;

    xdata_.Resize(npt_,ndim_);
    ydata_.Resize(npt_);


    int j=0;
    for(int i=0; i<nx; i++){
        for(int ie=0;ie<numeach;ie++){
            for(int id=0;id<ndim_;id++){
                xdata_(j,id)=xdata(i,id);
            }
            ydata_(j)=ydata(i,ie);
            j++;
        }
    }

    dataSetFlag_=true;

    return;
}

// Setting regression weights
void Lreg::SetRegWeights(Array1D<double>& weights)
{
    weights_=weights;
    CHECKEQ(this->nbas_,weights_.XSize());

  return;
}

// Building regression with BCS methods
void Lreg::BCS_BuildRegr(Array1D<int>& used, double eta)
{
    assert(dataSetFlag_);
    int ntot=ydata_.XSize();
    // Work variables
    Array1D<double> lambda_init,weights, errbars, basis, alpha ;
    Array2D<int> newmindex;
    Array2D<double> Sig;

    // Hardwired parameters
    int    adaptive = 0      ;
    int    optimal  = 1      ;
    double scale    = 0.1    ;
    int    verbose  = 0      ;
    // Initial variance 'rule-of-thumb'
    sigma2_(0) = max(1.e-12,get_var(ydata_)/1.0e6);

    // Compute the basis evaluation matrix
    this->EvalBases(xdata_,bdata_);

    // Run the BCS
    lambda_init=weights_;
    WBCS(bdata_,ydata_,sigma2_,eta,lambda_init,adaptive, optimal,scale,
                 verbose,coef_,used,coef_erb_,basis,alpha,coef_cov_);

    // Only retain the bases selected by WBCS
    this->StripBases(used);

    return;
}

// Building least-squares regression
void Lreg::LSQ_BuildRegr()
{
    residFlag_=false;
    diagPFlag_=false;

    // Make sure data is set
    assert(dataSetFlag_);

    // Fill in the basis evaluation matrix
    this->EvalBases(xdata_,bdata_);
    A_=MatTMat(bdata_) ;

    // Add the regularization weights to the diagonal
    for (int i=0;i<nbas_;i++){
        A_(i,i)+=weights_(i);
    }

    // Matrix inversion
    this->A_inv_ = INV(this->A_);

    // Compute via the classical least-square formula
    prodAlphaMatTVec(bdata_, ydata_, 1.0, Hty_) ;
    prodAlphaMatVec(A_inv_, Hty_, 1.0, coef_) ;


    // If errorbars are requested
    if (this->regMode_=="ms" or this->regMode_=="msc"){
        double s1=0.0,s2=0.0;

        for(int it=0;it<npt_;it++)
          s1 += ( ydata_(it)*ydata_(it) );
        for(int ip=0;ip<nbas_;ip++)
          s2 += ( Hty_(ip)*coef_(ip) );
        double betta_add=(s1-s2)/2.;



        // alpha=beta=0 for Jeffreys prior 1/sigma^2
        double alfa=0.0;
        double betta=0.0;
        // alpha>>1, beta>>1, beta/alpha=sigma for fixed sigma
        // \todo to be implemented

        sigma2_(0)=(betta+betta_add)/(alfa+0.5*(npt_-nbas_)-1.);
        if (sigma2_(0)<0.0){
            cout << "Negative (should be very small) data noise, set to zero. Sigma2=" << sigma2_(0) << endl;
            sigma2_(0)=0.0;
        }
        coef_cov_.Resize(nbas_,nbas_,0.e0);
        for(int ib=0;ib<nbas_;ib++)
          for(int jb=0;jb<nbas_;jb++)
            coef_cov_(ib,jb)=sigma2_(0)*A_inv_(ib,jb);//-1 due to eta/eta-2 factor in cov

        coef_erb_.Resize(nbas_,0.e0);
        for ( int i = 0; i<this->nbas_; i++) coef_erb_(i) = sqrt(coef_cov_(i,i));

    }

  return;
}


// Evaluate the pre-built regression at given values
void Lreg::EvalRegr(Array2D<double>& xcheck, Array1D<double>& ycheck,Array1D<double>& yvar,Array2D<double>& ycov)
{
    CHECKEQ(ndim_,xcheck.YSize());
    int ncheck=xcheck.XSize();
    Array2D<double> bcheck;

    this->EvalBases(xcheck,bcheck);
    prodAlphaMatVec(bcheck, coef_, 1.0, ycheck) ;

    // If more than the the means are requested
    if (this->regMode_=="ms" or this->regMode_=="msc"){
        if (sigma2_(0)==0.0){
            yvar.Resize(ncheck,0.0);
            ycov.Resize(ncheck,ncheck,0.0);
            return;
        }

        Array2D<double> L=coef_cov_;

        // Cholesky factorization of the covariance
        int nd=L.XSize();
        CHECKEQ(nd,this->nbas_);
        int chol_info=0;
        char lu='L';
        FTN_NAME(dpotrf)(&lu,&nd, L.GetArrayPointer(),&nd,&chol_info);

        // Catch the error in Cholesky factorization
        if (chol_info != 0 ){
            printf("Lreg::EvalRegr():Error in Cholesky factorization, info=%d\n", chol_info);
            exit(1);
        }

        for (int i=0;i<nd;i++)
            for (int j=i+1;j<nd;j++)
                L(i,j)=0.0;

        Array2D<double> A;
        prodAlphaMatMat(bcheck,L,1.0,A);

        yvar.Resize(ncheck);
        for (int i=0;i<ncheck;i++){
            double sum=0.0;
            for (int k=0;k<nd;k++)
                sum+=A(i,k)*A(i,k);
            yvar(i)=sum;

            // If covariances are requested
            if (this->regMode_=="msc"){

                ycov.Resize(ncheck,ncheck);
                ycov(i,i)=yvar(i);
                for (int j=0;j<i;j++){
                 double sum=0.0;
                    for (int k=0;k<ncheck;k++)
                        sum+=A(i,k)*A(j,k);
                    ycov(i,j)=sum;
                    ycov(j,i)=sum;
                }
            }
        }
    }

  return;
}


// Projection to the space spanned by the bases
void Lreg::Proj(Array1D<double>& array,Array1D<double>& proj_array)
{
    CHECKEQ(npt_,array.XSize());
    Array1D<double> tmp,tmp2;
    prodAlphaMatTVec(bdata_,array,1.0,tmp);
    prodAlphaMatVec(A_inv_,tmp,1.0,tmp2);
    prodAlphaMatVec(bdata_,tmp2,1.0,proj_array);

    for (int i=0;i<proj_array.XSize();i++){
        proj_array(i)=array(i)-proj_array(i);
    }

    return;
}

// Compute best values for regularization weight vector
Array1D<double> Lreg::LSQ_computeBestLambdas()
{

    Array1D<double> lambda(nbas_,0.1);
    int n=nbas_; //1;//ndim_;
    int m=5;
    Array1D<int> nbd(n,1);
    Array1D<double> l(n,0.e0);
    Array1D<double> u(n,0.e0);

    void* info=this;

    // Minimize leave-one-out error
    lbfgsDR(n,m,lambda.GetArrayPointer(),nbd.GetArrayPointer(),l.GetArrayPointer(),u.GetArrayPointer(),loo,NULL,info) ;


    return lambda;
}

// Compute best values for regularization weight parameter
double Lreg::LSQ_computeBestLambda()
{

    double lambda=0.1;
    int n=1;//ndim_;
    int m=5;
    Array1D<int> nbd(n,1);
    Array1D<double> l(n,0.e0);
    Array1D<double> u(n,0.e0);

    void* info=this;

    // Minimize leave-one-out error
    lbfgsDR(n,m,&lambda,nbd.GetArrayPointer(),l.GetArrayPointer(),u.GetArrayPointer(),loo,NULL,info) ;

    return lambda;
}

// Compute residual
void Lreg::getResid()
{

    if (!residFlag_){
       Array1D<double> tmp;
        prodAlphaMatVec(bdata_,coef_,1.0,tmp);
        resid_=subtract(ydata_,tmp);
        residFlag_=true;
    }
    return;
}

// COmpute the diagonal of the projection matrix
void Lreg::getDiagP()
{
    if (!diagPFlag_){
        diagP_.Resize(npt_,1.e0);
        for (int i=0;i<npt_;i++){
            for (int j=0;j<nbas_;j++){
                for (int k=0;k<j;k++)
                    diagP_(i)-=2.*bdata_(i,j)*A_inv_(j,k)*bdata_(i,k);
                diagP_(i)-=bdata_(i,j)*A_inv_(j,j)*bdata_(i,j);
            }
        }
        diagPFlag_=true;
    }

    return;
}

// Compute a few error metrics, for least-squares only
Array1D<double> Lreg::computeErrorMetrics(string method)
{
    Array1D<double> err(2,-999.);
    if (method=="lsq"){
        err(0)=this->LSQ_computeLOO();
        err(1)=this->LSQ_computeGCV();
    }
    else{
        printf("Computation of errors not implemented for %s method\n",method.c_str());
    }

    return err;
}


// Compute the validation error given a set of x-y data
double Lreg::computeRVE(Array2D<double>& xval,Array1D<double>& yval,Array1D<double>& yval_regr)
{
    int nval=xval.XSize();
    CHECKEQ(ndim_,xval.YSize());
    CHECKEQ(nval,yval.XSize());

    double sum=0.0;
    Array1D<double> dummy_var;     Array2D<double> dummy_cov;
    this->SetRegMode("m");
    this->EvalRegr(xval,yval_regr, dummy_var,dummy_cov);
    for (int i=0;i<nval;i++){
        double err=yval(i)-yval_regr(i);
        sum += (err*err);
    }
    sum /= nval;
    sum=sqrt(sum);

    return sum;
}


// Compute the leave-one-out error
double Lreg::LSQ_computeLOO()
{
    this->getResid();

    this->getDiagP();

    Array1D<double>   resid_scaled=dotdivide(resid_,diagP_);
    return norm(resid_scaled)/sqrt(npt_);

}

// Compute the generalized cross validation error
double Lreg::LSQ_computeGCV()
{
    this->getResid();
    this->getDiagP();

    double sum=0.0;
    for (int i=0;i<npt_;i++)
        sum+=diagP_(i);

    return  norm(resid_)*sqrt(npt_)/sum;

}


/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////

// Setting RBF centers
void RBFreg::SetCenters(Array2D<double>& centers)
{
    centers_=centers;
    CHECKEQ(ndim_,centers.YSize());
    CHECKEQ(nbas_,centers.XSize());
    return;
}

// Setting RBF widths
void RBFreg::SetWidths(Array1D<double>& widths)
{
    widths_=widths;
    CHECKEQ(ndim_,widths.XSize());
    return;
}

// Evaluate PC bases
void PCreg::EvalBases(Array2D<double>& xx,Array2D<double>& bb)
{
    CHECKEQ(ndim_,xx.YSize());

    PCSet currPCModel("NISPnoq",mindex_,pctype_,0.,1.);

    currPCModel.EvalBasisAtCustPts(xx,bb);

    return;

}

// Evaluate monomial bases
void PLreg::EvalBases(Array2D<double>& xx,Array2D<double>& bb)
{
    CHECKEQ(ndim_,xx.YSize());
    int nx=xx.XSize();
    bb.Resize(nx,nbas_,0.e0);

    for(int is=0;is<nx;is++){
        for(int ipc=0;ipc<nbas_;ipc++){
            bb(is,ipc)=1;
            for(int id=0;id<ndim_;id++){
                bb(is,ipc)*=pow(xx(is,id),mindex_(ipc,id));
            }
        }
    }

    return;

}


// Evaluate RBFs
void RBFreg::EvalBases(Array2D<double>& xx,Array2D<double>& bb)
{
    CHECKEQ(ndim_,xx.YSize());
    int nx=xx.XSize();

    bb.Resize(nx,nbas_,0.e0);
    for(int is=0;is<nx;is++){
        for(int ib=0;ib<nbas_;ib++){
            double norm2=0.0;
            for(int id=0;id<ndim_;id++){
                double diff=(xx(is,id)-centers_(ib,id))/widths_(id);
                norm2+=diff*diff;
            }
            bb(is,ib)=exp(-norm2/2.);
        }
   }

    return;

}

// Select given PC bases
void PCreg::StripBases(Array1D<int>& used)
{
    int nused=used.Length();
    Array2D<int> mindex_new(nused,ndim_,0);
	for (int i=0;i<nused;i++)
        for (int j=0;j<ndim_;j++)
            mindex_new(i,j)=mindex_(used(i),j);

    SetMindex(mindex_new);
    this->nbas_=nused;
    return;
}

// Select given monomial bases
void PLreg::StripBases(Array1D<int>& used)
{
    int nused=used.Length();
    Array2D<int> mindex_new(nused,ndim_,0);
	for (int i=0;i<nused;i++)
        for (int j=0;j<ndim_;j++)
            mindex_new(i,j)=mindex_(used(i),j);

    SetMindex(mindex_new);
    this->nbas_=nused;
    return;
}

// Select given RBFs
void RBFreg::StripBases(Array1D<int>& used)
{
    int nused=used.Length();
    Array2D<double> centers_new(nused,ndim_,0.);
    Array1D<double> widths_new(nused,0.);

	for (int i=0;i<nused;i++){
        for (int j=0;j<ndim_;j++)
            centers_new(i,j)=centers_(used(i),j);
        widths_new(i)=widths_(used(i));
    }
    SetCenters(centers_new);
    SetWidths(widths_new);
    this->nbas_=nused;

    return;
}


// Leave-one-out error estimator
double loo(int ndim, double* mm, void* classpointer)
{
    Lreg* thisClass=(Lreg*) classpointer;

    Array1D<double> lam(thisClass->GetNbas());
    if (ndim==1)
        lam.Resize(thisClass->GetNbas(),mm[0]);
    else{
        CHECKEQ(ndim,thisClass->GetNbas());
        for (int i=0;i<ndim;i++){
            lam(i)=mm[i];
        }
    }


    thisClass->SetRegWeights(lam);
    thisClass->LSQ_BuildRegr();

    Array1D<double> errors=thisClass->computeErrorMetrics("lsq");

    double err_loo=errors(0);
    return err_loo;

}
