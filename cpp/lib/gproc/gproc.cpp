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

     Questions? Contact the UQTk Developers at <uqtk-developers@software.sandia.gov>
     Sandia National Laboratories, Livermore, CA, USA
===================================================================================== */
/// \file gproc.cpp 
/// \author K. Sargsyan  2014 - 
/// \brief Gaussian Process class

#include <math.h>
#include <cfloat>
#include <iostream>

#include "gproc.h"
#include "error_handlers.h"
#include "gen_defs.h"
#include "arraytools.h"
#include "arrayio.h"

#include "tools.h"
#include "lbfgs_routines.h"

#include <assert.h>

/// Function to compute negative log posterior (needed to maximize with respect to roughness parameter)
/// \todo Find a more elegant way to do this within the class
double neglogPostParam(int ndim, double* m, void* classpointer);

// Constructor given PC basis and covariance type and parameters
Gproc::Gproc(const string covtype, PCSet *PCModel, Array1D<double>& param)
{
  // Set the PC regression
  PCModel_=PCModel;
  npc_=PCModel_->GetNumberPCTerms();  
  printf("gproc : Number of PC terms : %d\n",npc_);


  ndim_=PCModel_->GetNDim();
  param_=param;
  covType_=covtype;

  return;
}

// Setup the prior - currently hardwired to the common practice values
void Gproc::SetupPrior()
{
  // Inverse prior on regression coefficients is set to zero
  Vinv_.Resize(npc_,npc_,0.e0); //V^-1
  // Prior mean on coefficients is set to zero
  z_.Resize(npc_,0.e0); //z
  // Parameters of data noise prior: setting to zero corresponds to 1/sigma prior
  al_=0.e0;
  be_=0.e0;

  // If data noise is explicitely given, use this construction below
  //  sig2f_=0.0;
  //  if (sig2f_!=0.0)
  //    be_=al_*sig2f_;

  return;
}

// Setup x- and y- data, together with the variance
void Gproc::SetupData(Array2D<double>& xdata, Array1D<double>& ydata,Array1D<double>& datavar)
{
    xdata_=xdata;
    ydata_=ydata;
    npt_=xdata.XSize();
    assert(npt_=datavar.XSize());
    dataVar_=datavar;

    return;
}

// Building Gaussian Process 
void Gproc::BuildGP()
{
  cout << "Building GP" << endl;
  int npts=xdata_.XSize();
  int ndim=xdata_.YSize();
  assert(ndim==ndim_);

  // Fill in the projection matrix H_
  printf("Computing projection matrix\n");
  PCModel_->EvalBasisAtCustPts(xdata_,H_);
  transpose(H_,Ht_);

  // Fill in the covariance matrix A_
  printf("Filling in data covariance matrix of size %d x %d\n", npts,npts);
  computeDataCov_(xdata_,param_,A_);
  printf("Done.\n");

  printf("Solving Ax=y for A of size %d x %d, and y of size %d\n", npts,npts,npts);
  // Matrix manipulations to arrive to the answer(this should be done more efficiently)
  Ainvd_=Ainvb(A_,ydata_);
  printf("Done.\n");

  prodAlphaMatVec(Ht_,Ainvd_, 1.0, HtAinvd_) ;

  printf("Solving AX=H for A of size %d x %d, and H of size %d x %d\n", npts,npts,npts,H_.YSize());
  AinvH_=AinvH(A_, H_);
  printf("Done.\n");

  // Further linear algebra
  prodAlphaMatMat(Ht_,AinvH_, 1.0, HtAinvH_) ;
  Array2D<double> tmp;
  tmp=add(Vinv_,HtAinvH_);
  Vst_=INV(tmp);
  Array1D<double> temp;
  prodAlphaMatVec(Vinv_, z_, 1.0, Vinvz_) ;
  temp=add(Vinvz_,HtAinvd_);
  prodAlphaMatVec(Vst_, temp, 1.0, bhat_) ;
  prodAlphaMatVec(H_, bhat_, 1.0, Hbhat_) ;
  yHbhat_=subtract(ydata_,Hbhat_);

  AinvyHbhat_=Ainvb(A_,yHbhat_);

  // Get the sigma_^2_hat, i.e. best value for the data noise
  sig2hat_=2.*be_;
  sig2hat_+=prod_vecTmatvec(z_,Vinv_,z_);
  //sig2hat_+=prod_vecTmatvec(ydata_,Ainv_,ydata_);
  sig2hat_+=dot(ydata_,Ainvd_);

  Vstinv_=INV(Vst_);
  sig2hat_-=prod_vecTmatvec(bhat_,Vstinv_,bhat_);
  sig2hat_ /= (npts+2.*al_-npc_-2.);

  return;
}

// An older implementation that relies on explicit matrix inversion (is kept for further timing tests)
void Gproc::BuildGP_inv()
{
  cout << "Building GP" << endl;
  int npts=xdata_.XSize();
  int ndim=xdata_.YSize();
  assert(ndim==ndim_);

  // Fill in the projection matrix H_
  printf("Computing projection matrix\n");
  PCModel_->EvalBasisAtCustPts(xdata_,H_);
  transpose(H_,Ht_);

  // Fill in the covariance matrix A_
  printf("Computing data covariance matrix\n");
  computeDataCov_(xdata_,param_,A_);


  cout << "gproc : Computing the inverse of the covariance matrix at the data points...." << endl;
  Ainv_=INV(A_);
  cout << "gproc : Done computing the inverse." << endl;


  // Matrix manipulations to arrive to the answer (this should be done more efficiently)
  prodAlphaMatVec(Ainv_, ydata_, 1.0, Ainvd_) ;
  prodAlphaMatVec(Ht_,Ainvd_, 1.0, HtAinvd_) ;
  prodAlphaMatMat(Ainv_,H_, 1.0, AinvH_) ;
  prodAlphaMatMat(Ht_,AinvH_, 1.0, HtAinvH_) ;
  Array2D<double> tmp;
  tmp=add(Vinv_,HtAinvH_);
  Vst_=INV(tmp);
  Array1D<double> temp;
  prodAlphaMatVec(Vinv_, z_, 1.0, Vinvz_) ;
  temp=add(Vinvz_,HtAinvd_);
  prodAlphaMatVec(Vst_, temp, 1.0, bhat_) ;
  prodAlphaMatVec(H_, bhat_, 1.0, Hbhat_) ;
  yHbhat_=subtract(ydata_,Hbhat_);
  prodAlphaMatVec(Ainv_, yHbhat_, 1.0, AinvyHbhat_) ;
 
  // Get the sigma_^2_hat
  sig2hat_=2.*be_;
  sig2hat_+=prod_vecTmatvec(z_,Vinv_,z_);
  sig2hat_+=prod_vecTmatvec(ydata_,Ainv_,ydata_);
  Vstinv_=INV(Vst_);
  sig2hat_-=prod_vecTmatvec(bhat_,Vstinv_,bhat_);
  sig2hat_ /= (npts+2.*al_-npc_-2.);

  return;
}

// Evaluate GP after it is built, given x-data
void Gproc::EvalGP(Array2D<double>& xgrid, string msc, Array1D<double>& mst)
{
  assert(ndim_==xgrid.YSize());
  int totgrid=xgrid.XSize();
  int npts=xdata_.XSize();

  // Get the mean values at the grid points
  PCModel_->EvalPCAtCustPoints(mst,xgrid,bhat_);


  Array2D<double> ttmat(totgrid,npts);
  printf("Evaluating GP mean... \n");
  for (int igr=0;igr<totgrid;igr++){
    if ((igr+1)%10000==0)
      printf(" at the %d / %d point\n", igr+1,totgrid);

    Array1D<double> xcurr;
    getRow(xgrid, igr, xcurr);
    for (int ipts=0;ipts<npts;ipts++){
      Array1D<double> xdata_i;
      getRow(xdata_, ipts, xdata_i);
      ttmat(igr,ipts)=covariance(xcurr,xdata_i,param_);
    }
  }

  Array1D<double> mst_corr;
  prodAlphaMatVec(ttmat, AinvyHbhat_, 1.0, mst_corr);
  addinplace(mst,mst_corr);

  // If standard deviation is also requested
  if (msc=="ms"){
    var_.Resize(totgrid,0.e0);
    Array2D<double> ttmat_t;
    transpose(ttmat,ttmat_t);
    Array2D<double> tmp=AinvH(A_,ttmat_t);

    printf("Evaluating GP variance... \n");
    for(int it=0;it<totgrid;it++){
      if ((it+1)%1000==0)
        printf(" at the %d / %d point\n", it+1,totgrid);

      Array1D<double> xcurr;
      getRow(xgrid, it, xcurr);

      Array1D<double> tt;
      getRow(ttmat,it, tt);

      Array1D<double> tmpp;
      getCol(tmp,it,tmpp);
      double correction1=dot(tmpp,tt);

      Array1D<double> ttAinvH;
      prodAlphaMatTVec(AinvH_,tt,1.0,ttAinvH);
  
      // Fill in appropriate matrices
      Array2D<double> H_atgrid;
      PCModel_->EvalBasisAtCustPts(xgrid,H_atgrid);

      Array1D<double> ht;
      getRow(H_atgrid,it,ht);
  
      Array1D<double> httAinvH;
      httAinvH=subtract(ht,ttAinvH);

      double correction2=prod_vecTmatvec(httAinvH,Vst_,httAinvH);

    
      var_(it)=(sig2hat_*(covariance(xcurr,xcurr,param_)-correction1+correction2));

    }
  }
  
  // If covariances are also requested
  else  if (msc=="msc"){
    Array2D<double> ttmat_t;
    transpose(ttmat,ttmat_t);
    Array2D<double> tmp=AinvH(A_,ttmat_t);
    Array2D<double> cov_corr1;
    prodAlphaMatMat(ttmat, tmp, 1.0, cov_corr1); 

    // Fill in appropriate matrices
    Array2D<double> H_atgrid;
    PCModel_->EvalBasisAtCustPts(xgrid,H_atgrid);

    Array2D<double> ttmatAinvH;
    prodAlphaMatMat(ttmat,AinvH_,1.0,ttmatAinvH);

    Array2D<double> h_ttmatAinvH;
    h_ttmatAinvH=subtract(H_atgrid,ttmatAinvH);
    Array2D<double> h_ttmatAinvH_t;
    transpose(h_ttmatAinvH,h_ttmatAinvH_t);

    Array2D<double> tmpp;
    prodAlphaMatMat(Vst_,h_ttmatAinvH_t,1.0,tmpp);
    Array2D<double> cov_corr2;
    prodAlphaMatMat(h_ttmatAinvH,tmpp,1.0, cov_corr2);
  
    cov_.Resize(totgrid,totgrid,0.e0);
    var_.Resize(totgrid,0.e0);

    for(int it=0;it<totgrid;it++){
      Array1D<double> xcurr;
      getRow(xgrid, it, xcurr);

      for(int jt=0;jt<totgrid;jt++){
        Array1D<double> xcurra;
        getRow(xgrid, jt, xcurra);
        cov_(it,jt)=covariance(xcurr,xcurra,param_) ;
      }
    }
    subtractinplace(cov_,cov_corr1);
    addinplace(cov_,cov_corr2);
    scaleinplace(cov_,sig2hat_);

    for(int it=0;it<totgrid;it++)
      var_(it)=cov_(it,it);

  }

  return;
}

// An older implementation that relies on explicit matrix inversion (is kept for further timing tests)
void Gproc::EvalGP_inv(Array2D<double>& xgrid, string msc, Array1D<double>& mst)
{
  assert(ndim_==xgrid.YSize());
  int totgrid=xgrid.XSize();
  int npts=xdata_.XSize();

  // Get the mean values at the grid points
  PCModel_->EvalPCAtCustPoints(mst,xgrid,bhat_);

  // Fill in appropriate matrices
  Array2D<double> H_atgrid;
  PCModel_->EvalBasisAtCustPts(xgrid,H_atgrid);

  // If more than the mean is requested
  if (msc != "m"){
    cov_.Resize(totgrid,totgrid,0.e0);
    var_.Resize(totgrid,0.e0);
  }
  // Compute the covariance structure of the student-t process
  for(int it=0;it<totgrid;it++){
    Array1D<double> xcurr(ndim_,0.e0);
    for (int id=0;id<ndim_;id++)
      xcurr(id)=xgrid(it,id);

    Array1D<double> tt(npts,0.e0);
    for(int ipts=0;ipts<npts;ipts++){
      Array1D<double> xdata_i(ndim_,0.e0);
      for (int id=0;id<ndim_;id++)
        xdata_i(id)=xdata_(ipts,id);
      tt(ipts)=covariance(xcurr,xdata_i,param_);
    }
    double correction=dot(tt,AinvyHbhat_);
    mst(it)+=correction;

    if (msc != "m"){
      Array1D<double> ttAinvH;
      prodAlphaMatTVec(AinvH_,tt,1.0,ttAinvH);
    
      Array1D<double> ht(npc_,0.e0);
      for(int ipc=0;ipc<npc_;ipc++)
        ht(ipc)=H_atgrid(it,ipc);

      Array1D<double> httAinvH;
      httAinvH=subtract(ht,ttAinvH);

      for(int jt=it;jt<totgrid;jt++){
        Array1D<double> xcurra(ndim_,0.e0);
        for (int id=0;id<ndim_;id++)
          xcurra(id)=xgrid(jt,id);
      
        Array1D<double> tta(npts,0.e0);
        for(int ipts=0;ipts<npts;ipts++){
          Array1D<double> xdata_i(ndim_,0.e0);
          for (int id=0;id<ndim_;id++)
            xdata_i(id)=xdata_(ipts,id);
          
          tta(ipts)=covariance(xcurra,xdata_i,param_);
        }

        Array1D<double> ttaAinvH;
        prodAlphaMatTVec(AinvH_,tta,1.0,ttaAinvH);

        
        Array1D<double> hta(npc_,0.e0);
        for(int ipc=0;ipc<this->npc_;ipc++)
          hta(ipc)=H_atgrid(jt,ipc);

        
        Array1D<double> ht_t_AinvH;
        ht_t_AinvH=subtract(hta,ttaAinvH);

        double correction2=prod_vecTmatvec(httAinvH,Vst_,ht_t_AinvH);
        double correction1=prod_vecTmatvec(tt,Ainv_, tta);
        
        cov_(it,jt)=(sig2hat_*(covariance(xcurr,xcurra,param_)-correction1+correction2));
        if (it!=jt)
          cov_(jt,it)=cov_(it,jt);
    
        if (msc != "msc") break;

      }//jt=it
    }//if m
    var_(it)=cov_(it,it);
  }//it

  return;
}

// Parameters of the resulting student-t prediction
// The actual, sigma-integrated predictions are Student-t
void Gproc::getSttPars(Array1D<double>& sttmat)
{
  sttmat.Clear();
  /// \todo check that full cov_ already defined(i.e. msc) not just diagonal
  int npts=xdata_.XSize();
  int totgrid=cov_.XSize();


    for(int it=0;it<totgrid;it++)
      for(int jt=it;jt<totgrid;jt++)
	sttmat.PushBack(cov_(it,jt)*(npts+2.*al_-npc_-2.)/(npts+2.*al_-npc_)); 
  sttmat.PushBack(npts+2.*al_-npc_-2.);

  return;
}

// Get the full covariance at a given grid
void Gproc::getXYCov(Array2D<double>& xgrid,Array2D<double>& xycov)
{
  /// \todo check that full cov_ already defined(i.e. msc) not just diagonal
  int totgrid=cov_.XSize();

  xycov.Resize(totgrid*totgrid,2*ndim_+1,0.e0);
  int ii=0;
  for(int it=0;it<totgrid;it++){
    for(int jt=0;jt<totgrid;jt++){
      for (int id=0;id<ndim_;id++)
        xycov(ii,id)=xgrid(it,id);
      for (int id=0;id<ndim_;id++)
        xycov(ii,ndim_+id)=xgrid(jt,id);
      if (it<=jt)
        xycov(ii,2*ndim_)=cov_(it,jt);
      else 
        xycov(ii,2*ndim_)=cov_(jt,it);
      ii++;
    }
  }

  return;
}

// Compute covariance value given two data points
double Gproc::covariance(Array1D<double>& x1, Array1D<double>& x2,Array1D<double>& param)
{
  /// \todo put an 'if' check for covtype_
  int n=x1.XSize();
  int nn=x2.XSize();

  if (nn!=n ) {
    printf("Gproc:covariance() : Error message: covariance: matrix dimensions do not match! \n"); exit(1);
  }

  double cov=0.0;
  
  Array2D<double> B(nn,nn,0.e0);
  for(int i=0;i<nn;i++){
    B(i,i)=1./param(i);
  } 
  
  Array1D<double> x12;
  x12=subtract(x1,x2);
  
  Array1D<double> Bx12;
  prodAlphaMatVec (B, x12, 1.0, Bx12) ;
  
  cov=exp(-dot(x12,Bx12));

  return cov;
}

// Compute full covariance matrix given the data
void Gproc::computeDataCov_(Array2D<double>& xdata,Array1D<double>& param,Array2D<double>& A)
{
  int npts=xdata.XSize();
  int ndim=xdata.YSize();
  A.Resize(npts,npts,0.e0);
  for(int i=0;i<npts;i++){
    for(int j=0;j<npts;j++){
        Array1D<double> xdatai(ndim,0.e0);
        Array1D<double> xdataj(ndim,0.e0);
        for(int id=0;id<ndim;id++){
            xdatai(id)=xdata(i,id);
            xdataj(id)=xdata(j,id);
        }

        A(i,j)=this->covariance(xdatai,xdataj,param)+dataVar_(i)*(i==j);
    }
  }

  return;
}

// Find the best correlation length parameter according to max-posterior
void Gproc::findBestCorrParam()
{
  Array1D<double> logparam(ndim_,0.0);
  int n=ndim_;
  int m=5;
  Array1D<int> nbd(n,0);
  Array1D<double> l(n,0.e0);
  Array1D<double> u(n,0.e0);
  
  void* info=this;

  lbfgsDR(n,m,logparam.GetArrayPointer(),nbd.GetArrayPointer(),l.GetArrayPointer(),u.GetArrayPointer(),neglogPostParam,NULL,info) ;
  
  
  for (int i=0;i<ndim_;i++)
      param_(i)=exp(logparam(i));
  
  return;
}

// Function to compute negative log posterior (needed to maximize with respect to roughness parameter)
double neglogPostParam(int ndim, double* mm, void* classpointer)
{
    
    Gproc* thisClass=(Gproc*) classpointer;

    assert(ndim==thisClass->getNdim());
    Array1D<double> param(ndim,0.e0);
    
    for(int i=0;i<ndim;i++)
        param(i)=exp(mm[i]);

       
    thisClass->setCorrParam(param);
    thisClass->BuildGP();
    
    double lpp=-(thisClass->getNpt()+2.*thisClass->getAl()-thisClass->getNPC())/2.;
    
    lpp*=log(thisClass->getSig2hat());
    
    
    Array2D<double> vst;
    thisClass->getVst(vst);
    
    Array2D<double> acor;
    thisClass->getA(acor);

    lpp+=0.5*logdeterm(vst);
    lpp-=0.5*logdeterm(acor);
    
    return -lpp;
}



