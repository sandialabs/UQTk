/* =====================================================================================

                      The UQ Toolkit (UQTk) version 3.1.1
                          Copyright (2021) NTESS
                        https://www.sandia.gov/UQToolkit/
                        https://github.com/sandialabs/UQTk

     Copyright 2021 National Technology & Engineering Solutions of Sandia, LLC (NTESS).
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
/// \file PCBasis.cpp
/// \author B. Debusschere, C. Safta, K. Sargsyan, K. Chowdhary 2007 -
/// \brief Univariate PC class

#include "PCBasis.h"
#include "error_handlers.h"
#include "uqtkconfig.h"
#include <math.h>

#include "quad.h"
#include "arrayio.h"
#include "pcmaps.h"
#include "combin.h"

#include <iostream>
#include <string.h>
#include <stdio.h>
#include <sstream>
using namespace std; // needed for python string conversion

PCBasis::PCBasis( const string type, const double alpha, const double betta, const int maxord):
  type_(type), maxord_(maxord), alpha_(alpha), beta_(betta)
{

  // Make sure a valid basis type is requested
  if(!strcmp(type_.c_str(),"HG") || !strcmp(type_.c_str(),"WH") || !strcmp(type_.c_str(),"HER") || !strcmp(type_.c_str(),"HERMITE")){
    narg_=0;
    type_ = "HG"; // Hermite-Gaussian
  } else if (!strcmp(type_.c_str(),"LU") || !strcmp(type_.c_str(),"LEG") || !strcmp(type_.c_str(),"LEGENDRE")){
    narg_ = 0;
    type_ = "LU"; // Legendre-Uniform
  } else if (!strcmp(type_.c_str(),"LU_N") || !strcmp(type_.c_str(),"LEG_N") || !strcmp(type_.c_str(),"LEGENDRE_N")){
    narg_ = 0;
    type_ = "LU_N"; // Legendre-Uniform-Normalized
  } else if (!strcmp(type_.c_str(),"GLG") || !strcmp(type_.c_str(),"LG") || !strcmp(type_.c_str(),"LAG") || !strcmp(type_.c_str(),"LAGUERRE")){
    narg_ = 1;
    type_ = "LG"; // Laguerre-Gamma
  } else if (!strcmp(type_.c_str(),"JB") || !strcmp(type_.c_str(),"BJ") || !strcmp(type_.c_str(),"JAC") || !strcmp(type_.c_str(),"JACOBI")){
    narg_ = 2;
    type_ = "JB"; // Jacobi-Beta
  } else if (!strcmp(type_.c_str(),"SW") || !strcmp(type_.c_str(),"GW") || !strcmp(type_.c_str(),"WIG") || !strcmp(type_.c_str(),"WIGERT")){
    narg_ = 2;
    type_ = "SW"; // Gauss-Wigert, or Stieltjes-Wigert (lognormal)
  }
  else if (!strcmp(type_.c_str(),"pdf")){
    narg_ = 0;
    type_ = "pdf"; // Custom
  }
  else {
    string err_message = (string) "PCBasis::PCBasis(): requested basis type, " + type_ +
                         (string) ", is not supported";
    throw Tantrum(err_message);
  }

  // The default implementation relies on N_q=p+1 quadrature points, where p is the maximal order and N_q is the number of quadrature points
  this->Init1dQuadPoints(maxord+1);
  this->Eval1dBasisAtQuadPoints();
  this->Eval1dNormSq(maxord);

  // Seed for random number generator in case we need to sample
  // the basis functions
  int startSeed = 1;
  this->SeedRandNumGen(startSeed);

  return;
}

void PCBasis::Init1dQuadPoints(int nqdpts)
{
  // Move strings to char*, since Quad class needs char*
  char type_c[10];
  strcpy(type_c, type_.c_str());
  Quad spRule(type_c,(char *)"full",1,nqdpts,alpha_,beta_); // CHECK!
  spRule.SetRule();

  // Get the integration rule information
  spRule.GetRule(quadPoints_, quadWeights_,quadIndices_);

  return;
}

void PCBasis::Eval1dBasisAtQuadPoints()
{
  Array1D<double> quadPoints1d;
  for(int iq=0;iq<(int)quadPoints_.XSize();iq++)
    quadPoints1d.PushBack(quadPoints_(iq,0));

  this->Eval1dBasisAtCustPoints(psi1d_,maxord_,quadPoints1d);


  return;
}

void PCBasis::Eval1dBasisAtCustPoints(Array2D<double>& psi,int kord, const Array1D<double>& custPoints)
{

  int npts=custPoints.XSize();

  psi.Resize(npts,kord+1,0.e0);

  // Evaluate the basis functions at each of the quadrature points
  for(int isp=0; isp < npts; isp++){
    // Define a temporary array to store the basis
    // values in for this quadrature point per dimensions
    Array1D<double> basisVals(kord+1);

    this->EvalBasis(custPoints(isp),basisVals);

    // Store the results back in the psi_ array
    for(int iord=0; iord < kord+1; iord++)
      psi(isp,iord) =  basisVals(iord);

  }

  return;

}


double PCBasis::EvalBasis(const double &xi, Array1D<double> &basisEvals) const
{

  // Get the order, up to which the Basis is computed
  int kord=basisEvals.Length()-1;

  return (this->EvalBasis(xi, kord, basisEvals.GetArrayPointer()));

}

double PCBasis::EvalBasis(const double &xi, const int kord, double *basisEvals) const
{

  // Use Recursion formulas to evaluate polynomials at the requested xi value
  if(!strcmp(type_.c_str(),"HG")){ // Hermite-Gaussian

    basisEvals[0] = 1.e0;

    if(kord > 0){
      basisEvals[1] = xi;
      for(int iord=2; iord < kord+1; iord++){
        basisEvals[iord] = xi*basisEvals[iord-1]-(double)(iord-1)*basisEvals[iord-2];
      }
    } /* done if kord>0 for HG */

 } else if (!strcmp(type_.c_str(),"LU")){ // Legendre-Uniform

    double ca,cb;
    basisEvals[0] = 1.e0;
    if(kord > 0){
      basisEvals[1] = xi;
      for(int iord=2; iord < kord+1; iord++){
        ca = (2.e0*(double)(iord-1) + 1.e0)*xi;
        cb = (double)(iord-1);
        basisEvals[iord] = (ca*basisEvals[iord-1] - cb*basisEvals[iord-2])/(double)(iord);
      }
    } /* done if kord>0 for LU */

  } else if (!strcmp(type_.c_str(),"LU_N")){ // Legendre-Uniform-Normailzed

    double ca,cb;
    basisEvals[0] = 1.e0;
    if(kord > 0){
      basisEvals[1] = xi*sqrt(3.0);
      for(int iord=2; iord < kord+1; iord++){
        ca = sqrt(2.e0*(double)(iord-1) + 1.e0)*xi;
        cb = (double)(iord-1)/sqrt(2.e0*(double)(iord-2) + 1.e0);
        basisEvals[iord] = sqrt(2.e0*(double)(iord) + 1.e0)*(ca*basisEvals[iord-1] - cb*basisEvals[iord-2])/(double)(iord);
      }
    }  /* done if kord>0 for LU_N */

  } else if (!strcmp(type_.c_str(),"LG")){ // Laguerre-Gamma

    double ca,cb;
    basisEvals[0] = 1.e0;
    if(kord > 0){
      basisEvals[1] = -xi + alpha_ + 1.e0;
      for(int iord=2; iord < kord+1; iord++){
        ca = 2.e0*(double)(iord-1) + alpha_ + 1.e0 - xi;
        cb = (double)(iord-1) + alpha_;
        basisEvals[iord] = (ca*basisEvals[iord-1] - cb*basisEvals[iord-2])/(double)(iord);
      }
    }  /* done if kord>0 for LG*/

  } else if (!strcmp(type_.c_str(),"JB")){ // Jacobi-Beta

    double ca,cb,cc;
    basisEvals[0] = 1.e0;
    if(kord > 0){
      basisEvals[1] = (alpha_-beta_)/2.e0+xi*(alpha_+beta_+2.e0)/2.e0;
      for(int iord=2; iord < kord+1; iord++){
        ca = (2.e0*(double)(iord) + alpha_ + beta_ - 1.e0)*(alpha_*alpha_-beta_*beta_) + xi*(2.e0*(double)(iord)+alpha_+beta_-2.e0)*(2.e0*(double)(iord)+alpha_+beta_-1.e0)*(2.e0*(double)(iord)+alpha_+beta_);
        cb = 2.e0*((double)(iord-1) + alpha_)*((double)(iord-1) + beta_)*(2.e0*(double)(iord) + alpha_+beta_);
        cc = 2.e0*(double)(iord)*(double(iord)+alpha_+beta_)*(2.e0*(double)(iord-1)+alpha_+beta_);
        basisEvals[iord] = (ca*basisEvals[iord-1] - cb*basisEvals[iord-2])/cc;
      }
    }  /* done if kord>0 for JB */

  } else if (!strcmp(type_.c_str(),"SW")){ // Wigert-Lognormal

    double ca,cb,cc;
    double ee   = exp(beta_*beta_/2.);
    double eesq = ee*ee ;
    basisEvals[0] = 1.e0;
    if(kord > 0){
      basisEvals[1] = xi-ee*exp(alpha_);
      for(int iord=2; iord < kord+1; iord++){
        ca = xi-exp(alpha_)*pow(ee,2.e0*double(iord)-3.e0)*((eesq+1.0)*pow(eesq,(double)(iord)-1.e0)-1.e0);
        cb = exp(2.e0*alpha_)*pow(eesq,3.e0*(double)(iord)-5.e0)*(pow(ee,2.e0*(double)(iord-1))-1.e0);
        cc = 1.e0;
        basisEvals[iord] = (ca*basisEvals[iord-1] - cb*basisEvals[iord-2])/cc;
      }
    }  /* done if kord>0 for SW */

  } else if (!strcmp(type_.c_str(),"pdf")){ // Custom
    Array2D<double> albe;
    read_datafileVS(albe,"ab.dat");
    if ((int)albe.XSize()<maxord_+1)
    {
      std::cout<<"order+1 = "<<maxord_+1<<", nrows(ab.dat) = "<<albe.XSize()<<std::endl<<std::flush ;
      throw Tantrum("PCBasis:: EvalBasis: The input coefficient file 'ab.dat' has fewer rows than needed");
    }
    double ca,cb,cc;
    basisEvals[0] = 1.e0;
    if(kord > 0){
      basisEvals[1] = xi-albe(0,0);
      for(int iord=2; iord < kord+1; iord++){
        ca = xi-albe(iord-1,0);
        cb = albe(iord-1,1);
        cc = 1.e0;
        basisEvals[iord] = (ca*basisEvals[iord-1] - cb*basisEvals[iord-2])/cc;
      }
    }   /* done if kord>0 for pdf */
  }  else {
    string err_message = (string) "PCBasis:: EvalBasis: invalid basis type";
    throw Tantrum(err_message);
  }

  // Return the basis value for the highest order
  return basisEvals[kord];

}

/***********************************************************
Derivative methods
*************************************************************/

/*****************************************************
1st Derivatives of 1d Legendre Polynomials
******************************************************/
void PCBasis::EvalDerivBasis(const double& xi, Array1D<double>& basisDEvals){

  if (!strcmp(type_.c_str(),"LU")){
    // Only works for non-normalized Legendre uniform PCE

    int kord=basisDEvals.Length()-1;

    // get basis evals at xi (not derivative)
    Array1D<double> basisEvals(kord+1,0.0);
    EvalBasis(xi, basisEvals);

    basisDEvals(0) = 0.e0; // derivative of zeroth order leg poly

    if(kord > 0){

      basisDEvals(1) = 1; // derivative of first legendre poly.

      for(int iord=2; iord < kord+1; iord++){
        basisDEvals(iord) = basisDEvals(iord-2) + (2*iord - 1)*basisEvals(iord-1);
      }
    }
  }
  else{
    string err_message = (string) "PCBasis:: EvalBasis: invalid basis type";
    throw Tantrum(err_message);
  }

  return;

}

void PCBasis::Eval1dDerivBasisAtCustPoints(Array2D<double>& dpsi,int kord, const Array1D<double>& custPoints){

  if (!strcmp(type_.c_str(),"LU")){

    int npts=custPoints.XSize();
    dpsi.Resize(npts,kord+1,0.e0);

    // Evaluate the derivative of the basis functions at each of the quadrature points
    for(int isp=0; isp < npts; isp++){
      // Define a temporary array to store the basis
      // values in for this quadrature point per dimensions
      Array1D<double> basisDVals(kord+1);

      EvalDerivBasis(custPoints(isp),basisDVals);

      // Store the results back in the psi_ array
      for(int iord=0; iord < kord+1; iord++)
          dpsi(isp,iord) =  basisDVals(iord);

    }
  }
  else{
    string err_message = (string) "PCBasis:: EvalBasis: invalid basis type";
    throw Tantrum(err_message);
  }
}

/*****************************************************
2nd Derivatives of 1d Legendre Polynomials
******************************************************/
void PCBasis::Eval2ndDerivBasis(const double& xi,Array1D<double>& ddP) {

  if (!strcmp(type_.c_str(),"LU")) {
    int kord=ddP.Length()-1;
    Array1D<double> P(ddP.Length(),0);
    Array1D<double> dP(ddP.Length(),0);

    EvalBasis(xi,P);
    EvalDerivBasis(xi,dP);

    Array1D<double> temp(3,0);
    temp(0) = 0.0; temp(1) = 0.0; temp(2) = 3.0;

    if (kord > 2) {
      ddP(0)=0.0;
      ddP(1)=0.0;
      ddP(2)=3.0;
      for (int n=1; n<kord; n++) {
          ddP(n+1)=((2*n)+1)*dP(n)+ddP(n-1);
      }
    }
    else{
      for (int i = 0; i <= kord; i++){
        ddP(i) = temp(i);
      }
    }
  }
  else{
    string err_message = (string) "PCBasis:: EvalBasis: invalid basis type";
    throw Tantrum(err_message);
  }

  return ;

}

void PCBasis::Eval2ndDerivCustPoints(Array2D<double>& psi, int kord, Array1D<double>& custPoints) {

    if (!strcmp(type_.c_str(),"LU")){
      int npts=custPoints.XSize();

      psi.Resize(npts,kord+1,0.0);

      for (int i=0; i<npts; i++) {
          Array1D<double> ddbasisVals(kord+1);

          Eval2ndDerivBasis(custPoints(i),ddbasisVals);

          for(int iord=0; iord<kord+1; iord++) {
              psi(i,iord) = ddbasisVals(iord);
          }
      }
    }
    else{
      string err_message = (string) "PCBasis:: EvalBasis: invalid basis type";
      throw Tantrum(err_message);
    }
}
//***********************************************************

void PCBasis::GetRandSample(Array1D<double>& randSamples)
{
  int nSamples = randSamples.Length();

  // get the random variable samples
  this->GetRandSample(randSamples.GetArrayPointer(), nSamples);

  return;
}

void PCBasis::GetRandSample(double* randSamples, const int& nSamp)
{
  // Make local copy of nSamp to preserve const requirement on nSamp
  int nSamples = nSamp;


  // Depending on the basis type, pick the right random number generator
  if(!strcmp(type_.c_str(),"HG")){ // Hermite-Gaussian
     for (int i = 0 ; i < nSamples ; i++ )
      randSamples[i] = dsfmt_genrand_nrv(&rnstate_);
  } else if (!strcmp(type_.c_str(),"LU")){ // Legendre-Uniform
 for (int i = 0 ; i < nSamples ; i++ )
   randSamples[i] = dsfmt_genrand_urv_sm(&rnstate_,-1.0,1.0);
  } else if (!strcmp(type_.c_str(),"LU_N")){ // Legendre-Uniform-Normalized
      for (int i = 0 ; i < nSamples ; i++ )
          randSamples[i] = dsfmt_genrand_urv_sm(&rnstate_,-1.0,1.0);
  } else if (!strcmp(type_.c_str(),"LG")){ // Laguerre-Gamma
    for (int i = 0 ; i < nSamples ; i++ )
      randSamples[i] = dsfmt_genrand_urv_sm(&rnstate_,-1.0,1.0);
    // Map to gamma random variable
    for (int i = 0 ; i < nSamples ; i++ ) randSamples[i] = PCtoPC(randSamples[i], "LU", 0, 0, "LG", alpha_, 0) ;

  } else if (!strcmp(type_.c_str(),"JB")){ // Jacobi-Beta

 for (int i = 0 ; i < nSamples ; i++ )
   randSamples[i] = dsfmt_genrand_urv_sm(&rnstate_,-1.0,1.0);
    // Map to beta random variable
    for (int i = 0 ; i < nSamples ; i++ ) randSamples[i] = PCtoPC(randSamples[i], "LU", 0, 0, "JB", alpha_, beta_) ;

  } else if (!strcmp(type_.c_str(),"SW")){ // Stieltjes-Wigert

     for (int i = 0 ; i < nSamples ; i++ )
       randSamples[i] = exp(dsfmt_genrand_nrv_sm(&rnstate_,alpha_,beta_));

  } else if (!strcmp(type_.c_str(),"pdf")){ // Custom-PDF

    string err_message = (string) "PCBasis::GetRandSample(): Custom-PDF sampling is not implemented yet";
    throw Tantrum(err_message);

  }
   else {
    string err_message = (string) "PCBasis::GetRandSample(): requested basis type, " + type_ +
      (string) ", is not supported";
    throw Tantrum(err_message);
  }

  return;
}


void PCBasis::GetQuadRule(Array2D<double>& qPoints, Array1D<double>& qWeights, Array2D<int>& qIndices)
{
this->GetQuadPoints(qPoints);
this->GetQuadWeights(qWeights);
this->GetQuadIndices(qIndices);

return;
}


void PCBasis::SeedRandNumGen(const int& seed)
{

  // Update the stored seed
  this->rSeed_ = seed;

  // seed the appropriate random number generator
  dsfmt_init_gen_rand(&(this->rnstate_), (uint32_t) rSeed_ );

  return;
}

void PCBasis::Eval1dNormSq(int kord)
{
  // Norms^2 of basis functions
  psi1dSq_.Resize(kord+1,0.e0);

  // Take the norm of f^2 by numerical quadrature
  for(int iord=0;iord < kord+1; iord++){
    double tempSum = 0.e0;
    for(int iqp=0; iqp < (int) quadPoints_.XSize(); iqp++){
      tempSum += psi1d_(iqp,iord)*psi1d_(iqp,iord)*quadWeights_(iqp);
    }
    psi1dSq_(iord) = tempSum;
    //cout << "Norm(" << iord << ")=" << tempSum << endl;
  }
  return;
}

void PCBasis::Eval1dNormSq_Exact(int kord)
{
  // Norms^2 of basis functions
  psi1dSqExact_.Resize(kord+1,0.e0);

  // Take the norm of f^2 by numerical quadrature
  for(int iord=0;iord < kord+1; iord++){
    psi1dSqExact_(iord) = NormSq_Exact(iord);
  }
  return;
}


double PCBasis::NormSq_Exact(int kord)
{

  double normSq=1.e0;
  if(!strcmp(type_.c_str(),"HG")){
    for (int i=2;i<=kord;i++) normSq*=i;
  }
  else if(!strcmp(type_.c_str(),"LU"))
    normSq=1.e0/(2.e0*(double)(kord)+1.e0);
  else if(!strcmp(type_.c_str(),"LU_N"))
    normSq=1.e0;
  else if(!strcmp(type_.c_str(),"LG"))
    normSq=exp(lgamma((double)(kord)+alpha_+1.e0)-lgamma(alpha_+1.e0)-lgamma((double)(kord)+1.e0));
  else if(!strcmp(type_.c_str(),"JB")){
   if(kord>0){
     double dk   = (double) kord;
     double ap1  = alpha_+1.e0;
     double bp1  = beta_ +1.e0;
     double abp1 = alpha_+beta_ +1.e0;
     normSq = exp(lgamma(dk+ap1)+lgamma(dk+bp1)-lgamma(dk+1.e0)-lgamma(dk+abp1)-lgamma(ap1)-lgamma(bp1)+lgamma(ap1+bp1))
             /(2.e0*dk+abp1);
   }
   else
     normSq=1.e0;
  }
  else if(!strcmp(type_.c_str(),"SW")){
    if(kord>0){
      //normSq=exp((2.e0*alpha_+1.e0)*kord)*(3.e0*(double)(kord)-1.e0)*beta_*beta_/2.e0;
      normSq=exp((2.e0*alpha_)*kord)*exp((3.e0*(double)(kord)-1.e0)*beta_*beta_*kord/2.e0);
      for(int iord=1;iord<=kord;iord++)
        normSq *= (exp((double)(iord)*beta_*beta_)-1);
    }
    else
      normSq=1.e0;
  }
 else if (!strcmp(type_.c_str(),"pdf")){ // Custom-PDF

   string err_message = (string) "PCBasis::NormSq_Exact(): Custom-PDF norm-squared is not implemented yet";
   throw Tantrum(err_message);

  }
  else {
    string err_message = (string) "PCBasis::NormSq_Exact(): requested basis type, " + type_ +
      (string) ", is not supported";
    throw Tantrum(err_message);
  }

  return normSq;

}

