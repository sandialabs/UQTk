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
/// \file rosenblatt.cpp
/// \brief Tools related to Rosenblatt transformation
#include <math.h>
#include <iostream>
#include <float.h>
#include "Array1D.h"
#include "Array2D.h"
#include "tools.h"
#include "arraytools.h"

using namespace std;

// Implementation of inverse Rosenblatt map given dimension-specific bandwidths
void invRos(Array1D<double>& unif, Array2D<double>& xi, Array1D<double>& newXi, Array1D<double>& sig)
{
  // Accuracy, i.e. the stopping criterion in the bisection algorithm
  double xacc=1e-8;
  // The initial range of the bisection algorithm [-xmax,xmax]
  double xmax = 1.0e+8;


  // Dimensionality of the transformation
  int ndim = unif.XSize();
  int ns   = xi.YSize();

  // Output container
  newXi.Resize(ndim,0.e0);

  // Dimension check
  if (ndim != (int) xi.XSize()  || ndim != (int) sig.XSize())
    {printf("invRos: dimension error\n"); exit(1);}


  // Work arrays
  Array1D<double> kern(ns,1.e0);
  Array1D<double> numer(ndim,0.e0);
  Array1D<double> denom(ndim,0.e0);

  // Loop over dimensions
  for(int id=0; id < ndim; id++){

    ////  including this is equivalent to dropping | conditionals in Rosenblatt transformation
    //  kern.Resize(ndim,1.e0);

    // Get the corresponding uniform r.v. sample
    double Un=unif(id);

    // Starting point of bisection
    double xx = 0.0;

    // Cure against extremes
    if (Un==0.0)
      Un += DBL_EPSILON;
    if (Un==1.0)
      Un -= DBL_EPSILON;

    // Build denominator
    for(int is=0; is < ns; is++)
      denom(id) += kern(is);

    // Bisection iteration counter
    int ii=0;

    // Bisection looks in the interval [x1,x2]
    double x1 = -xmax;
    double x2 =  xmax;

    // Bisection loop
    do{

      numer(id) = 0 ;
      for(int is = 0; is < ns; is++)
        numer(id) += (kern(is)/denom(id)) *  ( 0.5+0.5*erf( (xx-xi(id,is))/(sqrt(2)*sig(id)) ) ) ;// (xx > xi(id,is));

      if (Un < numer(id)){ x1 = x1; x2 = xx; }
      else               { x1 = xx; x2 = x2; }
      xx = 0.5*(x1+x2);
      ii++;

    } while (fabs(x2-x1) > xacc);
    // End of bisection loop

    // Just-in-case warning
    if (ii>1000)
      printf( "Warning: Inverted CDF in really many (%d) iterations\n", ii);

    // Select the solution as a new sample
    newXi(id) = xx;

    // Update the kernel for the next dimension
    double sig2 = sig(id)*sig(id);
    for(int is=0; is < ns; is++)
      kern(is) = kern(is)*exp(-pow(newXi(id)-xi(id,is),2)/(2*sig2));  // (bw*sqrt(2*PI)) cancels;

  } // end of loop over dimensions

  return;

}

// Implementation of inverse Rosenblatt map given same bandwidth for all dimensions
void invRos(Array1D<double>& unif,  Array2D<double>& xi, Array1D<double>& newXi, double bw)
{
  // Sanity check
  if (bw<=0)
    {printf("invRos: bandwidth needs to be positive"); exit(1);}

  // Get dimensions
  int ndim = unif.XSize();

  // Dimension check
  if (ndim != (int) xi.XSize())
    {printf("invRos: dimension error"); exit(1);}

  // Populate bandwidth vector
  Array1D<double> sig(ndim,bw);

  // Perform inverse Rosenblatt
  invRos(unif,xi,newXi,sig);

  return;
}

// Implementation of inverse Rosenblatt map with bandwidths chosen according to a rule-of-thumb
void invRos(Array1D<double>& unif, Array2D<double>& xi, Array1D<double>& newXi)
{
  // Get dimensions
  int ndim = unif.XSize();
  const int nr_proj = xi.YSize();

  // Dimension check
  if (ndim != (int) xi.XSize())
    {printf("invRos: dimension error"); exit(1);}

  // Transpose the input for bandwidth selection code
  Array2D<double> xi_t(nr_proj,ndim,0.e0);
  transpose(xi,xi_t);

  // Get optimal bandwidths for all dimensions
  Array1D<double> sig;
  get_opt_KDEbdwth(xi_t,sig);

  // Perform inverse Rosenblatt
  invRos(unif,xi,newXi,sig);

  return;
}

// Implementation of inverse Rosenblatt map with bandwidths chosen according to a rule-of-thumb
// operating on a set of uniform samples.
void invRos(Array2D<double>& unif, Array2D<double>& xi, Array2D<double>& newXi)
{
  int npts = unif.XSize();
  int ndim = unif.YSize();
  int nspl = xi.YSize();

  // dimension check
  if (ndim != (int) xi.XSize())
    {printf("invRos: dimension error\n"); exit(1);}

  // Transpose the input for bandwidth selection code
  Array2D<double> xi_t(nspl,ndim,0.e0);
  transpose(xi,xi_t);

  // Get optimal bandwidths for all dimensions
  Array1D<double> sig;
  get_opt_KDEbdwth(xi_t,sig);

  // Inverse rosenblatt for each point in unif
  Array1D<double> uin(ndim),xiout(ndim);
  newXi.Resize(npts,ndim);
  for (int is=0; is<npts; is++) {
    for (int j=0;j<ndim;j++) uin(j) = unif(is,j);
    invRos(uin,xi,xiout,sig);
    for (int j=0;j<ndim;j++) newXi(is,j) = xiout(j);
  }

  return;
}

// A rule-of-thumb for optimal bandwidth selection
void get_opt_KDEbdwth(const Array2D<double>& data,Array1D<double>& bdwth)
{
  // Get dimensions
  int ndata=data.XSize();
  int ndim=data.YSize();

  // Bandwidth vector container
  bdwth.Resize(ndim);

  // Flag that measures proximity to boundary per dimension
  Array1D<double> flag(ndim,1.);

  // Initialize minimum and maximum per dimension
  Array1D<double> datamin(ndim,1000.0);
  Array1D<double> datamax(ndim,-1000.0);

  // Loop over dimensions
  for(int idim=0;idim<ndim;idim++){

    // Fnding min and max values
    for(int idata=0;idata<ndata;idata++){
      if (data(idata,idim)<datamin(idim)) {datamin(idim)=data(idata,idim);}
      if (data(idata,idim)>datamax(idim)) {datamax(idim)=data(idata,idim);}
    }
    // Define proximity to boundary
    double nearBorder=(datamax(idim)-datamin(idim))/20.;

    // Compute the number of samples near boundaries
    int numBorder=0;
    for(int idata=0;idata<ndata;idata++){
       if ( data(idata,idim)-datamin(idim)<nearBorder || datamax(idim)-data(idata,idim)<nearBorder ) {numBorder++;}
    }
    // Set the flag if too many samples near boundary
    if ( numBorder > ndata/20.) {flag(idim)=0.5;}
  } // end loop over dimensions

  // Standard deviation container
  Array1D<double> stdd(ndim,0.e0);

  // Auxiliary 1d vector for st-deviation computation
  Array1D<double> data_1d(ndata);

  // Loop over dimensions
  for(int idim=0;idim<ndim;idim++){
    // Get the 1-d slice
    for(int idata=0;idata<ndata;idata++)
      data_1d(idata)=data(idata,idim);

    // Compute standard deviation for each dimension
    stdd(idim)=get_std(data_1d)+1.e-16;
    // Compute the rule-of-thumb bandwidth with a factor (flag) that accounts for proximity to boundary
    bdwth(idim)=flag(idim)*pow(4./(ndim+2),1./(ndim+4.))*stdd(idim)*pow(ndata,-1./(ndim+4.));
  } // end loop over dimensions


  return;
}

///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////

// Implementation of Rosenblatt map given dimension-specific bandwidths
void Rosen(Array2D<double>& xi, Array2D<double>& xi_data, Array2D<double>& unif, Array1D<double>& sig)
{
  // Get the dimensions
  int nsd=xi_data.XSize();
  int ns=xi.XSize();
  int nd=xi.YSize();

  // Dimension check
  if (nd != (int) xi_data.YSize())
    {printf("Rosen: dimension error"); exit(1);}

  // Output samples container
  unif.Resize(ns,nd);

  // Loop over samples
  for (int js = 0; js < ns; js++) {

    // Setting the 'kernel' for summation
    Array1D<double> kern(nsd,1.0);

    // Loop over dimensions
    for (int id=0;id<nd;id++){

      // Initialize numerator and denominator
      double numer=0.0, denom=0.0;

      // Computing the numerator and denominator
      for (int is=0;is<nsd;is++){
        denom += kern(is);
        numer += (kern(is)*(0.5+0.5*erf((xi(js,id)-xi_data(is,id))/(sqrt(2.)*sig(id)))));
      }

      // Update the kernel
      double sig2 = sig(id)*sig(id);
      for (int is = 0; is < nsd; is++){
        kern(is) = kern(is)*exp(-pow(xi(js,id)-xi_data(is,id),2))/(2.*sig2);
      }

      // Compute the sample value
      unif(js,id) = numer / denom;

    } // end of loop over dimensions

  } // end of loop over samples

  return ;

}

// Implementation of Rosenblatt map given same bandwidth for all dimensions
void Rosen(Array2D<double>& xi, Array2D<double>& xi_data,Array2D<double>& unif, double bw)
{
  // Sanity check
  if (bw<=0)
    {printf("Rosen: bandwidth needs to be positive"); exit(1);}

  // Get the dimension
  int ndim = xi.YSize();

  // Populate dimension-specific bandwidth vector
  Array1D<double> sig(ndim,bw);

  // Perform Rosenblatt map
  Rosen(xi,xi_data,unif,sig);

  return;
}

// Implementation of Rosenblatt map given with optimal dimension-specific bandwidth selection
void Rosen(Array2D<double>& xi,Array2D<double>& xi_data, Array2D<double>& unif)
{
  // Get the dimension
  int ndim = xi.YSize();

  // Compute the optimal bandwidths
  Array1D<double> sig;
  get_opt_KDEbdwth(xi_data,sig);

  // Perform Rosenblatt map
  Rosen(xi,xi_data,unif,sig);

  return;

}
