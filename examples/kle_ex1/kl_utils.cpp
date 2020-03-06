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
#include "kl_utils.h"

double genCovAnl1D(const double x1, const double x2, const double clen, const double sigma, const string covtype) {

  double cov;

  if ( covtype == "SqExp" )
    cov = exp( - (x1-x2) * (x1-x2) / ( clen * clen ) ) * sigma * sigma ;
  else if ( covtype == "Exp" )
    cov = exp( - fabs(x1-x2) / clen ) * sigma * sigma ;
  else
    throw Tantrum("genCovAnl(): covariance type is not recognized!");

  return ( cov );

}

double genCovAnl2D( const Array2D<double> &xy, const int &i, const int &j,
		    const double &clen, const double &sigma, const string &covtype) {

  double cov;

  double dxy = sqrt((xy(i,0)-xy(j,0))*(xy(i,0)-xy(j,0))+(xy(i,1)-xy(j,1))*(xy(i,1)-xy(j,1)));
  if ( covtype == "SqExp" )
    cov = exp( - dxy*dxy / ( clen * clen ) ) * sigma * sigma ;
  else if ( covtype == "Exp" )
    cov = exp( - dxy / clen ) * sigma * sigma ;
  else
    throw Tantrum("genCovAnl(): covariance type is not recognized!");

  return ( cov );

}

double getMean(const Array2D<double> &y, const string &idir, const int &ij) {

  int nx = (int) y.XSize();
  int ny = (int) y.YSize();

  double dmean, dtmp = 0.0 ;

  if ( idir == "L" ) {

    int mst=ny%6;

    for ( int j = 0; j<mst; j++ )
      dtmp += y(ij,j);

    if ( ny >= 6 ) {
      for ( int j = mst; j < ny; j += 6)
        dtmp += (y(ij,j)+y(ij,j+1)+y(ij,j+2)+y(ij,j+3)+y(ij,j+4)+y(ij,j+5));
    }

    dmean = dtmp / ((double) ny );

  }
  else if ( idir == "C" ) {

    int mst=nx%6;

    for ( int i = 0; i<mst; i++ )
      dtmp += y(i,ij);

    if ( nx >= 6 ) {
      for ( int i = mst; i < nx; i += 6)
        dtmp += (y(i,ij)+y(i+1,ij)+y(i+2,ij)+y(i+3,ij)+y(i+4,ij)+y(i+5,ij));
    }

    dmean = dtmp / ( (double) nx );
  }

  return ( dmean ) ;

}

double dcomp(const double xi, const double a, const double b, const double L) {
  double dtmp = pow((b+1.0)/(b-1.0),(xi-a)/(1.0-a));
  return (
           L * ((2.0*a+b)*dtmp+2.0*a-b)
	      /((2.0*a+1.0)*(1.0+dtmp))
         );

}

double scomp(const double xi, const double b, const double L) {
  double dtmp = pow((b+1.0)/(b-1.0),1.0-xi);
  return (
	  L * (b+1.0-(b-1.0)*dtmp)
	     /(dtmp+1.0)
         );

}


void genGrid1D(Array1D<double> &xgrid, const int npts, const double L, const char *type,
             const double a, const double b) {

  if (npts<=0)
    throw Tantrum("kl_sample::genGrid() : number of grid points needs to be greater than 0") ;

  xgrid.Resize(npts,0.0);

  for ( int i=0; i<npts; i++ )
    xgrid(i) = (double) i / ((double) npts - 1.0);

  if ( string(type) == string("unif") )
    for ( int i=0; i<npts; i++ )
      xgrid(i) *= L ;
  else if ( string(type) == string("cl0") )
    for ( int i=0; i<npts; i++ )
      xgrid(i) = scomp(xgrid(i), b, L);
  else if ( string(type) == string("cl0L") )
    for ( int i=0; i<npts; i++ )
      xgrid(i) = dcomp(xgrid(i), a, b, L);
  else {
    cout<<"genGrid1D(): Unknown grid type: "<<string(type)<<endl<<flush;
    std::terminate() ;
  }

  return;

}

void genGrid2D(Array1D<double> &xgrid,Array1D<double> &ygrid, const
             int nx, const int ny, const double L, const char *type,
             const double a, const double b){

  if ( (nx<=0) || (ny<=0) )
    throw Tantrum("kl_sample::genGrid() : number of grid points needs to be greater than 0") ;

  xgrid.Resize(nx,0.0);
  for ( int i=0; i<nx; i++ ) xgrid(i) = (double) i / ((double) nx - 1.0);
  ygrid.Resize(ny,0.0);
  for ( int i=0; i<ny; i++ ) ygrid(i) = (double) i / ((double) ny - 1.0);

  if ( string(type) == string("unif") ) {
    for ( int i=0; i<nx; i++ ) xgrid(i) *= L ;
    for ( int i=0; i<ny; i++ ) ygrid(i) *= L ;
  }
  else if ( string(type) == string("cl0") ) {
    for ( int i=0; i<nx; i++ ) xgrid(i) = scomp(xgrid(i), b, L);
    for ( int i=0; i<ny; i++ ) ygrid(i) = scomp(ygrid(i), b, L);
  }
  else if ( string(type) == string("cl0L") ) {
    for ( int i=0; i<nx; i++ ) xgrid(i) = dcomp(xgrid(i), a, b, L);
    for ( int i=0; i<ny; i++ ) ygrid(i) = dcomp(ygrid(i), a, b, L);
  }
  else {
    cout<<"genGrid2D(): Unknown grid type: "<<string(type)<<endl<<flush;
    std::terminate() ;
  }

  return;

}

void getGrid1dEquiv(const Array1D<double> &xgrid, const Array1D<double> &ygrid, Array1D<double> &xg1d) {

  int nx = (int) xgrid.XSize();
  int ny = (int) ygrid.XSize();
  assert( nx*ny == (int) xg1d.XSize() );

  Array1D<double> w(nx*ny,0.0) ;

  double dxm,dxp,dym,dyp;
  for ( int k = 0; k < nx*ny; k++){
    int i=k%nx;
    int j=k/nx;
    if ( i==0 )
      dxm=0.0;
    else
      dxm=xgrid(i)-xgrid(i-1);
    if ( i==nx-1 )
      dxp=0.0;
    else
      dxp=xgrid(i+1)-xgrid(i);
    if ( j==0 )
      dym=0.0;
    else
      dym=ygrid(j)-ygrid(j-1);
    if ( j==ny-1 )
      dyp=0.0;
    else
      dyp=ygrid(j+1)-ygrid(j);
    w(k)=0.25*(dxm*(dym+dyp)+dxp*(dym+dyp));
  }

  xg1d(1)=w(0)*2.0+xg1d(0);
  for ( int k = 2; k < nx*ny; k++)
    xg1d(k)=2.0*w(k-1)+xg1d(k-2);

  write_datafile_1d(xg1d,"xg1d.dat");

  return ;

}


void getWeights2DU(const Array2D<double> &xgrid, const Array2D<int> &tgrid, Array1D<double> &w) {

  int nxy = xgrid.XSize() ;
  int nt  = tgrid.XSize() ;

  w.Resize(nxy,0.0) ;
  for ( int i = 0; i<nxy; i++)
    for ( int j = 0; j<nt; j++)
      if (( tgrid(j,0) == i ) || ( tgrid(j,1) == i ) || ( tgrid(j,2) == i ))
	w(i) += trArea(xgrid(tgrid(j,0),0),xgrid(tgrid(j,1),0),xgrid(tgrid(j,2),0),
		       xgrid(tgrid(j,0),1),xgrid(tgrid(j,1),1),xgrid(tgrid(j,2),1))/3.0;

  return ;

}

double trArea(const double x1,const double x2,const double x3,
              const double y1,const double y2,const double y3) {

  return (fabs(x1*(y2-y3)+x2*(y3-y1)+x3*(y1-y2))/2.0) ;

}
