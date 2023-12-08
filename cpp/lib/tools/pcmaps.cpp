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
/** \file pcmaps.cpp
 * \brief Suite of functions to help map one kind of a PC variable to another.
 */

#include "Array1D.h"
#include "Array2D.h"
#include "gen_defs.h"
#include "error_handlers.h"
#include "probability.h"
#include "arrayio.h"
#include "pcmaps.h"
#include "combin.h"
#include <assert.h>
#include <math.h>
#include <float.h>


#define EPS 1e-16
#define MXEPS 1e+100

using namespace std;

// Implementation of map from one PC germ to another
double PCtoPC(double x, const std::string pcIn, double in1, double in2, const std::string pcOut, double out1, double out2)
{
  // Output variable
  double y;

  ///////////////////////////////////////
  // Sanity check on PC inputs x
  ///////////////////////////////////////
  if ( pcIn=="HG" || pcIn=="pdf" ){
    // do nothing
  }
  else if (pcIn=="LU" or pcIn=="JB"){
    if (fabs(x)>1.0)
      throw Tantrum("pcmaps.cpp:: PCtoPC:: input x is out of the range!");
  }
  else if (pcIn=="SW" or pcIn=="LG"){
    if (x<0.0)
      throw Tantrum("pcmaps.cpp:: PCtoPC:: input x is out of the range!");
  }
  else if (pcIn=="TG"){
    if (x>in1)
      throw Tantrum("pcmaps.cpp:: PCtoPC:: input x is out of the range!");
  }
  else if (pcIn=="RB"){
    if (x>in2/(1.-in1) or x<0)
      throw Tantrum("pcmaps.cpp:: PCtoPC:: input x is out of the range!");
  }
  else {
    printf("pcIn = %s\n",pcIn.c_str()) ;
    throw Tantrum("pcmaps.cpp::Input PC is not recognized!\n");
  }

  ///////////////////////////////////////
  // Sanity check on input PC parameters
  ///////////////////////////////////////
  if (pcIn=="LU" or pcIn=="HG" or pcIn=="TG" or pcIn=="pdf"){
      // do nothing
  }
  else if (pcIn=="JB"){
    if (in1<=-1. or in2<=-1.)
      throw Tantrum("pcmaps.cpp:: PCtoPC:: input parameter is out of the range!");
  }
  else if (pcIn=="SW"){
     if (in2<=0.)
       throw Tantrum("pcmaps.cpp:: PCtoPC:: input parameter is out of the range!");
  }
  else if (pcIn=="LG"){
    if (in1<-1.)
      throw Tantrum("pcmaps.cpp:: PCtoPC:: input parameter is out of the range!");
  }
  else if (pcIn=="RB"){
    if (in2<=0. or in1<=0)
  	  throw Tantrum("pcmaps.cpp:: PCtoPC:: input parameter is out of the range!");
  }
  else {
    printf("pcIn = %s\n",pcIn.c_str()) ;
    throw Tantrum("pcmaps.cpp::Input PC is not recognized!\n");
  }

  ///////////////////////////////////////
  // Sanity check on output PC parameters
  ///////////////////////////////////////
  if (pcOut=="LU" or pcOut=="HG" or pcOut=="TG" or pcOut=="pdf"){
     // do nothing
  }
  else if (pcOut=="JB"){
    if (out1<=-1. or out2<=-1.)
      throw Tantrum("pcmaps.cpp:: PCtoPC:: output parameter is out of the range!");
  }
  else if (pcOut=="SW"){
    if (out2<=0.)
      throw Tantrum("pcmaps.cpp:: PCtoPC:: output parameter is out of the range!");
  }
  else if (pcOut=="LG") {
    if (out1<0.)
      throw Tantrum("pcmaps.cpp:: PCtoPC:: output parameter is out of the range!");
  }
  else if (pcOut=="RB"){
    if (out2<=0. or out1<=0)
      throw Tantrum("pcmaps.cpp:: PCtoPC:: output parameter is out of the range!");
  }
  else {
    printf("pcOut = %s\n",pcOut.c_str()) ;
    throw Tantrum("pcmaps.cpp::Output PC is not recognized!\n");
  }

  ///////////////////////////////////////
  // Maps when input and output PCs are the same
  ///////////////////////////////////////
  if (pcIn=="LU" && pcOut=="LU"){
    y=x;
  }
  else if (pcIn=="HG" && pcOut=="HG"){
     y=x;
  }
  else if (pcIn=="SW" && pcOut=="SW"){
    y=pow(x,out2/in2)*exp(out1-in1*out2/in2);
  }
  else if (pcIn=="JB" && pcOut=="JB"){
    if(fabs(x)==1) y=x;
    else
      y=PCtoPC(PCtoPC(x,"JB",in1,in2,"LU",0,0),"LU",0,0,"JB",out1,out2);
  }
  else if (pcIn=="LG" && pcOut=="LG"){
    if (x==0) y=0.;
    else y=PCtoPC(PCtoPC(x,"LG",in1,0,"LU",0,0),"LU",0,0,"LG",out1,0);
  }
  else if (pcIn=="TG" && pcOut=="TG"){
    if (x==in1) y=out1;
    else if (in1==out1) y=x;
    else y=PCtoPC(PCtoPC(x,"TG",in1,0,"LU",0,0),"LU",0,0,"TG",out1,0);
  }
  else if (pcIn=="RB" && pcOut=="RB"){
    y=PCtoPC(PCtoPC(x,"RB",in1,in2,"LU",0,0),"LU",0,0,"RB",out1,out2);
  }
  else if (pcIn=="pdf" && pcOut=="pdf"){
    y=x;
  }

  ///////////////////////////////////////
  // Maps when input and output PCs are not the same
  ///////////////////////////////////////

  // Gauss-Hermite //////////////////////
  else if (pcIn=="LU" && pcOut=="HG"){
    if(fabs(x)==1.)
      throw Tantrum("pcmaps.cpp::LU->HG: the value at the domain boundary (would map to infinity)!\n");
    y=sqrt(2.0)*inverf(x);
  }
  else if (pcIn=="HG" && pcOut=="LU"){
    y=erf(x/sqrt(2.0));
  }

  // Gamma-Laguerre //////////////////////
  else if (pcIn=="LG" && pcOut=="LU"){
    y=2.*gammai(in1+1,x)-1.;
  }
  else if (pcIn=="LU" && pcOut=="LG"){
    y=rtbis_mod(PCtoPC,1.e-12,1.e+12,EPS,x,"LG",out1,out2,"LU",in1,in2);
  }

  else if (pcIn=="LG" && pcOut=="HG"){
    y=PCtoPC(PCtoPC(x,"LG",in1,0,"LU",0,0),"LU",0,0,"HG",0,0);
  }
  else if (pcIn=="HG" && pcOut=="LG"){
    y=PCtoPC(PCtoPC(x,"HG",0,0,"LU",0,0),"LU",0,0,"LG",out1,0);
  }

  // Jacobi-Beta //////////////////////
  else if (pcIn=="JB" && pcOut=="LU"){
    y=2.*betai(in1+1,in2+1,(x+1.)/2.)-1.;
  }
  else if (pcIn=="LU" && pcOut=="JB"){
    if(fabs(x)==1) y=x;
    else
      y=rtbis_mod(PCtoPC,-1.,1.,EPS,x,"JB",out1,out2,"LU",in1,in2);
  }

  else if (pcIn=="JB" && pcOut=="HG"){
    y=PCtoPC(PCtoPC(x,"JB",in1,in2,"LU",0,0),"LU",0,0,"HG",0,0);
  }
  else if (pcIn=="HG" && pcOut=="JB"){
    y=PCtoPC(PCtoPC(x,"HG",0,0,"LU",0,0),"LU",0,0,"JB",out1,out2);
  }

  else if (pcIn=="LG" && pcOut=="JB"){
    y=PCtoPC(PCtoPC(x,"LG",in1,0,"LU",0,0),"LU",0,0,"JB",out1,out2);
  }
  else if (pcIn=="JB" && pcOut=="LG"){
    y=PCtoPC(PCtoPC(x,"JB",in1,in2,"LU",0,0),"LU",0,0,"LG",out1,0);
  }

  // Stieltjes-Wigart //////////////////////
  else if (pcIn=="LU" && pcOut=="SW"){
    y=PCtoPC(PCtoPC(x,"LU",0,0,"HG",0,0),"HG",0,0,"SW",out1,out2);
  }
  else if (pcIn=="SW" && pcOut=="LU"){
    y=PCtoPC(PCtoPC(x,"SW",in1,in2,"HG",0,0),"HG",0,0,"LU",0,0);
  }

  else if (pcIn=="HG" && pcOut=="SW"){
    y=exp(out1+out2*x);
  }
  else if (pcIn=="SW" && pcOut=="HG"){
    y=(log(x)-in1)/in2;
  }

  else if (pcIn=="LG" && pcOut=="SW"){
    y=PCtoPC(PCtoPC(x,"LG",in1,0,"LU",0,0),"LU",0,0,"SW",out1,out2);
  }
  else if (pcIn=="SW" && pcOut=="LG"){
    y=PCtoPC(PCtoPC(x,"SW",in1,in2,"LU",0,0),"LU",0,0,"LG",out1,0);
  }

  else if (pcIn=="JB" && pcOut=="SW"){
    y=PCtoPC(PCtoPC(x,"JB",in1,in2,"LU",0,0),"LU",0,0,"SW",out1,out2);
  }
  else if (pcIn=="SW" && pcOut=="JB"){
    y=PCtoPC(PCtoPC(x,"SW",in1,in2,"LU",0,0),"LU",0,0,"JB",out1,out2);
  }

  // Truncated Gaussian //////////////////////
  else if (pcIn=="TG" && pcOut=="LU"){
   y=normcdf(x)/normcdf(in1)*2.-1.;
  }
  else if (pcIn=="LU" && pcOut=="TG"){
    if(fabs(x)==1.)
       throw Tantrum("pcmaps.cpp::LU->TG: the value at the domain boundary (would map to infinity)!\n");
    y=invnormcdf(0.5*(x+1.)*normcdf(out1));
  }

  else if (pcIn=="TG" && pcOut=="HG"){
    y=PCtoPC(PCtoPC(x,"TG",in1,0,"LU",0,0),"LU",0,0,"HG",0,0);
  }
  else if (pcIn=="HG" && pcOut=="TG"){
    y=PCtoPC(PCtoPC(x,"HG",0,0,"LU",0,0),"LU",0,0,"TG",out1,0);
  }

  else if (pcIn=="TG" && pcOut=="LG"){
    y=PCtoPC(PCtoPC(x,"TG",in1,0,"LU",0,0),"LU",0,0,"LG",out1,0);
  }
  else if (pcIn=="LG" && pcOut=="TG"){
    y=PCtoPC(PCtoPC(x,"LG",in1,0,"LU",0,0),"LU",0,0,"TG",out1,0);
  }

  else if (pcIn=="TG" && pcOut=="JB"){
    y=PCtoPC(PCtoPC(x,"TG",in1,0,"LU",0,0),"LU",0,0,"JB",out1,out2);
  }
  else if (pcIn=="JB" && pcOut=="TG"){
    y=PCtoPC(PCtoPC(x,"JB",in1,in2,"LU",0,0),"LU",0,0,"TG",out1,0);
  }

  else if (pcIn=="TG" && pcOut=="SW"){
    y=PCtoPC(PCtoPC(x,"TG",in1,0,"LU",0,0),"LU",0,0,"SW",out1,out2);
  }
  else if (pcIn=="SW" && pcOut=="TG"){
    y=PCtoPC(PCtoPC(x,"SW",in1,in2,"LU",0,0),"LU",0,0,"TG",out1,0);
  }

  // Roe-Baker //////////////////////
  else if (pcIn=="RB" && pcOut=="LU"){
    y=PCtoPC(PCtoPC(x,"RB",in1,in2,"TG",in1,0),"TG",in1,0,"LU",0,0);
  }
  else if (pcIn=="LU" && pcOut=="RB"){
    if(fabs(x)==1.)
      throw Tantrum("pcmaps.cpp::LU->RB: the value at the domain boundary (would map to infinity)!\n");
    y=PCtoPC(PCtoPC(x,"LU",0,0,"TG",out1,0),"TG",out1,0,"RB",out1,out2);
  }

  else if (pcIn=="RB" && pcOut=="HG"){
     y=PCtoPC(PCtoPC(x,"RB",in1,in2,"LU",0,0),"LU",0,0,"HG",0,0);
  }
  else if (pcIn=="HG" && pcOut=="RB"){
     y=PCtoPC(PCtoPC(x,"HG",0,0,"LU",0,0),"LU",0,0,"RB",out1,out2);
  }

  else if (pcIn=="RB" && pcOut=="LG"){
    y=PCtoPC(PCtoPC(x,"RB",in1,in2,"LU",0,0),"LU",0,0,"LG",out1,0);
  }
  else if (pcIn=="LG" && pcOut=="RB"){
    y=PCtoPC(PCtoPC(x,"LG",in1,0,"LU",0,0),"LU",0,0,"RB",out1,out2);
  }

  else if (pcIn=="RB" && pcOut=="JB"){
    y=PCtoPC(PCtoPC(x,"RB",in1,in2,"LU",0,0),"LU",0,0,"JB",out1,out2);
  }
  else if (pcIn=="JB" && pcOut=="RB"){
    y=PCtoPC(PCtoPC(x,"JB",in1,in2,"LU",0,0),"LU",0,0,"RB",out1,out2);
  }

  else if (pcIn=="RB" && pcOut=="SW"){
    y=PCtoPC(PCtoPC(x,"RB",in1,in2,"LU",0,0),"LU",0,0,"SW",out1,out2);
  }
  else if (pcIn=="SW" && pcOut=="RB"){
    y=PCtoPC(PCtoPC(x,"SW",in1,in2,"LU",0,0),"LU",0,0,"RB",out1,out2);
  }

  else if (pcIn=="TG" && pcOut=="RB"){
    y=out2/(1.-PCtoPC(x,"TG",in1,0,"TG",out1,0));
  }
  else if (pcIn=="RB" && pcOut=="TG"){
    y=PCtoPC(1.-in2/x,"TG",in1,0,"TG",out1,0);
  }

  // Custom PDF given by cdf.dat /////////////////
  else if (pcIn=="LU" && pcOut=="pdf"){
    Array2D<double> cdf ;
    read_datafileVS(cdf,"cdf.dat");
    linint( cdf, 0.5*(x+1) , y, 1 ) ;
  }
  else if (pcIn=="pdf" && pcOut=="LU"){
    Array2D<double> cdf ;
    read_datafileVS(cdf,"cdf.dat");
    linint( cdf, x , y, 0 ) ;
    y = y*2.0-1.0 ;
  }

  else if (pcIn=="pdf" && pcOut=="HG"){
    y=PCtoPC(PCtoPC(x,"pdf",0,0,"LU",0,0),"LU",0,0,"HG",0,0);
  }
  else if (pcIn=="HG" && pcOut=="pdf"){
    y=PCtoPC(PCtoPC(x,"HG",0,0,"LU",0,0),"LU",0,0,"pdf",0,0);
  }

  else if (pcIn=="pdf" && pcOut=="LG"){
    y=PCtoPC(PCtoPC(x,"pdf",0,0,"LU",0,0),"LU",0,0,"LG",out1,0);
  }
  else if (pcIn=="LG" && pcOut=="pdf"){
    y=PCtoPC(PCtoPC(x,"LG",in1,0,"LU",0,0),"LU",0,0,"pdf",0,0);
  }

  else if (pcIn=="pdf" && pcOut=="JB"){
    y=PCtoPC(PCtoPC(x,"pdf",0,0,"LU",0,0),"LU",0,0,"JB",out1,out2);
  }
  else if (pcIn=="JB" && pcOut=="pdf"){
    y=PCtoPC(PCtoPC(x,"JB",in1,in2,"LU",0,0),"LU",0,0,"pdf",0,0);
  }

  else if (pcIn=="pdf" && pcOut=="SW"){
    y=PCtoPC(PCtoPC(x,"pdf",0,0,"LU",0,0),"LU",0,0,"SW",out1,out2);
  }
  else if (pcIn=="SW" && pcOut=="pdf"){
    y=PCtoPC(PCtoPC(x,"SW",in1,in2,"LU",0,0),"LU",0,0,"pdf",0,0);
  }

  else if (pcIn=="pdf" && pcOut=="TG"){
    y=PCtoPC(PCtoPC(x,"pdf",0,0,"LU",0,0),"LU",0,0,"TG",out1,0);
  }
  else if (pcIn=="TG" && pcOut=="pdf"){
    y=PCtoPC(PCtoPC(x,"TG",in1,0,"LU",0,0),"LU",0,0,"pdf",0,0);
  }

  else if (pcIn=="pdf" && pcOut=="RB"){
    y=PCtoPC(PCtoPC(x,"pdf",0,0,"LU",0,0),"LU",0,0,"RB",out1,out2);
  }
  else if (pcIn=="RB" && pcOut=="pdf"){
    y=PCtoPC(PCtoPC(x,"RB",in1,in2,"LU",0,0),"LU",0,0,"pdf",0,0);
  }

  else
    throw Tantrum("pcmaps.cpp::Input-Output PC pair is not recognized!\n");

   return y;
}

// Entrywise map from a 2d-array to a 2-d array
void PCtoPC(Array2D<double>& xx, const std::string pcIn, double in1, double in2, Array2D<double>& yy, const std::string pcOut, double out1, double out2)
{

  // Set the dimension
  int n=xx.XSize();
  int m=xx.YSize();

  // Output array
  yy.Resize(n,m);

  // Entrywise map
  for(int i=0;i<n;i++)
    for(int j=0;j<m;j++)
      yy(i,j)=PCtoPC(xx(i,j),pcIn,in1,in2,pcOut,out1,out2);

  return;
}

// Bisection method modified for PC maps to help compute inverse maps
double rtbis_mod(double func(double,const std::string,double,double,const std::string,double,double), const double x1, const double x2, const double xacc,double x, const std::string pcIn, double in1, double in2, const std::string pcOut, double out1, double out2)
{
	const int JMAX=400;
	int j;
	double dx,f,fmid,xmid,rtb;

	f=func(x1,pcIn,in1,in2,pcOut,out1,out2)-x;
	fmid=func(x2,pcIn,in1,in2,pcOut,out1,out2)-x;


	if (f*fmid > 0.0) {printf("Root must be bracketed for bisection in rtbis"); exit(1);}// made the inequality strict
	rtb = f < 0.0 ? (dx=x2-x1,x1) : (dx=x1-x2,x2);
	for (j=0;j<JMAX;j++) {
                xmid=rtb+(dx *= 0.5);
		fmid=func(xmid,pcIn,in1,in2,pcOut,out1,out2)-x;
		if (fmid <= 0.0) rtb=xmid;
		if (fabs(dx) < xacc || fmid == 0.0) return rtb;
	}
	printf("Too many bisections in rtbis"); exit(1);
	return 0.0;
}


// Linear interpolation according to first or second column of a given 2d array
void linint( Array2D<double> &xydata, const double x, double &y, int col )
{

  assert(col==0 || col==1) ;

  int n = xydata.XSize() ;

  if ( x < xydata(0,col) )
    y = 0.0 ;
  else if ( x > xydata(n-1,col) )
    y = 0.0 ;

  int i = 0 ;
  while ( ( x > xydata(i+1,col) ) && ( i < n-2 ) ) i++ ;

  y = xydata(i,1-col) + ( x - xydata(i,col) )
     *( xydata(i+1,1-col) - xydata(i,1-col) ) / ( xydata(i+1,col) - xydata(i,col) ) ;

  return ;

}

// Linear interpolation according to first column of a given 2d array
void linint( Array2D<double> &xydata, const double x, double &y )
{

  int n = xydata.XSize() ;

  if ( x < xydata(0,0) )
    y = 0.0 ;
  else if ( x > xydata(n-1,0) )
    y = 0.0 ;

  int i = 0 ;
  while ( ( x > xydata(i+1,0) ) && ( i < n-2 ) ) i++ ;

  y = xydata(i,1) + ( x - xydata(i,0) )
     *( xydata(i+1,1) - xydata(i,1) ) / ( xydata(i+1,0) - xydata(i,0) ) ;

  return ;

}
