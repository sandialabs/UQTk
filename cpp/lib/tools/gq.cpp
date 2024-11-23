/* =====================================================================================

                      The UQ Toolkit (UQTk) version 3.1.5
                          Copyright (2024) NTESS
                        https://www.sandia.gov/UQToolkit/
                        https://github.com/sandialabs/UQTk

     Copyright 2024 National Technology & Engineering Solutions of Sandia, LLC (NTESS).
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
/// \file gq.cpp
/// \brief Utilities to generate quadrature rules.

#include "stdio.h"
#include "stdlib.h"
#include <iostream>
#include <cmath>

#include "Array1D.h"
#include "Array2D.h"
#include "deplapack.h"
#include "gq.h"
#include "combin.h"
using namespace std;

//#define HERMITE_PROB
//#define HERMITE_PHYS

#define DPI 3.14159265358979323846

double lpol_gq (int n, double x) ;
double hpol_gq (int n, double x) ;
double hpol_phys_gq(int n, double x) ;
double jpol_gq (int n, double a, double b, double x) ;
double jpolp_gq(int n, double a, double b, double x) ;
double lgpol_gq(int n, double a, double x) ;
double fact_gq (int n) ;

/*
  Gauss Quadrature for:
  kind = 1 -> Legendre
         2 -> Chebyshev, 1st kind
         3 -> Chebyshev, 2nd kind
         4 -> Hermite
         5 -> Jacobi
         6 -> Laguerre
  These rules provide the mathematical engine for the quad class.
 */
void gq ( const int kind, const double a, const double b, Array1D<double> &x, Array1D<double> &w ) {

  int n = (int) x.XSize() ;
  gq(kind,n,a,b,x.GetArrayPointer(),w.GetArrayPointer()) ;
  return ;

}


void gq ( const int kind, const int n, const double a, const double b, double *x, double *w ) {

  if ( ( kind < 1 ) || ( kind > 6 ) ) {
    cout<<"ERROR in gq() : kind should be between 1 and 6 : "<<kind<<endl<<flush ;
    terminate() ;
  }

  /* call separate function for Chebyshev quadrature */
  if ( ( kind == 2 ) || ( kind == 3 ) ) {
    gchb(kind-1, n, x, w ) ;
    return ;
  }

  /* Test if n=1 */
  if ( n==1) {
    if ( kind == 1 ) { /* Legendre */
      x[0] = 0.0 ; w[0] = 2.0 ;
    }
    else if ( kind == 4 ) { /* Hermite */
      x[0] = 0.0 ;
      //#ifdef HERMITE_PROB
      //w[0]=1.0;
      //#else
      w[0] = sqrt(DPI) ;
      //#endif
    }
    else if ( kind == 5 ) { /* Jacobi */
      x[0] = (b-a)/(2.0+a+b);
      w[0] = -(4.0+a+b)/(2.0+a+b)*tgamma(2.0+a)*tgamma(2.0+b)/tgamma(2.0+a+b)/2.0*pow(2.0,a+b)
	/(jpolp_gq(1,a,b,x[0])*jpol_gq(2,a,b,x[0])) ;
    }
    else if ( kind == 6 ) { /* Laguerre */
      x[0] = a+1.0 ;
      w[0] = tgamma(2.0+a)*x[0]/(4.0*lgpol_gq(2,a,x[0]));
    }
    return ;
  }

  /* n > 1 -> compute quad points via Golub-Welch methodology */
  double *sdag = new double[n-1];
  if ( kind == 1 ) { /* Legendre */
    for ( int i=0; i<n  ; i++ ) x[i] = 0.0 ;
    for ( int i=0; i<n-1; i++ ) {
      double di = (double) i+1;
      sdag[i] = sqrt(di*di / ( ( 2.0*di-1.0 ) * ( 2.0*di+1.0 ) ));
    }
  }
  else if ( kind == 4 ) { /* Hermite */
    for ( int i=0; i<n  ; i++ ) x[i] = 0.0 ;
    for ( int i=0; i<n-1; i++ )
      //#ifdef HERMITE_PROB
      //sdag[i] = sqrt((double) i+1);
      //#else
      sdag[i] = sqrt(((double) i+1)/2.0);
    //#endif
  }
  else if ( kind == 5 ) { /* Jacobi */
    double ab2 = a*a-b*b;
    double apb = a+b;
    for ( int i=0; i<n  ; i++ ) {
      double dn = 2.0*((double) i + 1.0) ;
      if ( fabs((dn+apb-2.0)*(dn+apb)) < 1.e-10 )
        x[i] = 0.0 ;
      else
        x[i] = - ab2/((dn+apb-2.0)*(dn+apb)) ;
    }

    for ( int i=0; i<n-1; i++ ) {
      double sn = (double) i+1;
      double dn = 2.0*sn ;
      double numer = 4.0*sn*(sn+a)*(sn+b)*(sn+apb) ;
      double denom = pow(dn+apb,2)*(pow(dn+apb,2)-1.0) ;
      sdag[i] = sqrt( numer/denom );
    }
  }
  else if ( kind == 6 ) { /* Laguerre */
    for ( int i=0; i<n  ; i++ ) x[i] = 2.0*((double) i)+a+1.0 ;
    for ( int i=0; i<n-1; i++ ) {
      double di = (double) i+1;
      sdag[i] = sqrt((di+a) * di );
    }
  }

  int info ;
  int nloc = n, ldz=1 ;
  //cout<<"calling dstev"<<endl;
  //dstev_ ( (char *) "N", &nloc, x, sdag, NULL, &ldz, NULL, &info ) ;
  FTN_NAME(dsteqr)( (char *) "N", &nloc, x, sdag, NULL, &ldz, NULL, &info );

  if ( info != 0 ) {
    cout<<"Error in gq(): dstev returned error :"<<info<<endl<<flush ;
    terminate() ;
  }

  /* Compute weights */
  if ( kind == 1 ) { /* Legendre */
    double dn = (double) n;
    for ( int i=0; i<n  ; i++ ) {
      w[i] = 2.0*(1.0-x[i]*x[i])/pow((dn+1.0)*lpol_gq(n+1,x[i]),2) ;
    }
  }
  else if ( kind == 4 ) { /* Hermite */
    double dn = (double) n;
      for ( int i=0; i<n  ; i++ ){
      //#ifdef HERMITE_PROB
      //w[i] = fact_gq(n)*sqrt(2.0*DPI)/pow(dn*hpol_gq(n-1,x[i]),2);
      //#else
      //w[i] = pow(2.0,n-1)*fact_gq(n)*sqrt(DPI)/pow(dn*hpol_phys_gq(n-1,x[i]),2);
          // Can we have some comments on what the reason for these tests is?
          if (fabs(hpol_phys_gq(n-1,x[i])) != fabs(hpol_phys_gq(n-1,x[i])) or std::isinf(fabs(hpol_phys_gq(n-1,x[i]))))
              w[i]=0.0;
          else{
              w[i] = (n-1.)*log(2.)+logfactorial(n)+0.5*log(DPI)-2.*log(dn)-2.*log(fabs(hpol_phys_gq(n-1,x[i])));
              w[i]=exp(w[i]);
          }
      }
    //#endif
  }
  else if ( kind == 5 ) { /* Jacobi */
    double dn = (double) n;
    double wfac = -(2*dn+a+b+2.0)/(dn+a+b+1.0)*tgamma(n+a+1.0)*tgamma(n+b+1)/
      (fact_gq(n+1)*tgamma(n+a+b+1))*pow(2.0,a+b);
    for ( int i=0; i<n  ; i++ )
      w[i] = wfac/(jpolp_gq(n,a,b,x[i])*jpol_gq(n+1,a,b,x[i]));
  }
  else if ( kind == 6 ) { /* Laguerre */
    double dn = (double) n;
    double wfac = tgamma(n+a+1.0)/(fact_gq(n)*(dn+1.0)*(dn+1.0));
    for ( int i=0; i<n  ; i++ )
      w[i] = wfac*x[i]/pow(lgpol_gq(n+1,a,x[i]),2);
  }

  delete [] sdag ;

  return ;

}

void gchb(const int kind, const int n, double *x, double *w )
{

  if ( kind == 1 )
  {
    double dtheta = DPI/((double) n) ;
    double theta  = 0.5*dtheta ;
    for ( int i = 0 ; i < n ; i++ )
    {
      w[n-1-i] = dtheta ;
      x[n-1-i] = cos(theta) ;
      theta += dtheta ;
    }
  }
  else if ( kind == 2 )
  {
    double dtheta = DPI/((double) n + 1.0 ) ;
    double theta   = dtheta ;
    for ( int i = 0 ; i < n ; i++ )
    {
      w[n-1-i] = dtheta*pow(sin(theta),2) ;
      x[n-1-i] = cos(theta) ;
      theta += dtheta ;
    }
  }

  else
  {
    cout<<"ERROR in gchb() : kind should be either 1 or 2 : "<<kind<<endl<<flush ;
    terminate() ;
  }

  return ;

}

/*
  Gauss Quadrature for generic recursions
 */
void gq_gen(Array1D<double>& a, Array1D<double>& b, const double amu0,
            Array1D<double>& x, Array1D<double>& w)
{
  int i;
  int n = (int) a.XSize();
  int ldz = n;
  int info ;

  /* take sqrt of the off-diagonal */
  for (i=1;i<n;i++) b(i)=sqrt(b(i));
  for (i=0;i<n;i++) x(i)=a(i);

  double *e    = new double[n-1];
  for (i=0;i<n-1;i++) e[i]=b(i+1);

  double *evec = new double[n*n];
  double *work = new double[2*n];

  FTN_NAME(dstev)( (char *) "V", &n, x.GetArrayPointer(), e, evec, &ldz, work, &info ) ;
  if ( info != 0 ) {
    cout<<"Error in gq_gen(): dstev returned error: "<<info<<endl<<flush ;
    terminate() ;
  }

  for ( int i=0; i<n ; i++ )
    w(i) = amu0*evec[i*n]*evec[i*n] ;

  delete [] e;
  delete [] evec;
  delete [] work;

}

void vandermonde_gq(Array1D<double>& x, Array1D<double>& w, Array1D<double>& q) {

  int n = x.XSize() ;

  for (int i=0; i<n; i++) w(i) = q(i) ;

  for ( int k=0; k<n-1; k++ )
    for (int i=n-1; i>=k+1; i--)
      w(i) = w(i)-x(k)*w(i-1);

  for ( int k=n-2; k>=0; k-- ) {
    for (int i=k+1; i<n; i++)
      w(i) = w(i)/(x(i)-x(i-k-1));
    for (int i=k; i<n-1; i++)
      w(i) = w(i)-w(i+1);
  }

  return ;

}

/* Local recursions for orthogonal polynomials */

/* Legendre */
double lpol_gq(int n, double x) {
  if ( n==-1) return (0.0) ;
  if ( n==0 ) return (1.0) ;
  double pn2 = 0.0 ;
  double pn1 = 1.0 ;
  for ( int i=1; i<=n; i++ ) {
    double dn = (double) i ;
    double pn =  ((2.0*dn-1.0)*x*pn1-(dn-1.0)*pn2 )/dn ;
    pn2 = pn1 ;
    pn1 = pn  ;
  }
  //return ( ((2.0*dn-1.0)*x*lpol_gq(n-1,x)-(dn-1.0)*lpol_gq(n-2,x))/dn ) ;
  return ( pn1 ) ;
}

/* Hermite - probabilist */
double hpol_gq(int n, double x) {
  if ( n==-1) return (0.0) ;
  if ( n==0 ) return (1.0) ;
  double pn2 = 0.0 ;
  double pn1 = 1.0 ;
  for ( int i=1; i<=n; i++ ) {
    double pn = x*pn1-((double) (i-1))*pn2 ;
    pn2 = pn1 ;
    pn1 = pn  ;
  }
  //return (x*hpol_gq(n-1,x)-((double) (n-1))*hpol_gq(n-2,x)) ;
  return ( pn1 ) ;
}

/* Hermite - physicist */
double hpol_phys_gq(int n, double x) {
  if ( n==-1) return (0.0) ;
  if ( n== 0) return (1.0) ;
  double pn2 = 0.0 ;
  double pn1 = 1.0 ;
  for ( int i=1; i<=n; i++ ) {
    double pn = 2.0*x*pn1-2.0*((double) (i-1))*pn2 ;
    pn2 = pn1 ;
    pn1 = pn  ;
  }
  //return (2.0*x*hpol_phys_gq(n-1,x)-2.0*((double) (n-1))*hpol_phys_gq(n-2,x)) ;
  return ( pn1 ) ;
}



/* Jacobi */
double jpol_gq(int n, double a, double b, double x) {
  if ( n==-1) return (0.0) ;
  if ( n==0 ) return (1.0) ;
  double pn2 = 0.0 ;
  double pn1 = 1.0 ;
  for ( int i=1; i<=n; i++ ) {
    double dn = (double) i ;
    double denom = 2.0*dn*(dn+a+b)*(2.0*dn+a+b-2.0) ;
    double an = (2.0*dn+a+b-1.0)*(2.0*dn+a+b)*(2.0*dn+a+b-2.0)/denom ;
    double bn = (2.0*dn+a+b-1.0)*(a*a-b*b)/denom ;
    double cn = 2.0*(dn+a-1.0)*(dn+b-1.0)*(2.0*dn+a+b)/denom ;
    double pn = (an*x+bn)*pn1-cn*pn2 ;
    pn2 = pn1 ;
    pn1 = pn  ;
  }
  //return ((an*x+bn)*jpol_gq(n-1,a,b,x)-cn*jpol_gq(n-2,a,b,x)) ;
  return ( pn1 ) ;
}

/* Derivative of Jacobi polynomials */
double jpolp_gq(int n, double a, double b, double x) {
  if ( n==0 ) return (0.0) ;
  double dn = (double) n ;
  return ( 0.5*(dn+a+b+1.0)*jpol_gq(n-1,a+1.0,b+1.0,x));
}

/* Laguerre */
double lgpol_gq(int n, double a, double x) {
  if ( n==0 ) return (1.0) ;
  if ( n==1 ) return (-x+a+1.0) ;
  double pn2 = 0.0 ;
  double pn1 = 1.0 ;
  for ( int i=1; i<=n; i++ ) {
    double dn = (double) i ;
    double pn =  ((-x+2.0*dn+a-1.0)*pn1-(dn+a-1.0)*pn2)/dn ;
    pn2 = pn1 ;
    pn1 = pn  ;
  }
  //return ((-x+2.0*dn+a-1.0)*lgpol_gq(n-1,a,x)-(dn+a-1.0)*lgpol_gq(n-2,a,x))/dn ;
  return ( pn1 );
}

double fact_gq(int n) {
  return tgamma(n + 1);
}
