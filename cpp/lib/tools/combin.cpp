/* =====================================================================================
                     The UQ Toolkit (UQTk) version 3.0.4
                     Copyright (2017) Sandia Corporation
                     http://www.sandia.gov/UQToolkit/

     Copyright (2017) Sandia Corporation. Under the terms of Contract DE-AC04-94AL85000
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
/// \file combin.cpp
/// \brief Tools to evaluate combinatorial quantities.

#include "Array1D.h"
#include "Array2D.h"
#include "gen_defs.h"
#include "probability.h"
#include "combin.h"
#include <math.h>
#include <float.h>
#include "error_handlers.h"

// Calculates binomial formula n-choose-k
int choose(int n,int k) {

  int vmin=MIN(k,n-k);

  if (vmin< 0) return (0);
  if (vmin==0) return (1);

  int vmax = MAX(k,n-k);
  int vret = vmax+1;
  for (int i=2;i<vmin+1;i++)
    vret = ( vret * (vmax+i) ) / i;

  return ( vret );

}


// Calculates factorial of a number
int factorial(int number) {
	int temp;

	if(number <= 1) return 1;

	temp = number * factorial(number - 1);
	return temp;
}

// Calculates log-factorial of a number
double logfactorial(int number) {
  
	double temp;
    
	if(number <= 1) return 0;
    
	temp = log(number) + logfactorial(number - 1);
	return temp;
}

// Computes all combinations n-choose-k and lists them all in an array
void chooseComb(int n, int k,Array2D<int>& fullInd)
{
  int n_k=choose(n,k);
  fullInd.Resize(n_k,k,0);

  for(int ik=0;ik<k;ik++)
    fullInd(0,ik)=ik;

  fullInd(0,k-1)=k-2;

  int iii=0;
  int j=k-1;
  while(j>=0){
    if(fullInd(iii,j)<n-k+j){
        fullInd(iii,j)++;
        for(int jj=j+1;jj<k;jj++)
          fullInd(iii,jj)=fullInd(iii,j)+jj-j;
        
        j=k-1;
        iii++;
        if (iii==(int)fullInd.XSize()) {j=-1; break;}
        for(int ik=0;ik<k;ik++)
          fullInd(iii,ik)=fullInd(iii-1,ik);
    }
    else j--;

    
  }

  return;
}

// Permutes a given array in place
void get_perm(Array1D<int>& perm,int seed)
{

  int nn=perm.XSize();
  get_perm(nn,perm.GetArrayPointer(),seed);

  return;
}

// Permutes a given array in place
void get_perm(int nn, int* perm,int seed)
{

  dsfmt_gv_init_gen_rand(seed );
   
  int j,t;
  for(int is=0;is<nn;is++)
    perm[is]=is;

  for(int is=0;is<nn;is++){
    j=is+(int) (dsfmt_gv_genrand_urv()*(nn-is));
    t=perm[is];
    perm[is]=perm[j];
    perm[j]=t;
  }     

  return;
}



// Computes incomplete Gamma integral
double gammai ( const double p, const double x )
//****************************************************************************80
//
//  Purpose:
//
//    GAMMDS computes the incomplete Gamma integral.
//
//  Discussion:
//
//    The parameters must be positive.  An infinite series is used.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 January 2008
//
//  Author:
//
//    Original FORTRAN77 version by Chi Leung Lau.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Chi Leung Lau,
//    Algorithm AS 147:
//    A Simple Series for the Incomplete Gamma Integral,
//    Applied Statistics,
//    Volume 29, Number 1, 1980, pages 113-114.
//
//  Parameters:
//
//    Input, double X, P, the arguments of the incomplete
//    Gamma integral.  X and P must be greater than 0.
//
//    Output, int *IFAULT, error flag.
//    0, no errors.
//    1, X <= 0 or P <= 0.
//    2, underflow during the computation.
//
//    Output, double GAMMDS, the value of the incomplete
//    Gamma integral.
/////////////////////////////////////////////////////////////////////////////////////////
// The code is taken from http://people.sc.fsu.edu/~jburkardt/cpp_src/asa147/asa147.html
// UQTk group has modified the code in the following way:
// - Renamed this function from gammds to gammai for backward compatibility
// - Instead of integer fault indicator, we throw Tantrum and exit when parameters or arguments are outside bounds
// - Log of complete gamma function is computed by lgamma instead of being computed by another function from the same code-suite
// - Instead of an underflow error-message, implemented an approximation (i.e. returning 1) for large arguments when the output value is near 1
{
  double a;
  double arg;
  double c;
  double e = 1.0E-09;
  double f;
  //int ifault2;
  double uflo = 1.0E-37;
  double value;
//
//  Check the input.
//
  if ( x <= 0.0 )
  {
    throw Tantrum("gammai() error:: the argument of the incomplete gamma function is non-positive. Exit.");
    value = 0.0;
    return value;
  }

  if ( p <= 0.0 ) 
  {
    throw Tantrum("gammai() error:: the parameter of the incomplete gamma function is non-positive. Exit.");
    value = 0.0;
    return value;
  }
//
//
  arg = p * log ( x ) - lgamma ( p + 1.0 ) - x;

  

  f = exp ( arg );

  if ( f < uflo )
  {
    // This is added by UQTk group as an aproximation in the case of large arguments
    value = 1.0;
    return value;
  }

//
//  Series begins.
//
  c = 1.0;
  value = 1.0;
  a = p;

  for ( ; ; )
  {
    a = a + 1.0;
    c = c * x / a;
    value = value + c;

    if ( c <= e * value )
    {
      break;
    }
  }

  value = value * f;

  return value;
}

// Beta function
double beta(const double z, const double w)
{
  return exp(lgamma(z)+lgamma(w)-lgamma(z+w));
}



// Computes incomplete Beta function
double betai ( const double p, const double q, const double x )
//****************************************************************************
//
//  Purpose:
//
//    BETAIN computes the incomplete Beta function ratio.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 January 2008
//
//  Author:
//
//    Original FORTRAN77 version by KL Majumder, GP Bhattacharjee.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    KL Majumder, GP Bhattacharjee,
//    Algorithm AS 63:
//    The incomplete Beta Integral,
//    Applied Statistics,
//    Volume 22, Number 3, 1973, pages 409-411.
//
//  Parameters:
//
//    Input, double X, the argument, between 0 and 1.
//
//    Input, double P, Q, the parameters, which
//    must be positive.
//
//    Input, double BETA, the logarithm of the complete
//    beta function.
//
//    Output, int *IFAULT, error flag.
//    0, no error.
//    nonzero, an error occurred.
//
//    Output, double BETAIN, the value of the incomplete
//    Beta function ratio.
/////////////////////////////////////////////////////////////////////////////////////////
// The code is taken from http://people.sc.fsu.edu/~jburkardt/cpp_src/asa063/asa063.html
// UQTk group has modified the code in the following way:
// - Renamed this function from betain to betai for backward compatibility
// - Instead of integer fault indicator, we throw Tantrum and exit when parameters or arguments are outside bounds
// - Log of complete beta function is computed instead of being given as an argument
//
{
  double acu = 0.1E-14;
  double ai;
  //double betain;
  double cx;
  bool indx;
  int ns;
  double pp;
  double psq;
  double qq;
  double rx;
  double temp;
  double term;
  double value;
  double xx;




  value = x;
//
//  Check the input arguments.
//
  if ( p <= 0.0 || q <= 0.0 )
  {
    throw Tantrum("betai() error:: parameters of the incomplete beta function are negative. Exit.");

    return 0.0;
  }

  if ( x < 0.0 || 1.0 < x )
  {
    throw Tantrum("betai() error:: argument of the incomplete beta function are outside bounds [0,1]. Exit.");
    
    return 0.0;
  }
//
//  Special cases.
//
  if ( x == 0.0 || x == 1.0 )
  {
    return value;
  }

  // Added by Sandia UQTk group
  double lbeta=lgamma(p)+lgamma(q)-lgamma(p+q);

//
//  Change tail if necessary and determine S.
//
  psq = p + q;
  cx = 1.0 - x;

  if ( p < psq * x )
  {
    xx = cx;
    cx = x;
    pp = q;
    qq = p;
    indx = true;
  }
  else
  {
    xx = x;
    pp = p;
    qq = q;
    indx = false;
  }

  term = 1.0;
  ai = 1.0;
  value = 1.0;
  ns = ( int ) ( qq + cx * psq );
//
//  Use the Soper reduction formula.
//
  rx = xx / cx;
  temp = qq - ai;
  if ( ns == 0 )
  {
    rx = xx;
  }

  for ( ; ; )
  {
    term = term * temp * rx / ( pp + ai );
    value = value + term;;
    temp = fabs ( term );

    if ( temp <= acu && temp <= acu * value )
    {
      value = value * exp ( pp * log ( xx ) 
      + ( qq - 1.0 ) * log ( cx ) - lbeta ) / pp;

      if ( indx )
      {
        value = 1.0 - value;
      }
      break;
    }

    ai = ai + 1.0;
    ns = ns - 1;

    if ( 0 <= ns )
    {
      temp = qq - ai;
      if ( ns == 0 )
      {
        rx = xx;
      }
    }
    else
    {
      temp = psq;
      psq = psq + 1.0;
    }
  }

  return value;
}

// Calculates di-gamma function
double digama ( double x )
//****************************************************************************80
//
//  Purpose:
//
//    DIGAMA calculates DIGAMMA ( X ) = d ( LOG ( GAMMA ( X ) ) ) / dX
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 June 2013
//
//  Author:
//
//    Original FORTRAN77 version by Jose Bernardo.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Jose Bernardo,
//    Algorithm AS 103:
//    Psi ( Digamma ) Function,
//    Applied Statistics,
//    Volume 25, Number 3, 1976, pages 315-317.
//
//  Parameters:
//
//    Input, double X, the argument of the digamma function.
//    0 < X.
//
//    Output, int *IFAULT, error flag.
//    0, no error.
//    1, X <= 0.
//
//    Output, double DIGAMA, the value of the digamma function at X.
//
/////////////////////////////////////////////////////////////////////////////////////////
// The code is taken from http://people.sc.fsu.edu/~jburkardt/cpp_src/asa103/asa103.cpp
// UQTk group has modified the code in the following way:
// - Instead of integer fault indicator, we throw Tantrum and exit when parameters or arguments are outside bounds
//
{
  double euler_mascheroni = 0.57721566490153286060;
  double r;
  double value;
  double x2;
//
//  Check the input.
//
  if ( x <= 0.0 )
  {
    value = 0.0;
    throw Tantrum("digama() error:: argument of the digamma function is non-positive. Exit.");
    return value;
  }
//
//  Initialize.
//
  x2 = x;
  value = 0.0;
//
//  Use approximation for small argument.
//
  if ( x2 <= 0.00001 )
  {
    value = - euler_mascheroni - 1.0 / x2;
    return value;
  }
//
//  Reduce to DIGAMA(X + N).
//
  while ( x2 < 8.5 )
  {
    value = value - 1.0 / x2;
    x2 = x2 + 1.0;
  }
//
//  Use Stirling's (actually de Moivre's) expansion.
//
  r = 1.0 / x2;
  value = value + log ( x2 ) - 0.5 * r;
  r = r * r;
  value = value 
    - r * ( 1.0 / 12.0
    - r * ( 1.0 / 120.0 
    - r *   1.0 / 252.0 ) );

  return value;
}
//****************************************************************************


// K-center clustering of data
void clust(Array2D<double>& data_in, Array1D<double>& w,int ncl, Array1D<int>& numData,int *pClusterIndex)
{
  int nsample=data_in.XSize();
  int ndim=data_in.YSize();
  double *pSources;
  pSources=new double[ndim*nsample];

  for (int i=0;i<ndim*nsample;i++){
    int i_dim=i%ndim;
    int id=(int) (i-i_dim)/ndim;
 
    pSources[i]=data_in(id,i_dim)/w(i_dim);
  }

  KCenterClustering kcenter(ndim,nsample,pSources,pClusterIndex,ncl);
  kcenter.Cluster();

  numData.Resize(ncl,0);

  for (int i=0;i<nsample;i++)
    numData(pClusterIndex[i])++;
  

  delete []pSources;

  return;
}

// Multiple trials of K-center clustering and picking the best one according to explained variance criterion
double clust_best(Array2D<double>& data_in, Array1D<double>& w,int ncl, Array1D<int>& bestnumData,int *bestClusterIndex,int ntry)
{

  int nsample=data_in.XSize();
  int ndim=data_in.YSize();

  double bestExplVar=0;
  Array1D<double> mean(ndim,0.e0);
  double totVariance=getMean_Variance(data_in,w,mean);
  

  for (int itry=0;itry<ntry;itry++){
     
    Array1D<int> numData(ncl,0);   
    int *pClusterIndex;
    pClusterIndex=new int[nsample];
     
    clust(data_in,w,ncl,numData,pClusterIndex);

    double unexplVariance=0.0, explVariance=0.0;

    for (int icl=0;icl<ncl;icl++){
      Array2D<double> data_icl(numData(icl),ndim,0.e0);
      int jc=0;
      
      for(int i=0;i<nsample;i++){
        if (pClusterIndex[i]==icl) {
          for (int idim=0;idim<ndim;idim++)
            data_icl(jc,idim)=data_in(i,idim);
     
          jc++;
        }
      }
       
      Array1D<double> mean_cl(ndim,0.e0);
      unexplVariance += (numData(icl)*getMean_Variance(data_icl,w,mean_cl));
    }
     
    unexplVariance /= nsample;
    explVariance=1.-unexplVariance/totVariance;

    if (explVariance >= bestExplVar) {
      for(int i=0;i<nsample;i++){    bestClusterIndex[i]=pClusterIndex[i];      }
      for(int icl=0;icl<ncl;icl++){    bestnumData(icl)=numData(icl);     }
      bestExplVar=explVariance;
       //printf("Trial #%d, bestExplVarFrac=%lg\n",itry+1,bestExplVar);
    }
     
    delete []pClusterIndex;

  }
  
  return bestExplVar;
}

// Find the best number of clusters in a dataset according to one of three (hardcoded) criteria
int findNumCl(Array2D<double>& data_in,Array1D<double>& w,int ntry)// Below one can uncomment either a) or b) or c) for three differnet elbow-like criteria
{
  int maxNumCl=10;
  //double explThresh=0.5;
  int optNumCl=maxNumCl;
  double bestExplVar=0.;
  //double oldJump=0.;
  double newJump=0.;
  double prev_bestExplVar=0.;
  //double minSlopeDiff=0.;
  double maxJump=0.;

  //int ndim=data_in.YSize();//=w.XSize();
  int nsam=data_in.XSize();
 
  for(int tryNumCl=1; tryNumCl<=maxNumCl;tryNumCl++){
    Array1D<int> numData(tryNumCl,0);
    int *pClusterIndex;
    pClusterIndex=new int[nsam];
    bestExplVar=clust_best(data_in,w,tryNumCl,numData,pClusterIndex,ntry);
    
    delete []pClusterIndex;
    printf("bestExplVarFrac(%d)=%lg\n",tryNumCl,bestExplVar);


    //a)
    //    if (bestExplVar > explThresh) {optNumCl=tryNumCl; break;}

   //b)
   newJump=bestExplVar-prev_bestExplVar;
   if ( newJump > maxJump) {maxJump=newJump; optNumCl=tryNumCl;}
   prev_bestExplVar=bestExplVar;
    
   //c)
   //newJump=bestExplVar-prev_bestExplVar;
   //  if ( newJump-oldJump < minSlopeDiff) {minSlopeDiff=newJump-oldJump; optNumCl=tryNumCl-1;}
   //     oldJump=newJump;
   //     prev_bestExplVar=bestExplVar;
    
  }
  
  return optNumCl;
}

