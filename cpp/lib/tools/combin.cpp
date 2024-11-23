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



/*
Computes lower incomplete gamma integral, based on rewriting
https://www.cfa.harvard.edu/sma/miriad/wbcorrTest/downLoad/NewFormat/miriad-sma3.1.2/src/subs/gamma.f
in C.
*/
double gammai ( const double p, const double x )
{
//
//  Check the input.
//
  if ( x == 0.0 )
    return 0.0;

  if ( x < 0.0 )
  {

    throw Tantrum("gammai() error:: the argument of the incomplete gamma function is non-positive. Exiting.");
    exit(1);
  }

  if ( p <= 0.0 )
  {
    throw Tantrum("gammai() error:: the parameter of the incomplete gamma function is non-positive. Exiting.");
    exit(1);
  }

  if (x <= p + 1.0)
  {

    double ap = p;
    double sum = 1./p;
    double del = sum;
    while (fabs(del) > fabs(sum)*1.0e-12){
      ap += 1.0;
      del *= (x/ap);
      sum += del;
      }
    return sum * exp(-x + p*log(x) - lgamma(p));
  }

  else
  {
      double gold = 0.0;
      double a0 = 1.0;
      double a1 = x;
      double b0 = 0.0;
      double b1 = 1.0;
      double fac = 1.0;

      for (int i=1; i<=100; i++){
        double ana = float(i) - p;
        a0 = (a1 + a0*ana) * fac;
        b0 = (b1 + b0*ana) * fac;
        double anf = float(i) * fac;
        a1 = x*a0 + anf*a1;
        b1 = x*b0 + anf*b1;
        if (a1!=0){
          fac = 1./a1;
          double g = b1 * fac;
          if (abs((g-gold)/g) <= 1.0e-12)
            return 1.- g * exp(-x + p*log(x) - lgamma(p));
          else
            gold = g;
        }

      }


  }

}


// Beta function
double beta(const double z, const double w)
{
  return exp(lgamma(z)+lgamma(w)-lgamma(z+w));
}




/*
Implementation of incomplete beta function using https://codeplea.com/incomplete-beta-function-c (the code itself is in https://github.com/codeplea/incbeta/blob/master/incbeta.c)
*/
/*
 * zlib License
 *
 * Regularized Incomplete Beta Function
 *
 * Copyright (c) 2016, 2017 Lewis Van Winkle
 * http://CodePlea.com
 *
 * This software is provided 'as-is', without any express or implied
 * warranty. In no event will the authors be held liable for any damages
 * arising from the use of this software.
 *
 * Permission is granted to anyone to use this software for any purpose,
 * including commercial applications, and to alter it and redistribute it
 * freely, subject to the following restrictions:
 *
 * 1. The origin of this software must not be misrepresented; you must not
 *    claim that you wrote the original software. If you use this software
 *    in a product, an acknowledgement in the product documentation would be
 *    appreciated but is not required.
 * 2. Altered source versions must be plainly marked as such, and must not be
 *    misrepresented as being the original software.
 * 3. This notice may not be removed or altered from any source distribution.
 */
#define STOP 1.0e-8
#define TINY 1.0e-30

double betai(double a, double b, double x) {
//
//  Check the input arguments.
//
  if ( a <= 0.0 || b <= 0.0 )
  {
    throw Tantrum("betai() error:: parameters of the incomplete beta function are negative. Exit.");
    exit(1);

    return 0.0;
  }

  if ( x < 0.0 || 1.0 < x )
  {
    throw Tantrum("betai() error:: argument of the incomplete beta function are outside bounds [0,1]. Exit.");
    exit(1);

    return 0.0;
  }
 //   if (x < 0.0 || x > 1.0) return 1.0/0.0; //give error

    /*The continued fraction converges nicely for x < (a+1)/(a+b+2)*/
    if (x > (a+1.0)/(a+b+2.0)) {
        return (1.0-betai(b,a,1.0-x)); /*Use the fact that beta is symmetrical.*/
    }

    /*Find the first part before the continued fraction.*/
    const double lbeta_ab = lgamma(a)+lgamma(b)-lgamma(a+b);
    const double front = exp(log(x)*a+log(1.0-x)*b-lbeta_ab) / a;

    /*Use Lentz's algorithm to evaluate the continued fraction.*/
    double f = 1.0, c = 1.0, d = 0.0;

    int i, m;
    for (i = 0; i <= 200; ++i) {
        m = i/2;

        double numerator;
        if (i == 0) {
            numerator = 1.0; /*First numerator is 1.0.*/
        } else if (i % 2 == 0) {
            numerator = (m*(b-m)*x)/((a+2.0*m-1.0)*(a+2.0*m)); /*Even term.*/
        } else {
            numerator = -((a+m)*(a+b+m)*x)/((a+2.0*m)*(a+2.0*m+1)); /*Odd term.*/
        }

        /*Do an iteration of Lentz's algorithm.*/
        d = 1.0 + numerator * d;
        if (fabs(d) < TINY) d = TINY;
        d = 1.0 / d;

        c = 1.0 + numerator / c;
        if (fabs(c) < TINY) c = TINY;

        const double cd = c*d;
        f *= cd;

        /*Check for stop.*/
        if (fabs(1.0-cd) < STOP) {
            return front * (f-1.0);
        }
    }

    return 1.0/0.0; /*Needed more loops, did not converge.*/
}



/*
Calculates di-gamma function, also called psi-function, defined as the derivative of log-gamma function: psi(x) = d log (Gamma(x))/d x
The calculation relies on analytic derivative of the Lanczos approximation of the Gamma function https://en.wikipedia.org/wiki/Lanczos_approximation
Note that the function in principle can be defined for non-integer negative values (e.g. see https://www.codecogs.com/library/maths/special/gamma/psi.php), but here we constrain the domain to be the positive semi-axis only.
Also, inspired by this Matlab implementation https://www.mathworks.com/matlabcentral/fileexchange/978-special-functions-math-library, which did not work correctly as is.
*/
double digama ( double x )
{
  if ( x <= 0.0 )
  {
    throw Tantrum("digama() error:: argument of the digamma function is non-positive. Exit.");
    exit(1);
  }

  if (x<0.5)
      return digama(1.-x) - 4.*atan(1.)/tan(4.*atan(1.)*x);

  double g=607./128.; // best results when 4<=g<=5

double c[15] = {0.99999999999999709182,
      57.156235665862923517,
     -59.597960355475491248,
      14.136097974741747174,
      -0.49191381609762019978,
        .33994649984811888699e-4,
        .46523628927048575665e-4,
       -.98374475304879564677e-4,
        .15808870322491248884e-3,
       -.21026444172410488319e-3,
        .21743961811521264320e-3,
       -.16431810653676389022e-3,
        .84418223983852743293e-4,
       -.26190838401581408670e-4,
        .36899182659531622704e-5};

        double n=0.;
        double d=0.;
double dz = 0.0;
double dd = 0.0;

        for (int k=14; k>=1; k--){
          dz = 1./(x+float(k));
          dd = c[k]*dz;
          d += dd;
          n -= dd*dz;
        }
        d += c[0];
        double gg = x + g +0.5;
        double f = log(gg) + (n/d - g/gg) - 1./x;



        return f;
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
    explVariance=1.-unexplVariance/(totVariance+1.e-16);

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

