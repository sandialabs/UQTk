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
/** \file probability.cpp
 * \brief Probability and random number generation- related tools.
 */
#include <math.h>
#include <iostream>
#include <float.h>
#include <limits.h>

#include "Array1D.h"
#include "Array2D.h"
#include "probability.h"
#include "rosenblatt.h"
#include "combin.h"
#include "arraytools.h"
#include "gen_defs.h"

#include "minmax.h"

#ifndef M_PI
#define M_PI atan(1.0) * 4.0
#endif

using namespace std;

// Error function via incomplete gamma function
double erff(const double x)
{
  return x < 0.0 ? -gammai(0.5,x*x) : gammai(0.5,x*x);
}

// Inverse error function, input scaled to [-1,1]
double inverf(double y0)
{
  double result_;

    double expm2;
    double s2pi;
    double x;
    double y;
    double z;
    double y2;
    double x0;
    double x1;
    int code;
    double p0;
    double q0;
    double p1;
    double q1;
    double p2;
    double q2;


    y0=0.5*(y0+1);

    expm2 = 0.13533528323661269189;
    s2pi = 2.50662827463100050242;
    if( y0<=0 or y0 >=1)
    {
      cout << "Error in inverting erf, the argument should be -1<x<1:" << 2*y0-1 << endl;
      return 0;
      exit(1);
    }


    code = 1;
    y = y0;
    if( y>1.0-expm2 )
    {
        y = 1.0-y;
        code = 0;
    }
    if( y>expm2 )
    {
        y = y-0.5;
        y2 = y*y;
        p0 = -59.9633501014107895267;
        p0 = 98.0010754185999661536+y2*p0;
        p0 = -56.6762857469070293439+y2*p0;
        p0 = 13.9312609387279679503+y2*p0;
        p0 = -1.23916583867381258016+y2*p0;
        q0 = 1;
        q0 = 1.95448858338141759834+y2*q0;
        q0 = 4.67627912898881538453+y2*q0;
        q0 = 86.3602421390890590575+y2*q0;
        q0 = -225.462687854119370527+y2*q0;
        q0 = 200.260212380060660359+y2*q0;
        q0 = -82.0372256168333339912+y2*q0;
        q0 = 15.9056225126211695515+y2*q0;
        q0 = -1.18331621121330003142+y2*q0;
        x = y+y*y2*p0/q0;
        x = x*s2pi;
        result_ = x;
        return result_/sqrt(2.);
    }
    x = sqrt(-2.0*log(y));
    x0 = x-log(x)/x;
    z = 1.0/x;
    if( x<8.0 )
    {
        p1 = 4.05544892305962419923;
        p1 = 31.5251094599893866154+z*p1;
        p1 = 57.1628192246421288162+z*p1;
        p1 = 44.0805073893200834700+z*p1;
        p1 = 14.6849561928858024014+z*p1;
        p1 = 2.18663306850790267539+z*p1;
        p1 = -1.40256079171354495875*0.1+z*p1;
        p1 = -3.50424626827848203418*0.01+z*p1;
        p1 = -8.57456785154685413611*0.0001+z*p1;
        q1 = 1;
        q1 = 15.7799883256466749731+z*q1;
        q1 = 45.3907635128879210584+z*q1;
        q1 = 41.3172038254672030440+z*q1;
        q1 = 15.0425385692907503408+z*q1;
        q1 = 2.50464946208309415979+z*q1;
        q1 = -1.42182922854787788574*0.1+z*q1;
        q1 = -3.80806407691578277194*0.01+z*q1;
        q1 = -9.33259480895457427372*0.0001+z*q1;
        x1 = z*p1/q1;
    }
    else
    {
        p2 = 3.23774891776946035970;
        p2 = 6.91522889068984211695+z*p2;
        p2 = 3.93881025292474443415+z*p2;
        p2 = 1.33303460815807542389+z*p2;
        p2 = 2.01485389549179081538*0.1+z*p2;
        p2 = 1.23716634817820021358*0.01+z*p2;
        p2 = 3.01581553508235416007*0.0001+z*p2;
        p2 = 2.65806974686737550832*0.000001+z*p2;
        p2 = 6.23974539184983293730*0.000000001+z*p2;
        q2 = 1;
        q2 = 6.02427039364742014255+z*q2;
        q2 = 3.67983563856160859403+z*q2;
        q2 = 1.37702099489081330271+z*q2;
        q2 = 2.16236993594496635890*0.1+z*q2;
        q2 = 1.34204006088543189037*0.01+z*q2;
        q2 = 3.28014464682127739104*0.0001+z*q2;
        q2 = 2.89247864745380683936*0.000001+z*q2;
        q2 = 6.79019408009981274425*0.000000001+z*q2;
        x1 = z*p2/q2;
    }
    x = x0-x1;
    if( code!=0 )
    {
        x = -x;
    }
    result_ = x;

    return result_/sqrt(2.);
}

// Inverse of standard normal CDF
double invnormcdf(double y)
{
  // Scaled inverse error function, really
  return sqrt(2.)*inverf(2.*y-1.);
}

// Standard normal CDF
double normcdf(double y)
{
  // Scaled and shifted error function
  return 0.5*(1.+erf(y/sqrt(2.)));
}

// Complementary function for standard normal CDF
double normcdfc(double y)
{
  return 1.-normcdf(y);
}

// Generate an array of uniform random variables
void generate_uniform(double* rvar,int ns, int nd, int zSeed)
{

  dsfmt_gv_init_gen_rand(zSeed );


  for(int is = 0 ; is < ns*nd ; is++)
    rvar[is]=dsfmt_gv_genrand_urv();

  return;
}

// Generate an array of uniform random variables
void generate_uniform(Array2D<double>& rvar,int zSeed)
{
  int nsample = (int) rvar.XSize();
  int ndim    = (int) rvar.YSize();


  generate_uniform(rvar.GetArrayPointer(), nsample, ndim, zSeed);

  return;
}

// Generate an array of uniform random variables
void generate_uniform(double *rvar, int ns, int nd, dsfmt_t *rnstate)
{
  int i ;
  // Need to check allocation?
  for ( i=0; i<ns*nd; i++ )
      rvar[i] = dsfmt_genrand_open_open( rnstate ) ;

  return ;

}

// Generate an array of uniform random variables
void generate_uniform(Array2D<double> &rvar, dsfmt_t *rnstate)
{
  int nsample = (int) rvar.XSize() ;
  int ndim    = (int) rvar.YSize() ;

  generate_uniform(rvar.GetArrayPointer(),nsample, ndim, rnstate) ;

  return ;

}

// Generate an array of uniform random variables with LHS
void generate_uniform_lhs(double *rvar,int nsample, int ndim, int zSeed)
{

  int *perm = (int *) malloc(nsample*sizeof(int)) ;

  dsfmt_gv_init_gen_rand(zSeed );

  int ii=0;
  for(int id=0;id<ndim;id++){

    int seed=(int) (dsfmt_gv_genrand_urv()*INT_MAX);
    get_perm(nsample, perm, seed);

    for(int is=0;is<nsample;is++){
      double urv=dsfmt_gv_genrand_urv();
      rvar[ii]=(urv+perm[is])/nsample;
      ii++;
    }
  }

  free(perm);

  return;
}

// Generate an array of uniform random variables with LHS
void generate_uniform_lhs(Array2D<double>& rvar,int zSeed)
{
  int ndim=rvar.YSize();
  int nsample=rvar.XSize();

  generate_uniform_lhs(rvar.GetArrayPointer(),nsample,ndim,zSeed);

  return;
}

// Generate an array of uniform random variables with LHS
void generate_uniform_lhs(double *rvar,int nsample, int ndim, dsfmt_t *rnstate)
{

  int *perm = (int *) malloc(nsample*sizeof(int)) ;

  int ii=0;
  for(int id=0;id<ndim;id++){

    int seed = dsfmt_genrand_uint32(rnstate);
    get_perm(nsample,perm,seed);

    for(int is=0;is<nsample;is++){
      double urv=dsfmt_genrand_open_open( rnstate );;
      rvar[ii]=(urv+perm[is])/nsample;
      ii++;
    }
  }

  free(perm);

  return;
}

// Generate an array of uniform random variables with LHS
void generate_uniform_lhs(Array2D<double>& rvar,dsfmt_t *rnstate)
{
  int ndim=rvar.YSize();
  int nsample=rvar.XSize();

  generate_uniform_lhs(rvar.GetArrayPointer(),nsample, ndim,rnstate);

  return;
}

// Generate an array of standard normal random variables
void generate_normal(Array2D<double>& rvar,int zSeed)
{
  int ndim=rvar.YSize();
  int nsample=rvar.XSize();

  dsfmt_gv_init_gen_rand( zSeed );
  for(int is = 0 ; is < nsample ; is++){
    for(int id = 0 ; id < ndim ; id++){
      rvar(is,id)=dsfmt_gv_genrand_nrv();
    }
  }

  return;
}

// Generate an array of standard normal random variables with LHS
void generate_normal_lhs(Array2D<double>& rvar,int zSeed)
{
  int ndim=rvar.YSize();
  int nsample=rvar.XSize();

  generate_uniform_lhs(rvar,zSeed);
  for(int is = 0 ; is < nsample ; is++){
    for(int id = 0 ; id < ndim ; id++){
      rvar(is,id)=invnormcdf(rvar(is,id));
    }
  }

  return;
}

// Returns the median of an array of data
double get_median(const Array1D<double>& data)
{

  int ndata=data.XSize();
  Array1D<double> data_copy;
  data_copy=data;

  double median;
  if (ndata % 2 == 1){
    int k=(int) ndata/2;
    median=select_kth(k,data_copy);
  }
  else{
    int k=ndata/2;
    median=( select_kth(k,data_copy)+select_kth(k-1,data_copy) ) / 2.;
  }
  return median;

}

// Returns the mean of an array of data
double get_mean(const Array1D<double>& data) {

  double mean=0.0;
  int ndata=data.XSize();

  for(int i=0;i<ndata;i++) mean += data(i);
  mean=mean/ndata;

  return mean;

}

// Returns the mean of a 2d array of data
double get_mean(const Array2D<double>& data) {

  double mean=0.0;
  int nrows=data.XSize();
  int ncols=data.YSize();

  for(int i2=0;i2<ncols;i2++)
    for(int i1=0;i1<nrows;i1++)
      mean += data(i1,i2);

  mean = mean/(nrows*ncols);

  return mean;

}

// Returns the standard deviation of an array of data
double get_std(const Array1D<double>& data)
{
  double std;
  double var=get_var(data);

  if (var>=0)
    std=sqrt(var);
  else{
    cout << "probability.cpp:get_std(): negative variance!. Exiting." << endl;
    exit(1);
  }

  return std;

}

// Returns the unbiased estimator of variance of an array of data
double get_var(const Array1D<double>& data)
{
  double mean=get_mean(data);
  double mean2=0.0;

  int ndata=data.XSize();

  for(int i=0;i<ndata;i++) mean2 += (data(i)-mean)*(data(i)-mean);
  mean2=mean2/(ndata-1);

  double var=mean2;
  return var;

}

// Returns the vector mean and variance of a 2d array of data
double getMean_Variance(Array2D<double>& data_c, Array1D<double>& w, Array1D<double>& mean)
{
  double var=0.;
  int ndim=data_c.YSize();
  int nsam=data_c.XSize();

  getMean(data_c,mean);

  Array1D<double> data_1s(ndim,0.e0);

  for (int isam=0;isam<nsam;isam++){
    for (int idim=0;idim<ndim;idim++){
      data_1s(idim)=data_c(isam,idim);
    }
    var+=dist_sq(data_1s,mean,w);
  }
  var /= (nsam); //taking the biased estimator to avoid nsam=1 errors!

  return var;

}

// Computes vector mean, column by column
void getMean(Array2D<double>& data_c, Array1D<double>& mean)
{

  int ndim=data_c.YSize();//=mean.XSize()
  int nsam=data_c.XSize();

	CHECKEQ(mean.XSize(), ndim);

  for (int idim=0; idim<ndim; idim++){

    mean(idim) = 0.0;
    for (int isam=0;isam<nsam;isam++){
      mean(idim) += data_c(isam,idim);
    }
    mean(idim) /= nsam;

  }

  return;

}

// Computes vector mean, either column by column for RC="C" or
// row by row for RC="R"
void getMean(Array2D<double>& matrix, Array1D<double>& mean, char *RC)
{

  int nrows=matrix.XSize();
  int ncols=matrix.YSize();

  if ( std::string(RC) == std::string("C")) {
    mean.Resize(ncols,0.0);
    for (int i2=0; i2<ncols; i2++){
      for (int i1=0; i1<nrows; i1++){
        mean(i2) += matrix(i1,i2);
      }
      mean(i2) /= nrows;
    }
  } else if ( std::string(RC) == std::string("R")) {
		mean.Resize(nrows,0.0);
		for (int i1=0; i1<nrows;i1++){
      for (int i2=0; i2<ncols; i2++){
        mean(i1) += matrix(i1,i2);
      }
      mean(i1) /= ncols;
    }

  } else {
    cout << "probability.cpp:getMean(): unknown flag"
         << RC << endl << flush;
    exit(1);

  }
  return;

}

//  Returns a random permutation of 0..n-1
void rperm(int n, int *a, dsfmt_t *rnstate)
{
  int k;

  for (k = 0; k < n; k++) a[k] = k;

  for ( k = n-1; k > 0; k-- )
  {
    int j = dsfmt_genrand_uint32(rnstate) % (k+1) ;
    int temp = a[j];
    a[j] = a[k];
    a[k] = temp;
  }
    return ;
}

// KDE estimation of a PDF
void getPdf_figtree(Array2D<double>& source,Array2D<double>& target,Array1D<double>& sig, Array1D<double>& density, Array1D<double>& weight)
{
    int NSources=source.XSize();
    int Dim=source.YSize();

    int MTargets=target.XSize();

    Array2D<double> allpts(NSources+MTargets,Dim,0.e0);
    Array1D<double> aa(Dim,0.e0);
    Array1D<double> bb(Dim,0.e0);

    merge(source,target,allpts);
    getDomain(allpts,aa,bb);

    double alpha=sqrt(2.0)*sig(0)/(bb(0)-aa(0));
    for (int i_dim=1;i_dim<Dim;i_dim++) {if (sqrt(2.0)*sig(i_dim)/(bb(i_dim)-aa(i_dim)) < alpha) alpha=sqrt(2.0)*sig(i_dim)/(bb(i_dim)-aa(i_dim));}
    alpha=0.8*alpha;

    double *pSources;
    double Bandwidth=alpha;
    double *pWeights;
    double *pTargets;
    double epsilon = 1e-2;

    pSources=new double[Dim*NSources];
    pWeights=new double[NSources];
    pTargets=new double[Dim*MTargets];

    for (int i=0;i<Dim*NSources;i++){
      int i_dim=i%Dim;
      int id=(int) (i-i_dim)/Dim;

      pSources[i]=0.1+alpha*(source(id,i_dim)-aa(i_dim))/(sqrt(2.)*sig(i_dim));
      //printf("pSources(%d)=%lg\n",i,pSources[i]);
      if (pSources[i]<0) {printf("Source point < 0\n"); source(id,i_dim)=0.;}
      if (pSources[i]>1)  {printf("Source point > 1\n"); source(id,i_dim)=1.;}
    }

    for (int i=0;i<NSources;i++){
      pWeights[i]=weight(i);
    }

    for (int i=0;i<Dim*MTargets;i++){
      int i_dim=i%Dim;
      int id=(int) (i-i_dim)/Dim;

      pTargets[i]=0.1+alpha*(target(id,i_dim)-aa(i_dim))/(sqrt(2.)*sig(i_dim));
      //printf("pTargets(%d)=%lg\n",i,pTargets[i]);
      if (pTargets[i]<0 || pTargets[i]>1) {printf("Heads Up:Target outside [0,1] range\n");}
    }

    int W=1;
    double * g_auto = new double[W*MTargets];

    memset( g_auto, 0, sizeof(double)*W*MTargets );
    figtree( Dim, NSources, MTargets, W, pSources, Bandwidth, pWeights, pTargets, epsilon, g_auto,
			       FIGTREE_EVAL_AUTO, FIGTREE_PARAM_NON_UNIFORM, FIGTREE_TRUNC_CLUSTER, 0 );

    for (int i=0;i<MTargets;i++){
        density(i) = g_auto[i]/NSources;
        for (int i_dim=0;i_dim<Dim;i_dim++){
            density(i) *= 1./(sig(i_dim)*sqrt(2*M_PI));
        }

    }

    //  fclose(f_gt);

    delete []pSources;
    delete []pWeights;
    delete []pTargets;
    delete []g_auto;

    return;
}


// Compute the PDF of data at the given points using
// given number of clusters for cluster specific bandwidths
void getPdf_cl(Array2D<double>& data, Array2D<double>& points,
               Array1D<double>& dens, int ncl, double sfac)
{
  //double pdf_icl;

  const int ndata = data.XSize();
  const int ndim = data.YSize();
  const int npoints=points.XSize();//=dens.XSize();

  // dimension check
  if (ndim != (int) points.YSize())
    {printf("getPdf_cl: dimension error\n"); exit(1);}

  int *bClusterIndex;
  bClusterIndex=new int[ndata];
  Array1D<double> sig(ndim,0.e0);
  Array1D<double> we(ndata,1.e0);
  Array1D<double> point(ndim,0.e0);

  Array1D<double> data_1d(ndata,0.e0);
  Array1D<double> w(ndim,1.e0);
  for(int idim=0;idim<ndim;idim++){
    for(int idata=0;idata < ndata;idata++){
      data_1d(idata)=data(idata,idim);
    }
    w(idim)=data_1d(maxIndex(data_1d))-data_1d(minIndex(data_1d))+1.e-12;
  }

  if (ncl==0)  {ncl=findNumCl(data,w,10); printf("Optimal # clusters for PDF computation: %d\n",ncl);}

  Array1D<int> numData(ncl,0);
  double bestexpl ;
  bestexpl=clust_best(data,w,ncl,numData,bClusterIndex,10);
  //clust(data,w,ncl,numData,pClusterIndex);

  for (int icl=0;icl<ncl;icl++){
    Array2D<double> data_icl(numData(icl),ndim,0.e0);
    Array1D<double> dens_icl(npoints,0.e0);

    int jc=0;
    for(int i=0;i<ndata;i++){
      if (bClusterIndex[i]==icl) {
        for (int idim=0;idim<ndim;idim++){
          data_icl(jc,idim)=data(i,idim);
        }
        jc++;
      }
    }
    get_opt_KDEbdwth(data_icl,sig);

    for(int idim=0;idim<ndim;idim++) sig(idim) *= sfac;
    //a)
    //   for (int ip=0;ip<npoints;ip++){
    //       if (ip%10000 ==0) printf("PDF evaluation progress: %d/%d\n",ip,npoints);
    //   for (int idim=0;idim<ndim;idim++){
    //  point(idim)=points(ip,idim);
    //        }
    //     pdf_icl=getPdf_d(data_icl,point,sig);
    //    dens(ip)+=pdf_icl*numData(icl)/ndata;
    //     }

    // b)
    getPdf_figtree(data_icl,points,sig, dens_icl, we);

    for (int ip=0;ip<npoints;ip++){
      dens(ip)+=dens_icl(ip)*numData(icl)/ndata;
    }
#ifdef VERB
    printf("PDF evaluation: Done.\n");
#endif

  }

  delete []bClusterIndex;

  return;

}


// Covariance matrix computations
double covariance(Array1D<double>& x1, Array1D<double>& x2,Array1D<double>& param, string covtype)
{

  int n=x1.XSize();
  int nn=x2.XSize();
  if (nn!=n ) {printf("Gproc:covariance() : Error message: covariance: matrix dimensions do not match! \n"); exit(1);}

  double cov=0.0;

  if ( covtype == "SqExp" ){
    Array2D<double> B(nn,nn,0.e0);
    for(int i=0;i<nn;i++)
      B(i,i)=1./param(i);


    Array1D<double> x12;
    x12=subtract(x1,x2);

    Array1D<double> Bx12;
    prodAlphaMatVec (B, x12, 1.0, Bx12) ;

    cov=exp(-dot(x12,Bx12));
  }

  else if ( covtype == "Exp" ){

    double absdiff=0.0;
    for(int i=0;i<nn;i++)
      absdiff+=fabs(x1(i)-x2(i))/param(i);

    cov=exp(-absdiff);


  }

  else
    throw Tantrum("covariance(): covariance type is not recognized!");


  return cov;
}

void ihsU(Array2D<double> &rndnos, const int dfac, dsfmt_t *rnstate) {

  int ns   = (int) rndnos.XSize() ;
  int ndim = (int) rndnos.YSize() ;
  double *rndnosPNT = rndnos.GetArrayPointer() ;
  ihsU(ndim, ns, rndnosPNT, dfac, rnstate) ;

}

void ihsU(const int ndim, const int ns, double *rndnos, const int dfac, dsfmt_t *rnstate) {

  int *ipos = (int *) malloc( ns*ndim*sizeof(int)) ;

  ihsP(ndim, ns, ipos, dfac, rnstate) ;
  for ( int j = 0; j < ndim; j++)
    for ( int i = 0; i < ns; i++)
      rndnos[j*ns+i] = (((double) ipos[j*ns+i])+dsfmt_genrand_open_open(rnstate))*2.0/((double) ns)-1.0 ;

  free(ipos) ;

  return ;

}

void ihsP(const int ndim, const int ns, int *x, const int  dupl, dsfmt_t *rnstate) {

  double opt = ns / pow(ns,1.0/ndim) ;
  int nsdup = ns * dupl;
  int nsdim = ns * ndim;
  vector<int> avail(nsdim);
  vector<int> point(ndim*nsdup);

  for ( int i=0; i<nsdim; i++ ) x[i] = 0 ;
  for ( int i=0; i<ndim;  i++ )
    x[i*ns+ns-1] = dsfmt_genrand_uint32(rnstate) % ns ;

  for ( int i = 0; i < ndim; i++ )
    for ( int j = 0; j < ns; j++ )
      avail[i*ns+j] = j;
  for ( int i = 0; i < ndim; i++ )
    avail[i*ns+x[i*ns+ns-1]] = ns-1;

  for ( int i=0; i<ndim*nsdup; i++ ) point[i] = 0 ;

  for ( int count=ns-2; count>=1; count--) {
    vector<int> list1(count * dupl,0);
    for (int i=0; i < ndim; i++) {
      for (int k=0; k < dupl; k++)
        for (int j=0; j < count; j++)
    list1[count*k+j]=avail[i*ns+j];
      for (int k=count*dupl-1;k>=0;k--) {
  int ptidx = 0;
        if (k>0)
    ptidx=(dsfmt_genrand_uint32(rnstate) % (k+1));
        point[i*nsdup+k] = list1[ptidx];
        list1[ptidx]     = list1[k];
      }
    }

    double minall = 1.e30;
    int    best = 0;

    for (int k=0;k<dupl*count;k++) {
      double mincan = 1.e30;
      for (int j=count; j<ns; j++) {
  double pdist = 0.0;
        for (int i=0; i < ndim; i++) pdist += pow(point[i*ns*dupl+k]-x[i*ns+j],2);
  pdist  = sqrt(pdist);
  mincan = MIN(mincan,pdist);
      }
      if (fabs(mincan-opt) < minall) {
  minall = fabs(mincan-opt);
  best   = k ;
      }
    }
    for (int i=0; i < ndim; i++) x[i*ns+count] = point[i*ns*dupl+best] ;

    for (int i=0; i < ndim; i++)
      for (int j=0; j < ns; j++)
  if (avail[i*ns+j]==x[i*ns+count])
    avail[i*ns+j]=avail[i*ns+count];

  }

  for (int i=0; i < ndim; i++) x[i*ns]=avail[i*ns];

  return ;

}

void distCorr(const Array2D<double> &spl, Array2D<double> &dCor) {

  int nspl  = spl.XSize();
  int nvars = spl.YSize();
	//cout<<nspl<<" "<<nvars<<endl;

  if (nspl > 20000)
    printf("Warning ! This might be a lengthy calculation: nspl=%d\n",nspl);

  std::vector< Array2D<double> > As;
  for (int i=0; i<nvars; i++ ) {
    Array2D<double> Amat(nspl,nspl,0.0);
    for (int i1=1; i1<nspl; i1++ )
      for (int i2=0; i2<i1; i2++ ) {
        Amat(i1,i2) = fabs(spl(i1,i)-spl(i2,i));
        Amat(i2,i1) = Amat(i1,i2);
      }
    // Compute means
    Array1D<double> Arow, Acol;
    getMean(Amat, Arow, "R");
    getMean(Amat, Acol, "C");
    double Amn = get_mean(Arow);

    // subtract/add means (linewise, columnwise, overall)
    matPvec(Amat, Arow, -1.0, "R");
    matPvec(Amat, Acol, -1.0, "C");
    addVal(Amat,Amn);
    As.push_back(Amat);
  }

  Array1D<double> dVarX(nvars,0.0);
  for (int i=0; i<nvars; i++ ) {
    Array2D<double> Atmp = dotmult(As[i],As[i]);
    dVarX(i) = sqrt(get_mean(Atmp));
    //printf("%d: %e\n",i,dVarX(i));
  }

  dCor.Resize(nvars,nvars,0.0);
  for (int i1=0; i1<nvars; i1++ )
    for (int i2=0; i2<i1; i2++ ) {
      Array2D<double> Atmp = dotmult(As[i1],As[i2]);
      dCor(i1,i2) = sqrt(get_mean(Atmp))/sqrt(dVarX(i1)*dVarX(i2));
    }

  return ;

}


