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
#include "sampling.hpp"

// Permutes a given array in place
void Sampling::getPerm(const int nn, const int seed, int* perm)
{

  dsfmt_gv_init_gen_rand(seed);

  int j,t;
  for (int is = 0; is <nn; is++) perm[is]=is;

  for(int is = 0; is < nn; is++){
    j=is+(int) (dsfmt_gv_genrand_urv()*(nn-is));
    t       = perm[is];
    perm[is]= perm[j];
    perm[j] = t;
  }

  return;

}

// Generate an array of uniform random variables with LHS
void Sampling::unifLHS(const int nsample, const int ndim, const int zSeed, double *rvar) {

  int *perm = new int[nsample];

  dsfmt_gv_init_gen_rand(zSeed );

  int ii=0;
  for(int id=0;id<ndim;id++){

    int seed=(int) (dsfmt_gv_genrand_urv()*INT_MAX);
    this->getPerm(nsample, seed, perm);

    for(int is=0;is<nsample;is++){
      double urv=dsfmt_gv_genrand_urv();
      rvar[ii]=(urv+perm[is])/nsample;
      ii++;
    }
  }

  delete [] perm;

  return;

}

// Generate an array of uniform random variables with LHS
void Sampling::unifLHS(const int zSeed, Array2D<double> &rvar) {
  int ndim    = rvar.YSize();
  int nsample = rvar.XSize();

  this->unifLHS(nsample,ndim,zSeed,rvar.GetArrayPointer());

  return;

}

// Generate an array of uniform random variables with LHS
void Sampling::unifLHS(const int nsample, const int ndim, dsfmt_t *rnstate, double *rvar) {

  int *perm = new int[nsample];

  int ii=0;
  for(int id=0;id<ndim;id++) {

    int seed = dsfmt_genrand_uint32(rnstate);
    this->getPerm(nsample, seed, perm);

    for(int is=0; is < nsample; is++) {
      double urv=dsfmt_genrand_open_open( rnstate );
      rvar[ii]=(urv+perm[is])/nsample;
      ii++;
    }
  }

  delete []perm;

  return;

}

// Generate an array of uniform random variables with LHS
void Sampling::unifLHS(dsfmt_t *rnstate, Array2D<double> &rvar) {

  int nsample = rvar.XSize();
  int ndim    = rvar.YSize();

  this->unifLHS(nsample, ndim, rnstate, rvar.GetArrayPointer());

  return;

}

// Generate an array of standard normal random variables with LHS
void Sampling::normLHS(const int zSeed, Array2D<double> &rvar) {

  int nsample = rvar.XSize();
  int ndim    = rvar.YSize();

  this->unifLHS(zSeed, rvar);
  for (int is = 0 ; is < nsample ; is++)
    for (int id = 0 ; id < ndim ; id++)
      rvar(is,id)=invnormcdf(rvar(is,id));

  return;

}

// Generate an array of uniform random variables with LHS
void Sampling::normLHS(const int nsample, const int ndim, const int zSeed, double *rvar) {

  this->unifLHS(nsample, ndim, zSeed, rvar);
  for (int is = 0 ; is < nsample ; is++)
    for (int id = 0 ; id < ndim ; id++)
      rvar[id*nsample+is]=invnormcdf(rvar[id*nsample+is]);

  return;

}

// Generate an array of uniform random variables with IHS
void Sampling::unifIHS(const int dfac, dsfmt_t *rnstate, Array2D<double> &rndnos) {

  int ns   = (int) rndnos.XSize() ;
  int ndim = (int) rndnos.YSize() ;
  double *rndnosPNT = rndnos.GetArrayPointer() ;
  this->unifIHS(dfac, rnstate, ndim, ns, rndnosPNT) ;

}

// Generate an array of uniform random variables with IHS
void Sampling::unifIHS(const int dfac, dsfmt_t *rnstate, const int ndim, const int ns, double *rndnos) {

  int *ipos = (int *) malloc( ns*ndim*sizeof(int)) ;

  this->getIHSperm(ndim, ns, ipos, dfac, rnstate) ;
  for ( int j = 0; j < ndim; j++)
    for ( int i = 0; i < ns; i++)
      rndnos[j*ns+i] = (((double) ipos[j*ns+i])+dsfmt_genrand_open_open(rnstate))*2.0/((double) ns)-1.0 ;

  free(ipos) ;

  return ;

}

void Sampling::getIHSperm(const int ndim, const int ns, int *x, const int  dupl, dsfmt_t *rnstate) {

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
        if ( k > 0 )
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
        if (avail[i*ns+j] == x[i*ns+count])
          avail[i*ns+j] = avail[i*ns+count];

  }

  for (int i=0; i < ndim; i++) x[i*ns]=avail[i*ns];

  return ;

}
