/* =====================================================================================

                      The UQ Toolkit (UQTk) version 3.1.3
                          Copyright (2023) NTESS
                        https://www.sandia.gov/UQToolkit/
                        https://github.com/sandialabs/UQTk

     Copyright 2023 National Technology & Engineering Solutions of Sandia, LLC (NTESS).
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
/// \file generate_quad.cpp
/// \author K. Sargsyan 2013 -
/// \brief Command-line utility to generate quadrature points

#include <unistd.h>
#include "quad.h"
#include "tools.h"
#include "arrayio.h"

using namespace std;

/// default value of parameter (level for sparse quadrature, or number of grid points for full quadrature)
#define PARAM 3
/// default data dimensionality
#define DIM 2
/// default sparseness type (full or sparse)
#define FSTYPE "sparse"
/// default quadrature type
#define QUADTYPE "CC"
/// default alpha parameter for chaos
#define ALPHA 0.0
/// default beta parameter for chaos
#define BETA 1.0
/// default domain file
#define DOMAIN_FILE "param_domain.dat"
/// default verbosity
#define VERBOSITY 1

/******************************************************************************/
/// \brief Displays information about this program
int usage(){
  printf("usage: generate_quad [-h] [-r] [-d<ndim>] [-g<quadType>] [-x<fsType>] [-p<param>] [-a<alpha>] [-b<beta>] [-s<domain_file>] [-v <verb>]\n");
  printf(" -h                 : print out this help message \n");
  printf(" -r                 :  use if building the next quadrature level on top of existing rule\n");
  printf(" -d <ndim>          : define the data dimensionality (default=%d) \n",DIM);
  printf(" -g <quadType>      : define the quad type, implemented 'CC','CCO','NC','NCO','LU','HG','JB','LG','SW','pdf'.  (default=%s) \n",QUADTYPE);
  printf(" -x <fsType>        : define 'full' or 'sparse'  (default=%s) \n",FSTYPE);
  printf(" -p <param>         : define the level or nquad parameter(default=%d) \n",PARAM);
  printf(" -a <alpha>         : define the alpha parameter of the quadrature (default=%lg) \n",ALPHA);
  printf(" -b <beta>          : define the beta parameter of the quadrature (default=%lg) \n",BETA);
  printf(" -s <domain_file>   : define the domain file for compact-support quadratures (default=%s) \n",DOMAIN_FILE);
  printf(" -v <verb>          : define verbosity 0-no output/1-output info (default=%d) \n",VERBOSITY);
  printf("================================================================================\n");
  printf("Input  :  If -r flagged, files qdpts.dat, wghts.dat, indices.dat required as quadrature will be built on top of them\n");
  printf("Output : qdpts.dat, wghts.dat, indices.dat - quadrature points, weights, and indices w.r.t. default quadrature domain\n");
  printf("         xqdpts.dat, xwghts.dat            - quadrature points and weights w.r.t. given physical domain for compact domains,\n");
  printf("                                             if the domain is given by -s\n");
  printf("         *_new.dat                         - newly generated points/weights; if -r is not flagged these are the same as all points/wghts)\n");
  printf("--------------------------------------------------------------------------------\n");
  printf("Comments: -r flag may be activated only after a run with the SAME parameters, otherwise incremental addition does not make sense!\n");
  printf("================================================================================\n");
  exit(0);
  return 0;
}
/******************************************************************************/



///  Main program: Generates various kinds of quadrature points and weights
int main (int argc, char *argv[])
{
  /// Set the default values
  int    verb     = VERBOSITY ;
  int    ndim     = DIM ;
  char*  quadType = (char *) QUADTYPE;
  char*  fsType   = (char *) FSTYPE ;
  int    param    = PARAM ;
  double alpha    = ALPHA;
  double beta     = BETA;
  char*  domain_file = (char *) DOMAIN_FILE;

  /// Read the user input
  int c;

  bool rflag=false;
  bool aflag=false;
  bool bflag=false;
  bool sflag=false;

  while ((c=getopt(argc,(char **)argv,"hrd:g:x:p:a:b:s:v:"))!=-1){
    switch (c) {
    case 'h':
      usage();
      break;
    case 'r':
      rflag=true;
      break;
    case 'd':
      ndim = strtol(optarg, (char **)NULL,0);
      break;
    case 'g':
      quadType = optarg;
      break;
    case 'x':
      fsType = optarg;
      break;
    case 'p':
      param = strtol(optarg, (char **)NULL,0);
      break;
    case 'a':
      aflag=true;
      alpha = strtod(optarg, (char **)NULL);
      break;
    case 'b':
      bflag=true;
      beta = strtod(optarg, (char **)NULL);
      break;
    case 's':
      sflag=true;
      domain_file = optarg;
      break;
    case 'v':
      verb =  strtol(optarg, (char **)NULL,0);
      break;
    default :
      break;
    }
  }

  /// Print the input information on screen
  if ( verb > 0 ) {
    fprintf(stdout,"generate_quad() : parameters ================================= \n");
    fprintf(stdout," ndim     = %d \n",ndim);
    fprintf(stdout," quadType = %s \n",quadType);
    fprintf(stdout," fsType   = %s \n",fsType);
    fprintf(stdout," param    = %d \n",param);
    if (aflag)
      fprintf(stdout," alpha    = %lg \n",alpha);
    if (bflag)
      fprintf(stdout," beta     = %lg \n",beta);
    if (rflag)
      fprintf(stdout,"generate_quad() : building on top of existing quad points\n");
    if (sflag)
      fprintf(stdout,"generate_quad() : domain file %s is provided\n",domain_file);
  }
  /*----------------------------------------------------------------------------*/

  /// Parameter sanity checks
  if (rflag && string(fsType)=="full")
      throw Tantrum("Incremental addition makes sense only in the sparse mode!");
  if (sflag && string(quadType)!="CC"
            && string(quadType)!="CCO"
            && string(quadType)!="NC"
            && string(quadType)!="NCO"
            && string(quadType)!="LU"
            && string(quadType)!="JB")
      throw Tantrum("Input domain should be provided only for compact-support quadratures!");


  /// Declare the quadrature rule object
  Array1D<string> quadtypes(ndim,string(quadType));

  Array1D<double> alphas(ndim,alpha);
  Array1D<double> betas(ndim,beta);
  Array1D<int> params(ndim,param);

  Quad spRule(quadtypes,fsType,params,alphas, betas);
  spRule.SetVerbosity(verb);

  // Declare arrays
  Array1D<int> newPtInd;
  Array2D<double> qdpts;
  Array1D<double> wghts;

  spRule.SetRule();

  // DEBUG
  //Array1D<int> ind;
  //spRule.compressRule(ind);

  /// Extract the properties of the rule
  spRule.GetRule(qdpts,wghts);
  int nQdpts=qdpts.XSize();


  /// Write-out to files
  write_datafile(qdpts,"qdpts.dat");
  write_datafile_1d(wghts,"wghts.dat");

  /// Scale if domain is provided
  if (sflag){
    /// Set the domain
    Array1D<double> aa(ndim,-1.e0);
    Array1D<double> bb(ndim,1.e0);
    Array2D<double> aabb(ndim,2,0.e0);

    if(ifstream(domain_file)){
      read_datafile(aabb,domain_file);
      for (int i=0;i<ndim;i++){
        aa(i)=aabb(i,0);
        bb(i)=aabb(i,1);
      }
    }

    Array2D<double> xqdpts(nQdpts,ndim);
    //   Array2D<double> xqdpts_new(nNewQdpts,ndim);
    Array1D<double> xwghts(nQdpts);
    //    Array1D<double> xwghts_new(nNewQdpts);

    // Scale points according to the given domain
    for(int it=0;it<nQdpts;it++){
      double prod=1.;
      for(int id=0;id<ndim;id++){
        xqdpts(it,id)=(bb(id)+aa(id))/2.+qdpts(it,id)*(bb(id)-aa(id))/2.;
        prod*= ( (bb(id)-aa(id))/2. );
      }
      xwghts(it)=wghts(it)*prod;
    }

    // /// Scale weights according to the given domain
    // for(int iq=0;iq<nNewQdpts;iq++){
    //   double prod=1.;
    //   for(int id=0;id<ndim;id++){
    //     xqdpts_new(iq,id)=(bb(id)+aa(id))/2.+qdpts_new(iq,id)*(bb(id)-aa(id))/2.;
    //     prod*= ( (bb(id)-aa(id))/2. );
    //   }
    //   xwghts_new(iq)=wghts_new(iq)*prod;
    // }

    /// Write-out to files
    write_datafile(xqdpts,"xqdpts.dat");
    write_datafile_1d(xwghts,"xwghts.dat");
    // write_datafile(xqdpts_new,"xqdpts_new.dat");
    // write_datafile_1d(xwghts_new,"xwghts_new.dat");

  }

  if ( verb > 0 ) {
    //fprintf(stdout,"generate_quad() : generated %d new quadrature points\n",nNewQdpts);
    fprintf(stdout,"generate_quad() : total number of quadrature points: %d\n",nQdpts);
    fprintf(stdout,"generate_quad() : done ========================================\n");
  }

  return 0;
}


