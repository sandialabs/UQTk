/* =====================================================================================

                      The UQ Toolkit (UQTk) version 3.1.4
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

     Questions? Contact the UQTk Developers at https://github.com/sandialabs/UQTk/discussions
     Sandia National Laboratories, Livermore, CA, USA
===================================================================================== */
/// \file lr_regr.cpp
/// \author Prashant Rai 2016 -
/// \brief Command-line utility for low-rank regression

#include <iostream>
#include <getopt.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <string>

#include "Array1D.h"
#include "Array2D.h"

#include "lowrank.h"
// #include "lreg.h"
// #include "lr_regr.h"
#include "PCSet.h"
#include "error_handlers.h"
#include "ftndefs.h"
#include "gen_defs.h"
#include "assert.h"
//#include "alltools.h"
#include "quad.h"
//#include "deplapack.h"


#include "arrayio.h"
#include "tools.h"
#include "arraytools.h"
#include "dsfmt_add.h"

using namespace std;


// User can overwrite these default values per each run
#define XFILE      "xdata.dat"   // default x-file
#define YFILE      "ydata.dat"   // default y-file
#define BTYPE      "PC"          // default basis type
#define METH       "lsq"         // default method
//#define MSC        "m"           // default mode of computation
#define INTPAR     4     // default int parameter
#define ITERPAR    50     // default iter parameter
#define STRPAR     "LU"    // default str parameter


/******************************************************************************/
/// Displays information about this program
int usage(){
  printf("usage: lr_regr [-h] [-x<xfile>] [-y<yfile>] [-b<basistype>] [-r<rank>] [-t<xcheckfile>] ");
  printf("[-o<intpar>] [-i<maxiter] [-l<dblpar] [-s<strpar>]\n");

  printf(" -h                 : print out this help message \n");

  printf(" -x <xfile>         : xdata filename, matrix Nxd (default=%s) \n",XFILE);
  printf(" -y <yfile>         : ydata filename, matrix Nx1 (default=%s) \n",YFILE);
  printf(" -b <basistype>     : basistype - PC, POL (default=%s) \n",BTYPE);
  printf(" -r <rank>          : integer parameter - (max rank of approximation) (default=%d) \n",INTPAR);
  printf(" -t <xcheckfile>    : optional filename of x-values for validation/plotting \n");
  printf(" -o <intpar>        : integer parameter (order for PC and POL) (default=%d) \n",INTPAR);
  printf(" -i <maxiter>       : integer parameter (number of ALS iterations) (default=%d) \n",ITERPAR);
  printf(" -l <dblpar>        : optional double  parameter (for lsq: regularization lambda) \n");
  printf(" -s <strpar>        : string  parameter (PCtype for PC) (default=%s) \n",STRPAR);
  printf("================================================================================\n");
  printf("Input:: None \n");
  printf("Output:: Files 'ycheck_1.dat, ycheck_2.dat ,..., ycheck_r.dat \n");
  printf("      :: Files 'core_1.dat, core_2.dat ,..., core_r.dat \n");
  printf("      :: Files 'space_*_*.dat \n");
  //printf("--------------------------------------------------------------------------------\n");
  //printf("Comments: None yet.\n");
  //printf("Complexity: Not tested yet.\n");
  //printf("Todo: \n");
  printf("================================================================================\n");
  exit(0);
  return 0;
}

/******************************************************************************/

/******************************************************************************/
/// \brief  Main program of inferring PC expansion of a response curve given data
int main (int argc, char *argv[]) {


  /// Set the default values
  char* xfile = XFILE;
  char* yfile = YFILE;
  char* basistype = BTYPE;
  int intpar    = INTPAR;
  int rank = INTPAR;
  int maxiter = ITERPAR;
  char* strpar=STRPAR;
  char* meth=METH;

  char* xcheckfile;
  double dblpar;

  bool lflag=false;
  bool tflag=false;
  bool vflag=false;

  /// Read the user input
  int c;

  while ((c=getopt(argc,(char **)argv,"hx:y:b:r:t:o:i:s:l:"))!=-1){
    switch (c) {
    case 'h':
      usage();
      break;
    case 'x':
      xfile =  optarg;
      break;
    case 'y':
      yfile =  optarg;
      break;
    case 'b':
      basistype = optarg;
      break;
    case 'r':
      rank =  strtol(optarg, (char **)NULL,0);
      break;
    case 't':
      xcheckfile =  optarg;
      tflag=true;
      break;
    case 'v':
      vflag=true;
      break;
    case 'o':
      intpar =  strtol(optarg, (char **)NULL,0);
      break;
    case 'i':
      maxiter =  strtol(optarg, (char **)NULL,0);
      break;
    case 'l':
      dblpar =  strtod(optarg, (char **)NULL);
      lflag=true;
      break;
    case 's':
      strpar =  optarg;
      break;
    default :
      break;
    }
  }

  /// Print the input information on screen
  fprintf(stdout,"xfile      = %s \n",xfile);
  fprintf(stdout,"yfile      = %s \n",yfile);
  fprintf(stdout,"basistype  = %s \n",basistype);
  fprintf(stdout,"rank     = %d \n",rank);
  fprintf(stdout,"intpar     = %d \n",intpar);
  fprintf(stdout,"maxiter    = %d \n",maxiter);
  if (lflag)
    fprintf(stdout,"dblpar     = %lg \n",dblpar);
  if (tflag)
    fprintf(stdout,"xcheckfile  = %s \n",xcheckfile);

  /*----------------------------------------------------------------------------*/

  /// Read xdata as training and xcheck as test input data
  Array2D<double> xdata,xcheck;
  read_datafileVS(xdata,xfile);
  int nx = xdata.XSize();
  int ndim = xdata.YSize();
 /// Read true function evaluation
  Array1D<double> ydata(nx,0.e0);
  read_datafile_1d(ydata,yfile);
  int ny = ydata.XSize();

  /// Work arrays
  Array1D<double> yapprox,ytest,error;//(nx,0.e0);


  /// Sanity checks
  if(nx != ny){
    printf("Number of input samples is not the same as output samples. Exiting.\n");
    exit(1);
  }

  /// Declare tensor object
  CanonicalTensor f;

  /// Initialize an array of functional bases
  Array1D<FunctionalBases* > testFArr(1);
  int order = intpar;


  /// Set type of basis (PC or PL) in each dimension
  FunctionalBases *newBase;

  if (string(basistype)=="PC"){
    Array1D<int> baseOrder(ndim,order);
    Array1D<string> pctype(ndim,strpar);
    newBase = new PCBases(pctype,baseOrder);
    testFArr(0) = newBase;
  }
  else if(string(basistype)=="POL"){
    Array1D<int> baseOrder(ndim,order);
    newBase = new PLBases(baseOrder);
    testFArr(0) = newBase;
  }
  else{
    printf("Basistype %s is not recognized, should be PC or POL. Exiting\n",basistype);
    exit(0);
  }


  /// Initialize algorithm
  CanonicalTensorALSLeastSquares s;
  s.setBases(testFArr);
  s.setIteration(maxiter);


  /// Set the regularization parameter
  if (lflag){
    double lambda = dblpar;
    s.setLambda(lambda);
  }
  string t_filename;
  /// Loop over one to rank
  for (int r=1;r<=rank;r++){
    s.setRank(r);
    f = s.solveDirect(xdata,ydata,vflag);

    std::stringstream rank_str;
    rank_str << "core_" << r <<".dat";
    t_filename = rank_str.str();
    write_datafile_1d(f.getCore(),t_filename.c_str());
    for (int dim=0; dim<ndim; dim++){
      std::stringstream rank_dim_str;
      rank_dim_str << "space_" << r << "_" << (dim+1) <<".dat";
      t_filename = rank_dim_str.str();
      write_datafile(f.getSpace(dim),t_filename.c_str());
    }
    //f.printTensor(); // to be decided on the format to print low rank tensor

    FunctionalTensor fun(f,testFArr);
    if (tflag){ // If test data file exists
      read_datafileVS(xcheck,xcheckfile);
      ytest = fun.tensorEval(xcheck);
    }
    else{ // use xdata itself as test data
      ytest = fun.tensorEval(xdata);
    }
    char filename[20];


    sprintf(filename, "ycheck_%d.txt",r);
    write_datafile_1d(ytest,filename); // print evaluation of low rank approximant at rank = r
    yapprox = fun.tensorEval(xdata);
    error = subtract(ydata,yapprox);

    printf("Norm of the error for rank %d approximation = %e\n",r,(norm(error)/norm(ydata)));
  }
  delete newBase;
}
