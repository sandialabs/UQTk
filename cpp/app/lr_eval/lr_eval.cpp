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
/// \file lr_eval.cpp
/// \author Prashant Rai 2018 -
/// \brief Command-line utility for low-rank evaluation

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
#include "PCSet.h"
#include "error_handlers.h"
#include "ftndefs.h"
#include "gen_defs.h"
#include "assert.h"
#include "quad.h"


#include "arrayio.h"
#include "tools.h"
#include "arraytools.h"
#include "dsfmt_add.h"

using namespace std;


// User can overwrite these default values per each run
#define BTYPE      "PC"          // default basis type
#define METH       "lsq"         // default method
//#define MSC        "m"           // default mode of computation
#define STRPAR     "LU"    // default str parameter
#define INTPAR     1     // default int parameter


/******************************************************************************/
/// Displays information about this program
int usage(){
  printf("usage : lr_eval [-h] [-b<basistype>] [-r<rank>] [-t<xdatafile>] ");
  printf(" -h                 : print out this help message \n");
  printf(" -b <basistype>     : basistype - PC, POL (default=%s) \n",BTYPE);
  printf(" -r <rank>          : integer parameter - (max rank of approximation) (default=%d) \n",INTPAR);
  printf(" -t <xdatafile>    : optional filename of x-values for validation/plotting \n");

  printf("================================================================================\n");
  printf("Input:: core_<rank>.dat, space_<rank>_*.dat \n");
  printf("Output:: Files 'ydata_(rank).dat\n");
  printf("================================================================================\n");
  exit(0);
  return 0;
}

/******************************************************************************/

/******************************************************************************/
/// Main program for canonical low-rank tensor evaluation
int main (int argc, char *argv[]) {


  /// Set the default values
  char* basistype = BTYPE;
  int rank = INTPAR;
  char* strpar=STRPAR;


  char* xdatafile;
  double dblpar;

  /// Read the user input
  int c;

  while ((c=getopt(argc,(char **)argv,"hb:r:t:"))!=-1){
    switch (c) {
    case 'h':
      usage();
      break;
    case 'b':
      basistype = optarg;
      break;
    case 'r':
      rank =  strtol(optarg, (char **)NULL,0);
      break;
    case 't':
      xdatafile =  optarg;
      break;
    default :
      break;
    }
  }

  /// Print the input information on screen

  fprintf(stdout,"basistype  = %s \n",basistype);
  fprintf(stdout,"rank     = %d \n",rank);
  fprintf(stdout,"xdatafile  = %s \n",xdatafile);

  /*----------------------------------------------------------------------------*/

  /// Read xcheck as test input data
  Array2D<double> xcheck,space;
  read_datafileVS(xcheck,xdatafile);
  int nx = xcheck.XSize();
  int ndim = xcheck.YSize();

  /// Declare arrays to store cores and spaces for functional tensor
  Array1D<double> core,ytest;
  Array1D< Array2D<double> > spaces;
  Array1D<int> baseOrder(ndim,0);

  /// Read core data file

  char t_filename[20];
  sprintf(t_filename, "core_%d.dat",rank);

  read_datafileVS(core,t_filename);

  // Read space data file

  string tt_filename,space_file = "space";
  ostringstream convert;
  for (int dim=0; dim<ndim; dim++){
    convert << rank;
    tt_filename = space_file + "_" + convert.str();
    convert.str("");
    convert.clear();
    convert << dim+1;
    tt_filename = tt_filename + "_" + convert.str() + ".dat";
    convert.str("");
    convert.clear();
    read_datafileVS(space,tt_filename.c_str());
    tt_filename.clear();
    baseOrder(dim) = space.XSize()-1;
    spaces.PushBack(space);
  }

  // Declare tensor object

  CanonicalTensor f(core,spaces);

  /// Initialize an array of functional bases
  Array1D<FunctionalBases* > basisArray(1);


  /// Set type of basis (PC or PL) in each dimension
  FunctionalBases *newBase;

  if (string(basistype)=="PC"){
    Array1D<string> pctype(ndim,strpar);
    newBase = new PCBases(pctype,baseOrder);
    basisArray(0) = newBase;
  }
  else if(string(basistype)=="POL"){
    newBase = new PLBases(baseOrder);
    basisArray(0) = newBase;
  }
  else{
    printf("Basistype %s is not recognized, should be PC or POL. Exiting\n",basistype);
    exit(0);
  }

  /// Set functional tensor

  FunctionalTensor fun(f,basisArray);
  read_datafileVS(xcheck,xdatafile);

  /// Evaluate tensor
  ytest = fun.tensorEval(xcheck);
  char filename[10];
  sprintf(filename, "ydata_%d.dat",rank);
  write_datafile_1d(ytest,filename); // print evaluation of low rank approximant at rank = r
 }
