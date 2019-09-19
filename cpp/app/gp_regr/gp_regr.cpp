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
/// \file gp_regr.cpp 
/// \author K. Sargsyan 2015 - 
/// \brief Command-line utility for Gaussian Process regression

#include <iostream>
#include <getopt.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <string>

#include "Array1D.h"
#include "Array2D.h"

#include "PCSet.h"
#include "error_handlers.h"
#include "ftndefs.h"
#include "gen_defs.h"
#include "assert.h"
#include "quad.h"
#include "gproc.h"

#include "arrayio.h"
#include "tools.h"
#include "arraytools.h"
#include "dsfmt_add.h"

using namespace std;



/// default x-file
#define XFILE      "xdata.dat"  
/// default y-file
#define YFILE      "ydata.dat"   
/// default flag to output mean (m), mean+std (ms) or mean+std+cov (msc)
#define MSC "ms" 
/// default PC order 
#define ORD 3    

/******************************************************************************/

/// Displays information about this program
int usage(){
  printf("usage: gp_regr [-h]  [-x<xfile>] [-y<yfile>] [-m<msc>] [-t<xcheckfile>] [-o<nord>]  [-l<corlength>] [-w<corlength_file>] [-s<datavar_file>]\n");
  printf(" -h        : print out this help message \n");
  printf(" -x <xfile>         : xdata filename, matrix Nxd (default=%s) \n",XFILE);
  printf(" -y <yfile>         : ydata filename, matrix Nxe (default=%s) \n",YFILE);
  printf(" -m <msc>  : flag to determine whether only mean is needed or not (default=%s) \n",MSC);
  printf(" -t <xcheckfile>    : optional filename of x-values for validation/plotting \n");
  printf(" -o <nord> : define the PC order (default=%d) \n",ORD);
  printf(" -l <corlength>        : optional correlation length (isotropic)\n");
  printf(" -w <corlength_file>   : optional file name for correlation lengths, if not given and -l is not given, then finds the best.\n");
  printf(" -s <datavar_file>     : optional file name for data variance. If not given, set to a small nugget 1.e-6\n");
  printf("================================================================================\n");
  printf("Input:: \n");
  printf("Output:: ycheck.dat, ycheck_std.dat (if ms), cov.dat, xycov.dat, sttmat.dat (if msc)\n");

  printf("--------------------------------------------------------------------------------\n");
  printf("Comments: None yet.\n");
  printf("Complexity: Not tested yet.\n");
  printf("Todo: -Clean up the matrix manipulations\n");
  printf("Todo: -Investigate addition of a stochastic dimension with roughness parameter=infty\n");
  printf("================================================================================\n");
  exit(0);
  return 0;
}


/******************************************************************************/
///  Main program of building Gaussian Process response surface
int main (int argc, char *argv[]) {

  /// Set the default values
  char* xfile = XFILE;
  char* yfile = YFILE;
  int nord=ORD;
  char* msc=MSC;


  char* xcheckfile;
  double corlength;
  char* corlength_file;
  char* datavar_file;

  bool tflag=false;
  bool lflag=false;
  bool wflag=false;
  bool sflag=false;


  /// Read the user input
  int c;

  while ((c=getopt(argc,(char **)argv,"hx:y:m:t:o:l:w:s:"))!=-1){
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
     case 'm':
      msc =  optarg;
      break;
    case 't':
      xcheckfile =  optarg;
      tflag=true;
      break;
    case 'o':
      nord =  strtol(optarg, (char **)NULL,0);	
      break;
    case 'l':
      corlength =  strtod(optarg, (char **)NULL);
      lflag=true;
      break;
    case 'w':
      corlength_file =  optarg;
      wflag=true;
      break;
    case 's':
      datavar_file =  optarg;
      sflag=true;
      break;
     }
  }
  
  /// Sanity checks
  if (lflag and wflag){
    printf("gp_regr():: can not provide both -l and -w. Exiting\n") ;
    exit(1);
  }

  /// Print the input information on screen
  fprintf(stdout,"\ngp_regr() : parameters ================================= \n");
  fprintf(stdout,"xfile  = %s \n",xfile);
  fprintf(stdout,"yfile  = %s \n",yfile);
  fprintf(stdout,"msc    = %s \n",msc);  
  fprintf(stdout,"nord   = %d \n",nord);
  if (lflag)
    fprintf(stdout,"corlength     = %lg \n",corlength);
  if (wflag)
    fprintf(stdout,"corlength_file  = %s \n",corlength_file);
  if (sflag)
    fprintf(stdout,"datavar_file  = %s \n",datavar_file);
  /*----------------------------------------------------------------------------*/ 

  string msc_str(msc);

  /// Read data  
  Array2D<double> xdata,xcheck;
  read_datafileVS(xdata,xfile);
  int nx=xdata.XSize();
  int ndim=xdata.YSize();
  Array1D<double> ydata(nx,0.e0);
    
  read_datafile_1d(ydata,yfile);
  
  /// Set or read data variance
  Array1D<double> datavar(nx,1.e-6);
  if (sflag)
    read_datafile_1d(datavar,datavar_file);




  /// Read validation check data, if any
  if (tflag)
    read_datafileVS(xcheck,xcheckfile);
  else
    xcheck=xdata;
  int ncheck=xcheck.XSize();


  printf("gp_regr() : Number of training points : %d\n",nx);
  printf("gp_regr() : Dimensionality            : %d\n",ndim);
  printf("gp_regr() : Number of test points     : %d\n",ncheck);

  /// Set the correlation parameters
  Array1D<double> corlengths;

  if (lflag){
      assert (corlength>=0.0);
      corlengths.Resize(ndim,corlength);
  }
  else if (wflag){
      corlengths.Resize(ndim,0.0);
      read_datafile_1d(corlengths,corlength_file);
  }
  else{
    corlengths.Resize(ndim,0.0);
  }

  
  /// Set the PC trend
  PCSet PCModel("NISPnoq",nord,ndim,"LEG",0.0,1.0);

  /// Initialize a GP object
  Gproc gpr("SqExp",&PCModel,corlengths); 
  gpr.SetupPrior();
  gpr.SetupData(xdata,ydata,datavar);
  if (!lflag and !wflag){
    gpr.findBestCorrParam();
    /// Print out the roughness param
    Array1D<double> bestparam;
    gpr.getParam(bestparam);
    write_datafile_1d(bestparam,"best_corlengths.dat");
  }

  int npc=gpr.getNPC();
  int al=gpr.getAl();

  /// Sanity check to ensure the regression is well-defined
  if (nx+2*al<=npc+2){
    printf("gp_regr() : Error Message: [Number of input points + Prior constraint on sigma <= Number of PC terms + 2] Student-t process will have infinite variance.\n"); 
    exit(1);
  }

   /// Build the GP
  gpr.BuildGP();

  /// Evaluate the GP (actually, a Student-t process, see the UQTk Manual)
  Array1D<double> mst;
  gpr.EvalGP(xcheck,msc_str,mst);

  /// Write the mean
  write_datafile_1d(mst,"ycheck.dat");
  //write_datafile_1d(bhat,"PCcoeff.dat");

  /// If asked, compute and write standard deviation and covariance of the Student-t process
  Array2D<double> cov;
  Array1D<double> var;
  if (msc_str != "m"){
    gpr.getVar(var);
    //Array1D<double> std;
    //for(int it=0;it<ncheck;it++)
    //  std.PushBack(sqrt(var(it)));
    write_datafile_1d(var,"ycheck_var.dat");
  }

  if (msc_str == "msc"){
    gpr.getCov(cov);
    Array1D<double> sttmat;
    Array2D<double> xycov;
    
    gpr.getXYCov(xcheck,xycov);
    gpr.getSttPars(sttmat);
    
    write_datafile(cov,"cov.dat");
    write_datafile(xycov,"xycov.dat");
    write_datafile_1d(sttmat,"sttmat.dat");
  }
  
  /// Print out output information
  fprintf(stdout,"gp_regr() : mean values                       : ycheck.dat\n");
  if (msc_str != "m")  fprintf(stdout,"gp_regr() : standard deviation                : ycheck_std.dat\n");
  if (msc_str == "msc"){
    fprintf(stdout,"gp_regr() : covariance matrix                 : cov.dat\n");
    fprintf(stdout,"gp_regr() : covariance matrix to plot         : xycov.dat\n");
    fprintf(stdout,"gp_regr() : scale matrix values and d.o.f. for student-t : sttmat.dat\n");
  }
  fprintf(stdout,"gp_regr() : done ========================================\n");

  return 0;

}
