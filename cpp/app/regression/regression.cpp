/* =====================================================================================

                      The UQ Toolkit (UQTk) version 3.1.2
                          Copyright (2022) NTESS
                        https://www.sandia.gov/UQToolkit/
                        https://github.com/sandialabs/UQTk

     Copyright 2022 National Technology & Engineering Solutions of Sandia, LLC (NTESS).
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
/// \file regression.cpp
/// \author K. Sargsyan 2015 -
/// \brief Command-line utility for linear regression

#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <assert.h>

#include "Array1D.h"
#include "Array2D.h"
#include "PCSet.h"
#include "lreg.h"
#include <getopt.h>
#include "tools.h"
#include "arraytools.h"
#include "arrayio.h"

using namespace std;


/// default x-file
#define XFILE      "xdata.dat"
/// default y-file
#define YFILE      "ydata.dat"
/// default basis type
#define BTYPE      "PC"
/// default method
#define METH       "lsq"
/// default flag to output mean (m), mean+std (ms) or mean+std+cov (msc)
#define MSC        "m"
/// default int parameter
#define INTPAR     5
/// default string parameter
#define STRPAR     "LU"
/// default multiindex file
#define MINDEXFILE "mindex.dat"
/// default tolerance for bcs
#define ETADEFAULT 1.e-3


/******************************************************************************/
/// Displays information about this program
int usage(){
  printf("usage: regression [-h] [-x<xfile>] [-y<yfile>] [-b<basistype>] [-r<meth>] [-m<msc>] [-t<xcheckfile>] ");
  printf("[-o<intpar>] [-l<dblpar>] [-s<strpar>] [-p<mindexfile>] [-e<centerfile>] [-w<regparamfile>]\n");

  printf(" -h                 : print out this help message \n");

  printf(" -x <xfile>         : xdata filename, matrix Nxd (default=%s) \n",XFILE);
  printf(" -y <yfile>         : ydata filename, matrix Nxe (default=%s) \n",YFILE);
  printf(" -b <basistype>     : basistype - RBF, PC, PC_MI, POL, POL_MI (default=%s) \n",BTYPE);
  printf(" -r <meth>          : regression method - lsq, wbcs (default=%s) \n",METH);
  printf(" -m <msc>           : output mode - m (only mean), ms (mean and var), msc (mean and cov) (default=%s) \n",MSC);
  printf(" -t <xcheckfile>    : optional filename of x-values for validation/plotting \n");
  printf(" -o <intpar>        : integer parameter (order for PC and POL) (default=%d) \n",INTPAR);
  printf(" -l <dblpar>        : optional double  parameter (for lsq: regularization lambda, if not given then finds the best) \n");
  printf("                    :                            (for wbcs: regularization lambda, if not given then a file for weights has to be given in -w) \n");
  printf(" -s <strpar>        : string  parameter (PCtype for PC and PC_MI) (default=%s) \n",STRPAR);
  printf(" -c <etapar>        : tolerance for BCS algorithm (default=%e) \n",ETADEFAULT);
  printf(" -p <mindexfile>    : multiindex file name (relevant for PC_MI and POL_MI) (default=%s) \n",MINDEXFILE);
  printf(" -e <centerfile>    : optional file name for RBF centers, if not given then centers are at data\n");
  printf(" -w <regparamfile>  : optional file name for regularization weight vector, if not given then -l has to be given\n");
  printf("================================================================================\n");
  printf("Input:: \n");
  printf("Output:: Files 'coeff.dat', 'lambdas.dat', Sigma2.dat'(if -m is ms or msc), 'Sig.dat'(if -m is ms or msc), ycheck.dat', ycheck_var.dat'(if -m is ms or msc), 'errors.dat'(if -r is lsq), 'selected.dat', 'mindex_new.dat' (if -r is wbcs)].\n");
  printf("--------------------------------------------------------------------------------\n");
  //printf("Comments: None yet.\n");
  //printf("Complexity: Not tested yet.\n");
  //printf("Todo: \n");
  printf("================================================================================\n");
  exit(0);
  return 0;
}

/******************************************************************************/

/******************************************************************************/
/// Main program of linear regression given data
int main (int argc, char *argv[]) {


  /// Set the default values
  char* xfile = XFILE;
  char* yfile = YFILE;
  char* basistype = BTYPE;
  char* mindexfile = MINDEXFILE;
  int intpar    = INTPAR;
  char* strpar=STRPAR;
  char* msc = MSC;
  char* meth=METH;

  char* xcheckfile;
  char* centerfile;
  double dblpar;
  char* regparamfile;
  double eta = ETADEFAULT; //higher eta, fewer terms retained

  bool lflag=false;
  bool tflag=false;
  bool eflag=false;
  bool wflag=false;

  /// Read the user input
  int cc;

  while ((cc=getopt(argc,(char **)argv,"hx:y:b:r:m:t:o:l:c:s:p:e:w:"))!=-1){
    switch (cc) {
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
      meth =  optarg;
      break;
    case 'm':
      msc = optarg;
      break;
    case 't':
      xcheckfile =  optarg;
      tflag=true;
      break;
    case 'o':
      intpar =  strtol(optarg, (char **)NULL,0);
      break;
    case 'l':
      dblpar =  strtod(optarg, (char **)NULL);
      lflag=true;
      break;
    case 'c':
      eta =  strtod(optarg, (char **)NULL);
      break;
    case 's':
      strpar =  optarg;
      break;
    case 'p':
      mindexfile = optarg;
      break;
    case 'e':
      centerfile =  optarg;
      eflag=true;
      break;
    case 'w':
      regparamfile =  optarg;
      wflag=true;
      break;
    default :
      break;
    }
  }

  /// Sanity checks
  if (lflag and wflag){
    printf("regression:: can not provide both -l and -w. Exiting\n") ;
    exit(1);
  }

  if (string(meth)=="wbcs"){
    if (!lflag and !wflag){
      printf("regression:: please provide either -l or -w for wbcs method. Exiting.\n");
      exit(1);
    }
  }


  /// Print the input information on screen
  fprintf(stdout,"xfile      = %s \n",xfile);
  fprintf(stdout,"yfile      = %s \n",yfile);
  fprintf(stdout,"basistype  = %s \n",basistype);
  fprintf(stdout,"method     = %s \n",meth);
  fprintf(stdout,"intpar     = %d \n",intpar);
  if (lflag)
    fprintf(stdout,"dblpar     = %lg \n",dblpar);
  fprintf(stdout,"strpar     = %s \n",strpar);
  fprintf(stdout,"msc        = %s \n",msc);
  if (tflag)
    fprintf(stdout,"xcheckfile  = %s \n",xcheckfile);
  if (eflag)
    fprintf(stdout,"centerfile  = %s \n",centerfile);
  if (wflag)
    fprintf(stdout,"regparamfile  = %s \n",regparamfile);

  /*----------------------------------------------------------------------------*/



  /// Read data
  Array2D<double> xdata,xcheck;
  read_datafileVS(xdata,xfile);
  int nx=xdata.XSize();
  int ndim=xdata.YSize();
  Array2D<double> ydata;

  read_datafileVS(ydata,yfile);

  /// Read validation check data, if any
  if (tflag)
    read_datafileVS(xcheck,xcheckfile);
  else
    xcheck=xdata;

  /// Declare the 'parent' regression object
  Lreg* reg;

  /// Go through options, RBF, PC, PC_MI, POL, POL_MI
  if (string(basistype)=="RBF"){

    Array1D<double> a,b;
    getDomain(xdata,a,b);

    Array1D<double> widths(ndim);
    for (int i=0;i<ndim;i++)
      widths(i)=(b(i)-a(i))/sqrt(2*nx);
    //write_datafile_1d(widths,"sigmas.dat");

    Array2D<double> centers;
    if (eflag){
      read_datafileVS(centers,centerfile);
    }
    else
      centers=xdata;

    reg=new RBFreg(centers,widths);
  }

  else if (string(basistype)=="PC"){
    int order=intpar;
    reg=new PCreg(strpar,order,ndim);
  }

  else if (string(basistype)=="PC_MI"){
    //assert(!oflag);
    Array2D<int> mindex;
    read_datafileVS(mindex,mindexfile);
    reg=new PCreg(strpar,mindex);
  }

  else if (string(basistype)=="POL"){
    int order=intpar;
    reg=new PLreg(order,ndim);
  }

  else if (string(basistype)=="POL_MI"){
    //assert(!oflag);
    Array2D<int> mindex;
    read_datafileVS(mindex,mindexfile);
    reg=new PLreg(mindex);
  }
  else{
    printf("Basistype %s is not recognized, should be RBF, PC, PC_MI, POL or POL_MI. Exiting\n",basistype);
    exit(0);
  }


  int nbas=reg->GetNbas();
  cout << "Dimensionality " << reg->GetNdim() << endl;
  cout << "Number of bases " << reg->GetNbas() << endl;

  /// Initialize regression
  reg->InitRegr();
  /// Set the regression model
  reg->SetRegMode(string(msc));
  /// Setup data
  reg->SetupData(xdata,ydata);

  /// Set the regularization parameters
  Array1D<double> regweights;

  if (lflag){
      if (dblpar>=0.0)
        regweights.Resize(nbas,dblpar);
      else
        regweights.Resize(nbas,reg->LSQ_computeBestLambda());
  }
  else if (wflag){
      regweights.Resize(nbas,0.0);
      read_datafile_1d(regweights,regparamfile);
  }

  else{
      if (string(meth)=="lsq"){
        regweights=reg->LSQ_computeBestLambdas();
      }
      else if (string(meth)=="wbcs"){
        printf("For wbcs, either -l or -w has to be given. Exiting.\n");
        exit(1);
      }
      else{
        printf("Method not recognized, should be lsq or wbcs. Exiting.\n");
        exit(1);
      }
  }

  // TODO make sure all regweights are positive
  reg->SetRegWeights(regweights);
  write_datafile_1d(regweights,"lambdas.dat");



  /// Go through methods, lsq, wbcs
  Array1D<int> selected;
  if (string(meth)=="lsq"){
    reg->LSQ_BuildRegr();
    selected.Resize(nbas);
    for (int i=0;i<nbas;i++)
      selected(i)=i;
  }

  else if (string(meth)=="wbcs"){

    reg->BCS_BuildRegr(selected,eta);

  }
  else{
    printf("Method not recognized, should be lsq or wbcs. Exiting.\n");
    exit(1);
  }

  write_datafile_1d(selected,"selected.dat");
  if (string(basistype)!="RBF"){
    Array2D<int> mindex_new;
    reg->GetMindex(mindex_new);
    write_datafile(mindex_new,"mindex_new.dat");
  }


  /// Write out coefficients
  Array1D<double> coef;
  reg->GetCoef(coef);
  write_datafile_1d(coef,"coeff.dat");

  if (string(msc)!="m"){
    /// Get the data variance
    double sigma2=reg->GetSigma2();
    Array1D<double> sigma2_arr(1,sigma2);
    write_datafile_1d(sigma2_arr,"sigma2.dat");
    cout << "Sigma2 = " << sigma2 << endl;

    Array2D<double> coef_cov;
    reg->GetCoefCov(coef_cov);
    write_datafile(coef_cov,"Sig.dat");
  }

  /// Evaluate at validation points
  Array1D<double> ycheck,ycheck_var;
  Array2D<double> ycheck_cov;
  reg->EvalRegr(xcheck,ycheck,ycheck_var,ycheck_cov);
  write_datafile_1d(ycheck,"ycheck.dat");
  if (string(msc)!="m")
    write_datafile_1d(ycheck_var,"ycheck_var.dat");

  /// Print standard error measures for the lsq method
  if (string(meth)=="lsq"){
    Array1D<double> errors=reg->computeErrorMetrics(string(meth));
    cout << "LOO error : " << errors(0) << endl;
    cout << "GCV error : " << errors(1) << endl;
    write_datafile_1d(errors,"errors.dat");
  }

  delete reg;

  return 0;

}
