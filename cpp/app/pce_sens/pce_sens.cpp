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
/// \file pce_sens.cpp
/// \author K. Sargsyan 2014 -
/// \brief Command-line utility for Sobol sensitivity index computation given PC

#include "PCSet.h"
#include "tools.h"
#include "arrayio.h"
#include <unistd.h>

using namespace std;



/// default PC type
#define CHAOS "LU"
/// default multiindex file
#define MINDEX_FILE "mindex.dat"
/// default coefficient file
#define COEF_FILE "PCcoeff.dat"
///default alpha parameter for PC
#define ALPHA 0.0
/// default beta parameter for PC
#define BETA 1.0


/// Displays information about this program
int usage(){
  printf("This program parses the information contained in given pair of multiindex-coefficients\n");
  printf("usage: pce_sens [-h] [-m<mindex_file>] [-f<coef_file>] [-x<which_chaos>]\n");
  printf(" -h               : print out this help message \n");
  printf(" -x <which_chaos> : define the PC type (default=%s) \n",CHAOS);
  printf(" -m <mindex_file> : define multiindex filename (default=%s) \n",MINDEX_FILE);
  printf(" -f <coef_file>   : define the coefficient filename (default=%s) \n",COEF_FILE);
  printf(" -a <alpha>       : define the alpha parameter of the quadrature (default=%lg) \n",ALPHA);
  printf(" -b <beta>        : define the beta parameter of the quadrature (default=%lg) \n",BETA);
  printf("================================================================================\n");
  printf("Input  : None \n");
  printf("Output : sp_mindex.X.dat - sparse format of multiindices, has 2*X columns, \n");
  printf("                         - e.g. for X=3, each row has a form [i j k a_i a_j a_k] \n");
  printf("                         - and corresponds to a PC term with order a_i in the i-th dimension etc. \n");
  printf("       : varfrac.dat     - a column file of variance fractions corresponding to each PC term \n");
  printf("                         - the size is Nx1, where N is the number of PC terms \n");
  printf("       : mainsens.dat    - a columm file of main sensitivities. The size is dx1, where d is the dimensionality\n");
  printf("       : totsens.dat     - a columm file of total sensitivities. The size is dx1\n");
  printf("       : jointsens.dat   - a matrix file of joint sensitivities. The size is dxd\n");
  printf("--------------------------------------------------------------------------------\n");
  printf("Comments  : Sparse format is useful in high-dimensional problems \n");
  printf("Complexity: Linear in the number of PC terms \n");
  printf("================================================================================\n");

  exit(0);
  return 0;
}


/// Main program: parses the information contained in given multiindices and corresponding coefficients
int main(int argc, char *argv[])
{

  /// Set the defaults and parse the input arguments
  int c;

  char* mindex_file=(char *)MINDEX_FILE;
  char* coef_file=(char *)COEF_FILE;
  char* which_chaos=(char *)CHAOS;
  double alpha=ALPHA;
  double beta =BETA;

  bool aflag=false;
  bool bflag=false;

  while ((c=getopt(argc,(char **)argv,"hm:f:x:a:b:"))!=-1){
     switch (c) {
     case 'h':
       usage();
       break;
     case 'm':
       mindex_file =  optarg;
       break;
     case 'f':
       coef_file =  optarg;
       break;
     case 'x':
       which_chaos =  optarg;
       break;
    case 'a':
      aflag=true;
      alpha = strtod(optarg, (char **)NULL);
      break;
    case 'b':
      bflag=true;
      beta = strtod(optarg, (char **)NULL);
      break;
     default :
       break;
     }
  }

  /// Print out input information
  fprintf(stdout,"---------------------------------\n") ;
  fprintf(stdout,"pce_sens() parameters : \n") ;
  fprintf(stdout,"mindex_file = %s \n",mindex_file);
  fprintf(stdout,"coef_file = %s \n",coef_file);
  fprintf(stdout,"which_chaos = %s \n",which_chaos);
  if ( aflag )
    fprintf(stdout,"alpha     = %lg \n",alpha);
  if ( bflag )
    fprintf(stdout,"beta     = %lg \n",beta);
  fprintf(stdout,"---------------------------------\n") ;


  /// Read the multiindex and coefficients' files
  Array2D<int> mindex;
  read_datafileVS(mindex,mindex_file);
  int npc=mindex.XSize();
  int ndim=mindex.YSize();
  Array1D<double> coef(npc,0.e0);
  read_datafile_1d(coef,coef_file);


  /// Declare PC in NISP formulation with no quadrature
  string which_chaos_str(which_chaos);
  PCSet PCModel("NISPnoq",mindex,which_chaos_str,alpha, beta);

  /// Encode the multiindex in a sparse format and print to files
  Array1D< Array2D<int> > sp_mindex;
  PCModel.EncodeMindex(sp_mindex);

  char filename[25];
  for (int i=1;i<(int)sp_mindex.XSize();i++){
    int nn=sprintf(filename,"sp_mindex.%d.dat",i);
    write_datafile(sp_mindex(i),filename);
  }


  /// Compute mean and variance of PC and variance fractions for each term
  Array1D<double> varfrac;
  double mean, var;

  mean=PCModel.ComputeMean(coef);
  var=PCModel.ComputeVarFrac(coef,varfrac);
  cout << "Mean = " << mean << endl;
  cout << "Var  = " << var << endl;
  write_datafile_1d(varfrac,"varfrac.dat");

  /// Compute main sensitivities
  Array1D<double> mainsens;
  PCModel.ComputeMainSens(coef,mainsens);
  write_datafile_1d(mainsens,"mainsens.dat");

  /// Compute total sensitivities
  Array1D<double> totsens;
  PCModel.ComputeTotSens(coef,totsens);
  write_datafile_1d(totsens,"totsens.dat");

  /// Compute joint sensitivities
  Array2D<double> jointsens;
  PCModel.ComputeJointSens(coef,jointsens);
  for(int id=0;id<ndim;id++)
    jointsens(id,id)=mainsens(id);
  write_datafile(jointsens,"jointsens.dat");

  return 0;
}


