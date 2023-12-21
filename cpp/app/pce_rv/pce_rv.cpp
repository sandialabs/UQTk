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
/// \file pce_rv.cpp 
/// \author K. Sargsyan 2014 - 
/// \brief Command-line utility for PC-related random variable generation

#include <math.h>
#include <unistd.h>

#include "PCSet.h"
#include "tools.h"
#include "arrayio.h"

using namespace std;



/// default r.v. type
#define RVTYPE "PCvar"  
/// default r.v. dimensionality
#define DIM 1	          
/// default number of samples
#define SAMPLE 1000	  
/// default PC dimensionality
#define PCDIM 1           
/// default PC order
#define PCORD 3         
/// default parameter file
#define PARAMFILE "pccf.dat" 
/// default first parameter of PC, if needed
#define AA 0.0 
/// default second parameter of PC, if needed
#define BB 1.0 
/// default multiindex file
#define MIFILE "mi.dat"
/// default seed 
#define SEED 1 
/// default string parameter
#define STRPARAM "LEG"    



/// \brief Displays information about this program
/// \todo  Add more detailed information on options. E.g. what are the different options for type of random variable?
///        When does the order need to be specified?
int usage(){
  printf("usage: pce_rv [-h]  [-w<type>] [-d<ndim>] [-n<nsample>] [-p<pcdim>] [-o<pcord>] [-f<paramfile>] [-a<a>] [-b<b>] [-s<seed>] [-x<xstr>] [-m<mstr>]\n");
  printf(" -h                : print out this help message \n");
  printf(" -w <type>         : define the type of r.v.(default=%s) \n",RVTYPE);
  printf(" -d <ndim>         : define the r.v. dimensionality (default=%d) \n",DIM);
  printf(" -n <nsample>      : define number of samples (default=%d) \n",SAMPLE);
  printf(" -p <pcdim>        : define the PC dimensionality (default=%d) \n",PCDIM);
  printf(" -o <pcord>        : define the PC order (default=%d) \n",PCORD);
  printf(" -f <paramfile>    : define a parameter file (default=%s) \n",PARAMFILE);
  printf(" -a <a>            : define the double parameter #1 (default=%lg) \n",AA);
  printf(" -b <b>            : define the double parameter #2 (default=%lg) \n",BB);
  printf(" -s <seed>         : define the seed (default=%d) \n",SEED);
  printf(" -x <xstr>         : define string a parameter(default=%s) \n",STRPARAM);
    printf(" -m <mstr>         : define multiindex file name(default=%s) \n",MIFILE);
    
  printf("================================================================================\n");
  printf("Input  : None \n");
  printf("Output : rvar.dat  -  random samples with format nsample x ndim\n");
  printf("================================================================================\n");
  exit(0);
  return 0;
}


/// Main program: generates PC-related random variables
int main(int argc, char *argv[])
{
  
  /// Set the defaults and parse input arguments
  int c;

  int ndim=DIM;
  int nsample=SAMPLE;
  int seed=SEED;
  int pcdim=PCDIM;
  int pcord=PCORD;
  char *type       = (char *)RVTYPE;;
  char *paramfile  = (char *)PARAMFILE;
  char *xstr       = (char *)STRPARAM;
    char *mstr       = (char *)MIFILE;

  double a = AA;
  double b = BB;
  
  bool oflag=false;
  bool fflag=false;
  bool aflag=false;
  bool bflag=false;
    bool mflag=false;

  while ((c=getopt(argc,(char **)argv,"hw:d:n:p:o:f:a:b:x:m:s:"))!=-1){
     switch (c) {
     case 'h':
       usage();
       break;
     case 'w':
       type =  optarg;	
       break;
     case 'd':
       ndim =  strtol(optarg, (char **)NULL,0);	
       break;
     case 'n':
       nsample =  strtol(optarg, (char **)NULL,0);	
       break;
     case 'p':
       pcdim =  strtol(optarg, (char **)NULL,0);	
       break;
     case 'o':
       oflag=true;
       pcord =  strtol(optarg, (char **)NULL,0);	
       break;
     case 'f':
       fflag=true;
       paramfile =  optarg;	
       break;
     case 'a':
       aflag=true;
       a = strtod(optarg, (char **)NULL);	
       break;
     case 'b':
       bflag=true;
       b = strtod(optarg, (char **)NULL);	
       break;
     case 'x':
       xstr = optarg;	
       break;
    case 'm':
        mflag=true;
        mstr = optarg;
        break;
     case 's':
       seed =  strtol(optarg, (char **)NULL,0);	
       break;
     default :
       break;
     }
  }

  /// Print out input information
  fprintf(stdout,"---------------------------------\n") ;
  fprintf(stdout,"pce_rv() parameters : \n") ;
  fprintf(stdout,"type      = %s \n",type);
  fprintf(stdout,"ndim      = %d \n",ndim);
  fprintf(stdout,"nsample      = %d \n",nsample);
  fprintf(stdout,"pcdim        = %d \n",pcdim);
  if (oflag)
    fprintf(stdout,"pcord        = %d \n",pcord);
  if (fflag)
    fprintf(stdout,"paramfile  = %s \n",paramfile);  
  if (aflag)
    fprintf(stdout,"a     = %lg \n",a);
  if (bflag)
    fprintf(stdout,"b     = %lg \n",b);
  fprintf(stdout,"xstr     = %s \n",xstr);  
    if (mflag)
        fprintf(stdout,"mstr     = %s \n",mstr);
    
    // Declare the samples container
    Array2D<double> rvar;
    // The r.v. type cast as a string
    string type_str(type);
    // PC type as a string
    string pcType(xstr);
    // Multiindex file name
    string miFile(mstr);
    

 
////////////////////////////////////////////////////////////////////////

  /// Go through options for R.V. generation
  if (type_str=="PC")
    {
      // Resize the r.v. array
      rvar.Resize(nsample,ndim,0.e0);
      // Make sure the requested dimnesionality is equal to 1
      // In the future, shoud generalize this
      assert(ndim==1);

      // Declare the PC object 
      PCSet currPCModel("NISPnoq",pcord,pcdim,pcType,a, b);

      // The number of PC basis terms
      int npc=currPCModel.GetNumberPCTerms();
     
      // Read the coefficient file
      // Note: it has to have the correct number of entries
      Array1D<double> c_k(npc,0.e0); 
      read_datafile_1d(c_k,paramfile);
     
      // A temporary container for the samples
      Array1D<double> rvar_1d(nsample,0.e0);
      // Draw sample set
      currPCModel.SeedBasisRandNumGen(seed);
      currPCModel.DrawSampleSet(c_k,rvar_1d);

      // Cast the 1d array as a 2d one
      for(int is=0;is<nsample;is++)
	       rvar(is,0)=rvar_1d(is);
     
    }
    ////////////////////////////////////////////////////////////////////////

  else if (type_str=="PCmi")
  {
      // Resize the r.v. array
      rvar.Resize(nsample,ndim,0.e0);
      // Make sure the requested dimnesionality is equal to 1
      // In the future, shoud generalize this
      assert(ndim==1);
      
      // Read the multiindex given the file
      Array2D<int> mindex;
      read_datafileVS(mindex,mstr);
      assert(pcdim==mindex.YSize());

      // Declare the PC object
      PCSet currPCModel("NISPnoq",mindex,pcType,a, b);
      
      // The number of PC basis terms
      int npc=currPCModel.GetNumberPCTerms();
      
      // Read the coefficient file
      // Note: it has to have the correct number of entries
      Array1D<double> c_k(npc,0.e0);
      read_datafile_1d(c_k,paramfile);
      
      // A temporary container for the samples
      Array1D<double> rvar_1d(nsample,0.e0);
      // Draw sample set
      currPCModel.SeedBasisRandNumGen(seed);
      currPCModel.DrawSampleSet(c_k,rvar_1d);
      
      // Cast the 1d array as a 2d one
      for(int is=0;is<nsample;is++)
          rvar(is,0)=rvar_1d(is);
      
  }
  ////////////////////////////////////////////////////////////////////////

  else if (type_str=="PCvar")
    {
      // Resize the r.v. array
      rvar.Resize(nsample,pcdim,0.e0);

      // Sanity check for dimensionalities to match
      assert(ndim==pcdim);

      // Dummy variable as order is actually not needed here
      int dummy_order = 1;
      
      // Declare the PC object
      PCSet currPCModel("NISPnoq",dummy_order,pcdim,pcType,a, b);
      
      // Draw sample set
      currPCModel.SeedBasisRandNumGen(seed);
      currPCModel.DrawSampleVar(rvar);
      
    }
  
  ////////////////////////////////////////////////////////////////////////
  
  
  else  {
    printf("pce_rv: unknown random variable type.\n"); 
    exit(1);
  }
  
  
  ////////////////////////////////////////////////////////////////////////
 
  
  /// Write out to a file
  write_datafile(rvar,"rvar.dat");
  
 


  return 0;
}

