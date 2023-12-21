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
/// \file gen_mi.cpp 
/// \author K. Sargsyan 2014 - 
/// \brief Command-line utility to generate multiindex

#include "tools.h"
#include "arrayio.h"
#include "arraytools.h"
#include <unistd.h>

using namespace std;


/// default multiindex type
#define MI_TYPE "TO" 
/// default multiindex sequence
#define MI_SEQ "NONE" 
/// default order
#define ORD 1
/// default dimensionality		
#define DIM 3 
/// default parameter filename
#define PARAM_FILE "mi_param.dat" 
/// default verbosity
#define VERBOSITY 1



/******************************************************************************/
/// Displays information about this program
int usage(){
  printf("This program to generate multiindex files given rules.\n");
  printf("usage: gen_mi [-h]  [-x<mi_type>] [-s<mi_seq>] [-p<nord>] [-q<ndim>] [-f<param_file>] [-v <verb>]\n");
  printf(" -h               : print out this help message \n");
  printf(" -x <mi_type>     : define the multiindex type (default=%s) \n",MI_TYPE);
  printf(" -s <mi_seq>      : define the multiindex sequence (default=%s) \n",MI_SEQ);
  printf(" -p <nord>        : define the first parameter (default=%d) \n",ORD);
  printf(" -q <ndim>        : define the second parameter (default=%d) \n",DIM);
  printf(" -f <param_file>  : define the parameter filename for multiindex (default=%s) \n",PARAM_FILE);
  printf(" -v <verb>        : define verboosity 0-no output/1-output info (default=%d) \n",VERBOSITY);
  printf("================================================================================\n");
  printf("Input  : None \n");
  printf("Output : File 'mindex.dat'\n");
  printf("--------------------------------------------------------------------------------\n");
  printf("================================================================================\n");
  exit(0);
  return 0;
}


/// Main program: Generates multiindex of requested type with given parameters
int main(int argc, char *argv[])
{

  
  /// Set the default values
  int nord = ORD;
  int ndim = DIM;
  int verb = VERBOSITY;
  char* param_file= (char *)PARAM_FILE;
  char* mi_type   = (char *)MI_TYPE;
  char* mi_seq    = (char *)MI_SEQ;
    
  bool pflag = false;
  bool sflag = false;
  bool qflag = false;
  bool fflag = false;

  /// Read the user input
  int c;

  while ((c=getopt(argc,(char **)argv,"hx:s:p:q:f:v:"))!=-1){
     switch (c) {
     case 'h':
       usage();
       break;
     case 'x':
       mi_type = optarg;
       break;
     case 's':
       mi_seq  = optarg;
       sflag=true;
       break;
     case 'p':
       nord =  strtol(optarg, (char **)NULL,0);	
       pflag=true;
       break;
     case 'q':
       ndim =  strtol(optarg, (char **)NULL,0);
       qflag=true;
       break;
     case 'f':
       param_file =  optarg;
       fflag=true;
       break;
     case 'v':
       verb =  strtol(optarg, (char **)NULL,0);
       break;
     default :
       break;
     }
  }

  /*----------------------------------------------------------------------------*/ 
  /// Print the input information on screen 
  if ( verb > 0 ) {
    fprintf(stdout,"mi_type    = %s \n",mi_type);
    if (sflag) fprintf(stdout,"mi_seq     = %s \n",mi_seq);
    if (qflag) fprintf(stdout,"ndim       = %d \n",ndim);
    if (pflag) fprintf(stdout,"nord       = %d \n",nord);
    if (fflag) fprintf(stdout,"param_file = %s \n",param_file);
  }
  /*----------------------------------------------------------------------------*/

  if(fflag && (pflag)){
    printf("gen_mi(): Can not specify both parameter file and order. Exiting.\n");
    exit(1);
  }
 
  // Cast multiindex type as string
  string mi_type_str(mi_type);

  int npc;
  Array2D<int> mindex;

  // Choose between TO, TP or HDMR
  
  // Total order
  if (mi_type_str == "TO") {
    if ( not sflag )
      npc=computeMultiIndex(ndim,nord,mindex);
    else
      npc=computeMultiIndex(ndim,nord,mindex,string(mi_seq));
  }

  else if (mi_type_str == "TP") {
    Array1D<int> maxorders;
    Array2D<int> maxorders2d;
    read_datafileVS(maxorders2d,param_file);
    getCol(maxorders2d, 0, maxorders);
    npc=computeMultiIndexTP(maxorders, mindex);
  }

  // HDMR ordering
  else if(mi_type_str=="HDMR"){
    Array1D<int> maxorders;
    Array2D<int> maxorders2d;
    read_datafileVS(maxorders2d,param_file);
    getCol(maxorders2d, 0, maxorders);
    npc=computeMultiIndexHDMR(ndim, maxorders, mindex);
  }
   
  else {
    printf("gen_mi():: Multiindex type %s is not recognized. \n", mi_type);
    exit(1);
  }
   
  /// Write to file mindex.dat
  write_datafile(mindex, "mindex.dat");
  if ( verb > 0 )
    cout << "Generated multiindex of size " << npc 
         << " and stored in mindex.dat" << endl;
 
  return 0;

}


