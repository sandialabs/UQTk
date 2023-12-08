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
/// \file pce_quad.cpp 
/// \author K. Sargsyan 2012 - 
/// \brief Command-line utility for PC construction given samples

#include "PCSet.h"
#include "tools.h"
#include "arraytools.h"
#include "arrayio.h"
#include <unistd.h>

using namespace std;



/// default PC order
#define ORD 1	
/// default bandwidth for Rosenblatt transformation
#define BWDTH 0.01      
/// default input filename
#define FILEIN "data_in.dat" 
/// default chaos type
#define CHAOS "LU" 
/// default PC parameter #1
#define ALPHA 0.0  
/// default PC parameter #2
#define BETA 1.0   


/// diagnostic: frequency of showing progress of Galerkin projection
#define DIAG_GP 1000 

/******************************************************************************/
/// Displays information about this program
int usage(){
  printf("This program to find PC coefficients of a given set of data.\n");
  printf("usage: pce_quad [-h] [-o<nord>] [-f<file_in>] [-w<bw>] [-x<which_chaos>] [-a<alpha>] [-b<beta>]\n");
  printf(" -h               : print out this help message \n");
  printf(" -o <nord>        : define the PC order (default=%d) \n",ORD);
  printf(" -f <file_in>     : define the filename with input data (default=%s) \n",FILEIN);
  printf(" -w <bw>          : define sigma for Rosenblatt, if this is non-positive it finds the rule-of-thumb automatically(default=%lg) \n",BWDTH);
  printf(" -x <which_chaos> : define the PC type (default=%s) \n",CHAOS);
  printf(" -a <alpha>       : chaos parameter #1 (default=%lg) \n",ALPHA);
  printf(" -b <beta>        : chaos parameter #2 (default=%lg) \n",BETA);
  printf("================================================================================\n");
  printf("Input  : None \n");
  printf("Output : PCcoeff.dat  -  PC coefficient file \n");
  printf("Output : mindex.dat  -  Multi-indices file \n");
  printf("--------------------------------------------------------------------------------\n");
  printf("Comments  :\n");
  printf("Complexity: linear wrt number of samples, exponential wrt dimensionality.\n");
  printf("================================================================================\n");
  exit(0);
  return 0;
}

/// Program to find PC coefficient of a random variable given a set of data samples
int main(int argc, char *argv[])
{

  /// Set the defaults and parse the input arguments
  int nord = ORD;
  char* file_in=(char *)FILEIN;
  double bw = BWDTH;
  char* which_chaos=(char *)CHAOS;
  double alpha = ALPHA;
  double beta  = BETA;

  int c;

  while ((c=getopt(argc,(char **)argv,"ho:f:w:x:a:b:"))!=-1){
     switch (c) {
     case 'h':
       usage();
       break;
     case 'o':
       nord =  strtol(optarg, (char **)NULL,0);	
       break;
     case 'f':
       file_in =  optarg;	
       break;
     case 'w':
       bw =  strtod(optarg, (char **)NULL);	
       break;
     case 'x':
       which_chaos =  optarg;
       break;
     case 'a':
       alpha =  strtod(optarg, (char **)NULL);	
       break;
     case 'b':
       beta =  strtod(optarg, (char **)NULL);	
       break;
     default :
       break;
     }
  }

  /// Print the input information
  fprintf(stdout,"---------------------------------\n") ;
  fprintf(stdout,"pce_quad() parameters : \n") ;
  fprintf(stdout,"nord = %d \n",nord);
  fprintf(stdout,"file_in = %s \n",file_in);
  fprintf(stdout,"bw = %lg \n",bw);
  fprintf(stdout,"which_chaos = %s \n",which_chaos);
  fprintf(stdout,"---------------------------------\n") ;

  /// Read the input datafile and get its size
  Array2D<double> ydata;
  read_datafileVS(ydata,file_in); 
  int nsample=ydata.XSize();
  int nvar=ydata.YSize();

  /// For the projection based PC (unlike the inference), the data dimensionality has to coincide with the stochastic dimension
  int ndim=nvar;

  /// Set the chaos
  string which_chaos_str(which_chaos);
  PCSet PCModel("NISP",nord,ndim,which_chaos_str,alpha,beta);
  int npc=PCModel.GetNumberPCTerms();

  /// Get the default quadrature points
  Array2D<double> qdpts;
  PCModel.GetQuadPoints(qdpts);
  int totquad=PCModel.GetNQuadPoints();
  write_datafile(qdpts,"quadpts.dat");

  /// Transpose the data array to prepare for invRos()
  Array2D<double> ydata_t(nvar,nsample,0.e0); 
  transpose(ydata,ydata_t);
  
  /// Array to contain inverse-Rosenblatt transformed points
  Array2D<double> invRosData(totquad,ndim,0.e0); 

  /// Frequency of showing Galerkin projection progress
  int iiout=DIAG_GP;

  /// Begin Loop over all quadrature points
  for(int it=0;it<totquad;it++){
    
    // Working arrays
    Array1D<double> quadunif(ndim,0.e0);
    Array1D<double> invRosData_1s(ndim,0.e0);
    
    /// Map quadrature points to uniform[0,1]
    for(int id=0;id<ndim;id++)
      quadunif(id)= (PCtoPC(qdpts(it,id),which_chaos_str,alpha,beta,"LU",0,0)+1.)/2.;
    
    /// Map uniform[0,1] to the distribution given by the original samples via inverse Rosenblatt
    if (bw>0)
      invRos(quadunif,ydata_t,invRosData_1s,bw);
    else
      invRos(quadunif,ydata_t,invRosData_1s);
      
    // Diagnostic output
    if ((it/iiout)*iiout==it) 
      fprintf(stdout,"invRos for Galerkin projection: %d/%d=%d%% completed\n",it,totquad,it*100/totquad);

    // Store the results
    for(int idim=0;idim<ndim;idim++)
      invRosData(it,idim)=invRosData_1s(idim);
         
  }  // End Loop over all quadrature points
  write_datafile(invRosData,"quadpts_mapped.dat");
  fprintf(stdout,"invRos for Galerkin projection: done\n");

  /// Get PC coefficients by Galerkin Projection (xxik(ip,idim)=<xi_idim * psi_ip>)
  Array2D<double> c_k(npc,nvar,0.e0); // remember, above we defined nvar=ndim 

  for(int idim=0;idim<ndim;idim++){
    Array1D<double> c_k_1d(npc,0.e0); 
    Array1D<double> invRosData_1d(totquad,0.e0);
    for(int it=0;it<totquad;it++)
      invRosData_1d(it)=invRosData(it,idim);
    PCModel.GalerkProjection(invRosData_1d,c_k_1d);
    for(int ip=0;ip<npc;ip++)
      c_k(ip,idim)=c_k_1d(ip);
  }
  
  /// Write PC coefficients to the output file
  write_datafile(c_k,"PCcoeff.dat");

  /// Write multi-indices to the output file
  Array2D<int> mindex;
  PCModel.GetMultiIndex(mindex);
  write_datafile(mindex,"mindex.dat");

  return 0;

}


