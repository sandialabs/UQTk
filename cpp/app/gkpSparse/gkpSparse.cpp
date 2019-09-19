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
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <unistd.h>
#include <math.h>
#include <assert.h>

#include "arrayio.h"
#include "gkplib.h"
#define NDIM      2
#define NLEV      2
#define VERBOSITY 1
#define PDFTYPE   "unif"

/******************************************************************************/
/// \brief Displays information about this program
int usage(){
  printf("usage: gkpSparse [-h] [-d<ndim>] [-l<lev>] [-v <verb>]\n");
  printf(" -h         : print out this help message \n");
  printf(" -d <ndim>  : define the data dimensionality (default=%d) \n",NDIM);
  printf(" -l <param> : define the level or nquad parameter(default=%d) \n",NLEV);
  printf(" -t <pdf>   : pdf type 'unif'/'norm'/'cc' (default=%s) \n",PDFTYPE);
  printf(" -v <verb>  : define verboosity 0-no output/1-output info (default=%d) \n",VERBOSITY);
  printf("================================================================================\n");
  printf("Output : qdpts.dat, wghts.dat - quadrature points, weights w.r.t. default quadrature domain\n");
  printf("================================================================================\n");
  exit(0);
  return (0);
}

int main(int argc, char *argv[]) {

  double *qpts=NULL, *w=NULL;
  int dim=NDIM, lev=NLEV, verb=VERBOSITY, nqpts;
  char *pdftype = (char *) PDFTYPE;
  bool anisFlag = false;

  /// Read the user input
  int c; 
  while ((c=getopt(argc,(char **)argv,"had:l:v:t:"))!=-1){
     switch (c) {
     case 'h':
       usage();
       break;
     case 'd':
       dim =  strtol(optarg, (char **)NULL,0);  
       break;
     case 'l':
       lev =  strtol(optarg, (char **)NULL,0);
       break;
     case 'v':
       verb =  strtol(optarg, (char **)NULL,0);
       break;
    case 't':
      pdftype = optarg;
      break;
    case 'a':
      anisFlag = true;
      break;
    default :
      break;
     }
  }

  if ( verb > 0 ) {
    fprintf(stdout,"gkpSparse : parameters ================================= \n");
    fprintf(stdout," ndim = %d \n",dim);
    fprintf(stdout," nlev = %d \n",lev);
    fprintf(stdout," pdf  = %s \n",pdftype);
    fprintf(stdout," anis = %d \n",anisFlag);
    fprintf(stdout," verb = %d \n",verb);
  }
  
  
  Array2D<int> levList;
  if (anisFlag) {
    read_datafileVS(levList,"levList.dat");
    assert(levList.XSize() == dim);
    assert(levList.YSize() == 1  );
  }

  /* Get sparse quad */
  if ( std::string(pdftype) == std::string("unif") ) {
    if (anisFlag)
      getSpgAnisQW ( getGKPunif, getOrderGKPunif, dim, levList.GetArrayPointer(), &nqpts, &qpts, &w ) ;
    else
      getSpgQW ( getGKPunif, getOrderGKPunif, dim, lev, &nqpts, &qpts, &w ) ;
  }
  else if ( std::string(pdftype) == std::string("norm") ) {
    if (anisFlag)
      getSpgAnisQW ( getGKPnorm, getOrderGKPnorm, dim, levList.GetArrayPointer(), &nqpts, &qpts, &w ) ;
    else
      getSpgQW ( getGKPnorm, getOrderGKPnorm, dim, lev, &nqpts, &qpts, &w ) ;
  }
  else if ( std::string(pdftype) == std::string("cc") ) {
    if (anisFlag)
      getSpgAnisQW ( getCC, getOrderCC, dim, levList.GetArrayPointer(), &nqpts, &qpts, &w ) ;
    else
      getSpgQW ( getCC, getOrderCC, dim, lev, &nqpts, &qpts, &w ) ;
  }
  else {
    std::cout<<"Unknown quadrature type: "<<pdftype<<std::endl;
    exit(1);
  }

  /* Output sparse grid */
  FILE *fout=fopen("qdpts.dat","w");
  for (int j = 0; j < nqpts; j++){
    for (int i = 0; i < dim; i++) fprintf(fout,"%14.8e  ",qpts[j*dim+i]);
    fprintf(fout,"\n");
  }
  fclose(fout);

  /* Normalize weights for uniform density */
  if ( std::string(pdftype) == std::string("unif") || 
       std::string(pdftype) == std::string("cc")) {
    double wnorm = pow(2.0,dim);
    for (int j = 0; j < nqpts; j++)
      w[j] /= wnorm;
  }

  /* Output weights */
  fout=fopen("wghts.dat","w");
  double wsum = 0.0;
  for (int j = 0; j < nqpts; j++) {
    fprintf(fout,"%14.8e\n",w[j]);
    wsum += w[j];
  }
  fclose(fout);

  if ( verb > 0 ) {
    std::cout<<"No. of quadrature points:  "<<nqpts<<std::endl;
    std::cout<<"Sum of quadrature weights: "<<wsum <<std::endl;
  }

  return (0);

}
