/* =====================================================================================

                      The UQ Toolkit (UQTk) version 3.1.0
                          Copyright (2020) NTESS
                        https://www.sandia.gov/UQToolkit/
                        https://github.com/sandialabs/UQTk

     Copyright 2020 National Technology & Engineering Solutions of Sandia, LLC (NTESS).
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
/// \file pdf_cl.cpp
/// \author K. Sargsyan, C. Safta 2013 -
/// \brief Command-line utility for KDE given samples

#include <cstdio>
#include <stddef.h>
#include <fstream>
#include <string>
#include <math.h>
#include <iostream>


//#include "ftndefs.h"
#include <getopt.h>
#include "Array1D.h"
#include "Array2D.h"

#include "assert.h"

#include "tools.h"
#include "arrayio.h"
#include "arraytools.h"

using namespace std;



/// default input file
#define FILE_IN    "data_in.dat"
/// default number of grid pts in each dimension
#define GRID       100
/// default number of clusters
#define N_CL       1
/// default bandwidth scale factor
#define BFAC       1.0

/// Displays information about this program
int usage(){
  printf("usage: pdf_cl [-h]  [-i <filename>] [-k <ncl>] [-g <ngrid>] [-l <limitsfile>] [-f <sigfactor>] [-x<targetfile>] \n");
  printf(" -h              : print out this help message \n");
  printf(" -i <filename>   : data file (default=%s) \n",FILE_IN);
  printf(" -k <ncl>        : define number of clusters, 0 means find the optimal (default=%d) \n",N_CL);
  printf(" -g <ngrid>      : define number of grid pts in each dimension (default=%d) \n",GRID);
  printf(" -l <limitsfile> : define the file where the grid limits are specified \n");
  printf(" -f <sigfactor>  : bandwidth scaling factor (default=%e) \n",BFAC);
  printf(" -x <targetfile> : define the file where the target points are specified \n");
  printf("================================================================================\n");
  printf("Input:: --\n");
  printf("Output:: File 'dens.dat'(makes sense to plot only for ndim=1,2)\n");
  printf("--------------------------------------------------------------------------------\n");
  printf("Comments:\n a) pdf computation is based on clustering and KDE.\n");
  printf("================================================================================\n");
  exit(0);
  return 0;
}


/// Program to compute PDF via KDE given samples
int main(int argc, char *argv[])
{

  /// Read the user input
  int ndim,nsample,ngrid,ncl;
  int c;
  double bfac;

  char *filename   = FILE_IN    ;
  char *targetfile;
  char *limsfile;

  bool gflag = false;
  bool xflag = false;
  bool lflag = false;

  ngrid = GRID;
  ncl   = N_CL;
  bfac  = BFAC;

  while ((c=getopt(argc,(char **)argv,"hi:g:l:x:k:f:"))!=-1){
    switch (c) {
      case 'h':
        usage();
        break;
      case 'i':
        filename =  optarg;
        break;
      case 'g':
        ngrid =  strtol(optarg, (char **)NULL,0);
        gflag = true;
        break;
      case 'l':
        limsfile =  optarg;
        lflag = true;
        break;
      case 'x':
        targetfile =  optarg;
        xflag = true;
        break;
      case 'k':
        ncl =  strtol(optarg, (char **)NULL,0);
        break;
      case 'f':
        bfac = strtod(optarg, (char **)NULL);
        break;
      default :
        break;
    }
  }

  /// Input checks
  if (lflag && xflag){
    printf("Please do not specify both grid limits(-l) and a filename of the target points(-x)\n");
    exit(1);
  }

  if (gflag && xflag){
    printf("Please do not specify both the number of grid points(-g) and a filename of the target points(-x)\n");
    exit(1);
  }


  Array2D<double> data;
  read_datafileVS(data,filename);
  nsample = data.XSize();
  ndim    = data.YSize();

  int totpts;
  Array2D<double> points;


  if (xflag) {
    read_datafileVS(points,targetfile);
    assert(ndim == (int)points.YSize());
  }

  /// Prepare grid
  else{

    Array1D<double> data_1d(nsample,0.e0);
    double min,max;

    Array2D<double> grid(ngrid,ndim,0.e0) ;


    Array2D<double> lims ;
    if (lflag) {
      read_datafileVS(lims,limsfile);
      assert(ndim == (int)lims.XSize());
      assert(2    == (int)lims.YSize());
    }

    for(int idim=0; idim<ndim;idim++){

      if (lflag) {
        min = lims(idim,0);
        max = lims(idim,1);
      }
      else {
        for(int is=0; is<nsample;is++){
          data_1d(is) = data(is,idim);
        }
        double minn = data_1d(minIndex(data_1d));
        double maxn = data_1d(maxIndex(data_1d));
        min = minn - (maxn-minn) * 0.1;
        max = maxn + (maxn-minn) * 0.1;
      }
      for(int igrid = 0; igrid < ngrid; igrid++){
        grid(igrid,idim) = min + (max-min) * igrid / ngrid;
      }
    }

    //totpts = (int) pow(ngrid,ndim);
    //points.Resize(totpts,ndim,0.e0);
    generate_multigrid(points,grid);

	}


  totpts = (int) points.XSize();

  /// Compute densities
  Array1D<double> dens(totpts,0.e0);
  getPdf_cl(data,points,dens,ncl,bfac);

  /// Write PDF to file
  FILE* fdens;
  fdens=fopen("dens.dat","w");
  for (int it = 0; it < totpts; it++){
    for(int idim = 0; idim < ndim; idim++){
      fprintf(fdens,"%24.16lg ",points(it,idim));
    }
    fprintf(fdens,"%24.16lg \n",dens(it));
    //if (gflag && (it+1)%ngrid==0) fprintf(fdens,"\n");
  }
  fclose(fdens);

  return (0);

}



