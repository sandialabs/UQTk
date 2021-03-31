/* =====================================================================================

                      The UQ Toolkit (UQTk) version 3.1.1
                          Copyright (2021) NTESS
                        https://www.sandia.gov/UQToolkit/
                        https://github.com/sandialabs/UQTk

     Copyright 2021 National Technology & Engineering Solutions of Sandia, LLC (NTESS).
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
#include <cstdio>
#include <string>
#include <iostream>
#include "Utils.h"
#include "Array1D.h"

void WriteToFile(Array1D<double>& data, char* filename)
{
  int nx=data.XSize();
  
  FILE* f_out;
  if(!(f_out = fopen(filename,"w"))){ 
    printf("WriteToFile: could not open file '%s'\n",filename); 
    exit(1); 
  }

 for(int ix = 0 ; ix < nx ; ix++){
   fprintf(f_out, "%lg\n", data(ix));
 }

 if(fclose(f_out)){ 
   printf("WriteToFile: could not close file '%s'\n",filename); 
   exit(1); 
 }

 printf("Data written to '%s'\n", filename);

 return; 
}

void WriteModesToFilePtr(const double tym, const double* u, const double* v, const double* w, const int n, FILE* f_dump)
{
  // Write time
  fprintf(f_dump, "%lg ", tym);

  // Write modes
  for (int ip=0; ip < n; ip++){
    fprintf(f_dump, "%lg ", u[ip]);
  }

  for (int ip=0; ip < n; ip++){
    fprintf(f_dump, "%lg ", v[ip]);
  }

  for (int ip=0; ip < n; ip++){
    fprintf(f_dump, "%lg ", w[ip]);
  }
  fprintf(f_dump, "\n");
  
  return;
}

void WriteMeanStdDevToFilePtr(const double tym, const double u0, const double v0, const double w0, 
                                                const double u_std, const double v_std, const double w_std, FILE* f_dump)
{
  // Write time
  fprintf(f_dump, "%lg ", tym);

  // Write mean and std. dev.
  fprintf(f_dump, "%lg %lg ", u0, u_std);
  fprintf(f_dump, "%lg %lg ", v0, v_std);
  fprintf(f_dump, "%lg %lg ", w0, w_std);

  // New line
  fprintf(f_dump, "\n");

  return;
}

void WriteMeanStdDevToStdOut(const int step, const double tym, const double u0, const double v0, const double w0, 
                                                       const double u_std, const double v_std, const double w_std)
{
  // write time, u, v, w (mean and standard deviation) to screen 
  cout << "Time Step: " << step << ", time: " << tym;
  cout << ", u: " << u0 << ", u_StDev: " << u_std;
  cout << ", v: " << v0 << ", v_StDev: " << v_std;
  cout << ", w: " << w0 << ", w_StDev: " << w_std;
  cout << endl;

  return;
}

