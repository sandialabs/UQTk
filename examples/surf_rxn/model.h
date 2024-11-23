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
#include "PCSet.h"

struct dumpInfo
{
  int dumpInt;
  int fdumpInt;
  string dumpfile;
};

// Function declarations



/// \brief Compute RHS of monomer-dimer reaction model, given the parameters a through f, and
/// the current state u, v, w, z. Return result in dudt, dvdt, dwdt
void GetRHS(const double& a, const double& b, const double & c, const double & d, const double & e, const double & f,
            const double& u, const double& v, const double& w, const double& z,
            double& dudt, double& dvdt, double& dwdt);


/// \brief Compute RHS of monomer-dimer reaction model, given the parameters a through f, and
/// the current state u, v, w, z. Return result in dudt, dvdt, dwdt
/// The parameter b is considered to be uncertain and is represented with a PC expansion. As a result
/// u, v, w, z, dudt, dvdt, and dwdt are also PCEs.
/// aPCSet contains all parameters of the PC type used
void GetRHS(const PCSet &aPCSet, const double* a, const double* b, const double* c, 
            const double* d, const double* e, const double* f,
            const double* u, const double* v, const double* w, const double* z,
            double* dudt, double* dvdt, double* dwdt);


/// \brief Move the outputs forward in time with time step dTym, given input parameter vector
void forwardFunctionDt(double* inpParams, double dTym, double* outValues);

/// \brief Forward function from t0 to tf, with time step dTym, given input parameter vector
void forwardFunction(double* inpParams, double t0, double tf, double dTym, double* outValues,dumpInfo* outPrint);

void forwardFunction(double* inpParams, double t0, double tf, double dTym, double* outValues, dumpInfo* outPrint)
{
  bool silent=false;
  if (!outPrint)
    silent=true;
  

 // Current time
  double tym = t0;
  int step=0;
  double u = outValues[0];
  double v = outValues[1];
  double w = outValues[2];


  FILE* f_dump;

 if(!silent){
   


 // Separator in output
  cout << endl << string(70,'=') << endl;
  cout << "Integration" << endl;
  cout << string(70,'=') << endl << endl;
 // write time, u, v, w to file
  if(!(f_dump = fopen(outPrint->dumpfile.c_str(),"w"))){ 
    printf("Could not open file '%s'\n",outPrint->dumpfile.c_str()); 
    exit(1); 
  }
  fprintf(f_dump, "%lg %lg %lg %lg\n", tym, u, v, w );

  // write u, v, w to screen
  cout << "Time Step: " << step << ", time: " << tym;
  cout << ", u: " << u;
  cout << ", v: " << v;
  cout << ", w: " << w;
  cout << endl;
}

  // Forward run
  while(tym < tf) {
    forwardFunctionDt(inpParams,dTym,outValues);
    tym += dTym;
    step+=1;

    if (!silent){
      u = outValues[0];
      v = outValues[1];
      w = outValues[2];
 
 // write time, u, v, w to file
    if(step % outPrint->fdumpInt == 0){
      fprintf(f_dump, "%lg %lg %lg %lg\n", tym, u, v, w );
    }
    
    // write u, v, w to screen
    if(step % outPrint->dumpInt == 0){
      cout << "Time Step: " << step << ", time: " << tym;
      cout << ", u: " << u;
      cout << ", v: " << v;
      cout << ", w: " << w;
      cout << endl;
    }
    }

  }

  if(!silent){
  // write time, u, v, w to file if not already done so
  if(step % outPrint->fdumpInt != 0){
    fprintf(f_dump, "%lg %lg %lg %lg\n", tym, u, v, w );
  }

  // write u, v, w to screen if not already done so
  if(step % outPrint->dumpInt != 0){
    cout << "Time Step: " << step << ", time: " << tym;
    cout << ", u: " << u;
    cout << ", v: " << v;
    cout << ", w: " << w;
    cout << endl;
  }

  // Close output file
  if(fclose(f_dump)){ 
    printf("Could not close file '%s'\n",outPrint->dumpfile.c_str()); 
    exit(1); 
  }
}

 
  return;
}





void forwardFunctionDt(double* inpParams, double dTym, double* outValues)
{
  // variables to store right hand side
  double dudt = 0.e0;
  double dvdt = 0.e0;
  double dwdt = 0.e0;
  
  // parse input parameters
  const double a=inpParams[0];
  const double b=inpParams[1];
  const double c=inpParams[2];  
  const double d=inpParams[3];
  const double e=inpParams[4];
  const double f=inpParams[5];

  // Save solution at current time step
  double u_o = outValues[0];
  double v_o = outValues[1];
  double w_o = outValues[2];
  double z_o = 1.0 - u_o - v_o - w_o;

  // Integrate with 2nd order Runge Kutta
  
  // Compute right hand sides
  GetRHS(a,b,c,d,e,f,u_o,v_o,w_o,z_o,dudt,dvdt,dwdt);

  // Advance to mid-point
  double u = u_o + 0.5*dTym*dudt;
  double v = v_o + 0.5*dTym*dvdt;
  double w = w_o + 0.5*dTym*dwdt;
  double z = 1.0 - u - v - w;
  
  // Compute right hand sides
  GetRHS(a,b,c,d,e,f,u,v,w,z,dudt,dvdt,dwdt);
  
  // Advance to next time step
  u = u_o + dTym*dudt;
  v = v_o + dTym*dvdt;
  w = w_o + dTym*dwdt;
  z = 1.0 - u - v - w;
  
  // Store the output
  outValues[0]=u;
  outValues[1]=v;
  outValues[2]=w;
  
  return;
}

void GetRHS(const double& a, const double& b, const double & c, const double & d, const double & e, const double & f,
            const double& u, const double& v, const double& w, const double& z,
            double& dudt, double& dvdt, double& dwdt)
{  
  dudt = a*z - c*u - 4.e0*d*u*v;
  dvdt = 2.e0*b*z*z - 4.e0*d*u*v;
  dwdt = e*z - f*w;

  return;
}




void GetRHS(const PCSet& aPCSet, const double* a, const double* b, const double* c, const double* d, const double* e, const double* f,
            const double* u, const double* v, const double* w, const double* z,
            double* dudt, double* dvdt, double* dwdt)
{ 
  // Number of PC terms
  const int nPCTerms = aPCSet.GetNumberPCTerms();

  // Work variable
  double* dummy1 = new double[nPCTerms];
  double* dummy2 = new double[nPCTerms];
  double* dummy3 = new double[nPCTerms];

  // Build du/dt = a*z - c*u - 4.0*d*u*v
  aPCSet.Prod(z,a,dummy1);           // dummy1 = a*z
  aPCSet.Prod(u,c,dummy2);           // dummy2 = c*u
  aPCSet.SubtractInPlace(dummy1,dummy2); // dummy1 = a*z - c*u
  aPCSet.Prod(u,v,dummy2);               // dummy2 = u*v
  aPCSet.Prod(dummy2,d,dummy3);               // dummy3 = d*u*v
  aPCSet.MultiplyInPlace(dummy3,4.e0); // dummy3 = 4.0*d*u*v
  aPCSet.Subtract(dummy1,dummy3,dudt);   // dudt = a*z - c*u - 4.0*d*u*v

  // Build dv/dt = 2.0*b*z*z - 4.0*d*u*v
  aPCSet.Prod(z,z,dummy1);               // dummy1 = z*z
  aPCSet.Prod(dummy1,b,dummy2);          // dummy2 = b*z*z
  aPCSet.Multiply(dummy2,2.e0,dummy1);   // dummy1 = 2.0*b*z*z
  aPCSet.Prod(u,v,dummy2);               // dummy2 = u*v
  aPCSet.Prod(dummy2,d,dummy3);               // dummy3 = d*u*v
  aPCSet.MultiplyInPlace(dummy3,4.e0); // dummy3 = 4.0*d*u*v
  aPCSet.Subtract(dummy1,dummy3,dvdt);   // dvdt = 2.0*b*z*z - 4.0*d*u*v

  // Build dw/dt = e*z - f*w
  aPCSet.Prod(z,e,dummy1);           // dummy1 = e*z
  aPCSet.Prod(w,f,dummy2);           // dummy2 = f*w
  aPCSet.Subtract(dummy1,dummy2,dwdt);   // dwdt = e*z - f*w

  // Memory clean-up
  delete [] dummy1;
  delete [] dummy2;
  delete [] dummy3;

  return;
}
