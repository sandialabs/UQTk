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

#include <math.h>
#include "XMLUtils.h"
#include "tools.h"
#include "mcmc.h"
#include "arrayio.h"
#include "arraytools.h"
#include "PCSet.h"

#include "model.h"
#include "XMLreader.h"
#include "Utils.h"



/// \brief Example that solves deterministic 3-equation model for heterogeneous surface
/// reaction involving a monomer, dimer, and inert species adsorbing
/// onto a surface out of gas phase. This model mimics some aspects
/// of CO oxidation. Any number of the parameters are considered to be uncertain.
///
/// For more details on the reaction model, see
/// [1] R. Vigil and F. Willmore, “Oscillatory dynamics in a heterogeneous surface
///     reaction: Breakdown of the mean-field approximation.,” Phys Rev E,
///     vol. 54, no. 2, pp. 1225–1231, Aug. 1996.
/// [2] A. G. Makeev, D. Maroudas, and I. G. Kevrekidis, “‘Coarse’ stability and
///     bifurcation analysis using stochastic simulators: Kinetic Monte Carlo examples,”
///     J. Chem. Phys., vol. 116, no. 23, p. 10083, 2002.
/// The equations solved are:
///     du/dt = az - cu - 4duv (coverage fraction of monomer)
///     dv/dt = 2bz^2 - 4duv   (coverage fraction of dimer)
///     dw/dt = ez - fw        (coverage fraction of inert species)
///       z   = 1 - u - v - w  (vacant fraction)



/// Main program of uncertainty propagation of the ODE model parameters via intrusive spectral projection (ISP)
int main()
{
  // Model parameters
  Array1D<double> modelparams;
  // Model parameter names
  Array1D<string> modelparamnames;
  // Auxiliary parameters: final time and time step of integration
  Array1D<double> modelauxparams;
  
  // Read the xml tree
  RefPtr<XMLElement> xmlTree=readXMLTree("surf_rxn.in.xml");
  // Read the model-specific input
  readXMLModelInput(xmlTree,modelparams, modelparamnames, modelauxparams);
  // Total nuber of input parameters
  int fulldim=modelparams.XSize();
  // Read the output preferences
  dumpInfo* outPrint=new dumpInfo;
  readXMLDumpInfo( xmlTree, &(outPrint->dumpInt), &(outPrint->fdumpInt), &(outPrint->dumpfile) );

  // Output PC order
 int order;
 // PC type
 string pcType;

 // A 2d array (each row is an array of coefficients for the corresponding uncertain input parameter)
 Array2D<double> allPCcoefs;
 // The indices of the uncertain model parameters in the list of model parameters
 Array1D<int> uncParamInd;
 // Read the UQ-specific information from the xml tree
 readXMLUncInput(xmlTree,allPCcoefs,uncParamInd , &order, &pcType);

 // Stochastic dimensionality
 int dim=uncParamInd.XSize();
 
 // Instantiate a PC object for ISP computations
 PCSet myPCSet("ISP",order,dim,pcType,0.0,1.0); 
 
 // The number of PC terms
 const int nPCTerms = myPCSet.GetNumberPCTerms();
 cout << "The number of PC terms in an expansion is " << nPCTerms << endl;

  // Print the multiindices on screen
  myPCSet.PrintMultiIndex();

 // Initial time
 double t0 = 0.0;
 // Final time
 double tf = modelauxparams(0);
 // Time step
 double dTym = modelauxparams(1);
 // Number of steps
 int nStep=(int) tf / dTym;
 
  // Initial conditions of zero coverage (based on Makeev:2002)
  Array1D<double> u(nPCTerms,0.e0);
  Array1D<double> v(nPCTerms,0.e0);
  Array1D<double> w(nPCTerms,0.e0);
  Array1D<double> z(nPCTerms,0.e0);

  // Array to hold the PC representation of the number 1
  Array1D<double> one(nPCTerms,0.e0);
  one(0)=1.0;

  // The z-species is described as z=1-u-v-w
  z=one;
  myPCSet.SubtractInPlace(z,u);
  myPCSet.SubtractInPlace(z,v);
  myPCSet.SubtractInPlace(z,w);

  // Right-hand sides
  Array1D<double> dudt(nPCTerms,0.e0);
  Array1D<double> dvdt(nPCTerms,0.e0);
  Array1D<double> dwdt(nPCTerms,0.e0);

  // Array of arrays to hold the input parameter PC representations in the output PC
  // Each element is an array of coefficients for the corresponding input parameter, whether deterministic or uncertain
  // The size of the array is the total number input parameters
  
  Array1D< Array1D<double> > modelparamPCs(fulldim);
  

  printf("\nInput parameter PC coefficients are given below\n");
  for (int i=0; i<fulldim; i++){
    printf("%s: ",modelparamnames(i).c_str());
    modelparamPCs(i).Resize(nPCTerms,0.e0);
    for (int j=0; j<nPCTerms; j++){
      modelparamPCs(i)(j)=allPCcoefs(j,i);
      printf(" %lg ",modelparamPCs(i)(j));
    }
    printf("\n");
  }
  printf("\n");


  // Initial time and time step counter
  int step=0;
  double tym=t0;

  // Work arrays for integration
  Array1D<double> u_o(nPCTerms,0.e0);
  Array1D<double> v_o(nPCTerms,0.e0);
  Array1D<double> w_o(nPCTerms,0.e0);

  Array1D<double> tmp_u(nPCTerms,0.e0);
  Array1D<double> tmp_v(nPCTerms,0.e0);
  Array1D<double> tmp_w(nPCTerms,0.e0);

  // File to write the mean and stdev, name read from xml
  FILE *f_dump,*modes_dump;
  if(!(f_dump = fopen(outPrint->dumpfile.c_str(),"w"))){ 
    printf("Could not open file '%s'\n",outPrint->dumpfile.c_str()); 
    exit(1); 
  }
  
  // File to dump the PC modes, name hardwired
  string modes_dumpfile = "solution_ISP_modes.dat";
  if(!(modes_dump = fopen(modes_dumpfile.c_str(),"w"))){ 
    printf("Could not open file '%s'\n",modes_dumpfile.c_str()); 
    exit(1); 
  }

  // write time, u, v, w (all modes) to file
  WriteModesToFilePtr(tym, u.GetArrayPointer(), v.GetArrayPointer(), w.GetArrayPointer(), nPCTerms, modes_dump);
    
  // Write out initial step
  // Get standard deviations
  double uStDv = myPCSet.StDv(u);
  double vStDv = myPCSet.StDv(v);
  double wStDv = myPCSet.StDv(w);

  // write u, v, w (mean and standard deviation) to file
  WriteMeanStdDevToFilePtr(tym, u(0), v(0), w(0), uStDv, vStDv, wStDv, f_dump);
    
  // write u, v, w (mean and standard deviation) to screen
  WriteMeanStdDevToStdOut(step, tym, u(0), v(0), w(0), uStDv, vStDv, wStDv);
    
  // Forward run
  while(tym < tf) {
    // Integrate with 2nd order Runge Kutta

    // Save solution at current time step
    myPCSet.Copy(u_o,u);
    myPCSet.Copy(v_o,v);
    myPCSet.Copy(w_o,w);
    
    // Compute right hand sides
    GetRHS(myPCSet,modelparamPCs(0).GetArrayPointer(),modelparamPCs(1).GetArrayPointer(),modelparamPCs(2).GetArrayPointer(),modelparamPCs(3).GetArrayPointer(),modelparamPCs(4).GetArrayPointer(),modelparamPCs(5).GetArrayPointer(),u.GetArrayPointer(),v.GetArrayPointer(),w.GetArrayPointer(),z.GetArrayPointer(),dudt.GetArrayPointer(),dvdt.GetArrayPointer(),dwdt.GetArrayPointer());

    // Advance u, v, w to mid-point
    myPCSet.Multiply(dudt,0.5*dTym,tmp_u); // 0.5*dTym*dudt
    myPCSet.Multiply(dvdt,0.5*dTym,tmp_v); // 0.5*dTym*dvdt
    myPCSet.Multiply(dwdt,0.5*dTym,tmp_w); // 0.5*dTym*dwdt

    myPCSet.Add(u_o,tmp_u,u); // u = u_o + 0.5*dTym*dudt
    myPCSet.Add(v_o,tmp_v,v); // v = v_o + 0.5*dTym*dvdt
    myPCSet.Add(w_o,tmp_w,w); // w = w_o + 0.5*dTym*dwdt
    
    // Compute z = 1 - u - v - w
    z=one;
    myPCSet.SubtractInPlace(z,u);
    myPCSet.SubtractInPlace(z,v);
    myPCSet.SubtractInPlace(z,w);
    
    // Compute right hand sides
    GetRHS(myPCSet,modelparamPCs(0).GetArrayPointer(),modelparamPCs(1).GetArrayPointer(),modelparamPCs(2).GetArrayPointer(),modelparamPCs(3).GetArrayPointer(),modelparamPCs(4).GetArrayPointer(),modelparamPCs(5).GetArrayPointer(),u.GetArrayPointer(),v.GetArrayPointer(),w.GetArrayPointer(),z.GetArrayPointer(),dudt.GetArrayPointer(),dvdt.GetArrayPointer(),dwdt.GetArrayPointer());
    
    // Advance u, v, w to next time step
    myPCSet.Multiply(dudt,dTym,tmp_u); // dTym*dudt
    myPCSet.Multiply(dvdt,dTym,tmp_v); // dTym*dvdt
    myPCSet.Multiply(dwdt,dTym,tmp_w); // dTym*dwdt
    
    
    myPCSet.Add(u_o,tmp_u,u); // u = u_o + dTym*dudt
    myPCSet.Add(v_o,tmp_v,v); // v = v_o + dTym*dvdt
    myPCSet.Add(w_o,tmp_w,w); // w = w_o + dTym*dwdt
    
    // Compute z = 1 - u - v - w
    z=one;
    myPCSet.SubtractInPlace(z,u);
    myPCSet.SubtractInPlace(z,v);
    myPCSet.SubtractInPlace(z,w);


    // Advance time and step counter
    tym += dTym;
    step+=1;
   

    // write time, u, v, w (all modes) to file
    if(step % outPrint->fdumpInt == 0){
      WriteModesToFilePtr(tym, u.GetArrayPointer(), v.GetArrayPointer(), w.GetArrayPointer(), nPCTerms, modes_dump);
    }

    // Get standard deviations
    uStDv = myPCSet.StDv(u);
    vStDv = myPCSet.StDv(v);
    wStDv = myPCSet.StDv(w);
    
    // write u, v, w (mean and standard deviation) to file
    if(step % outPrint->fdumpInt == 0){
      WriteMeanStdDevToFilePtr(tym, u(0), v(0), w(0), uStDv, vStDv, wStDv, f_dump);
    }
    // write u, v, w (mean and standard deviation) to screen
    if(step % outPrint->dumpInt == 0){
      WriteMeanStdDevToStdOut(step, tym, u(0), v(0), w(0), uStDv, vStDv, wStDv);
    }
    
  }
  
  // Close output file
  if(fclose(f_dump)){ 
    printf("Could not close file '%s'\n",outPrint->dumpfile.c_str()); 
    exit(1); 
  }


  // Close output file
  if(fclose(modes_dump)){ 
    printf("Could not close file '%s'\n",modes_dumpfile.c_str()); 
    exit(1); 
  }
  
  return 0;
}


