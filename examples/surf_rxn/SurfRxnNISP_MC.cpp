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

#include <math.h>
#include "XMLUtils.h"
#include "tools.h"
#include "mcmc.h"
#include "arrayio.h"
#include "arraytools.h"
#include "PCSet.h"
#include "ss.h"
#include "model.h"
#include "XMLreader.h"
#include "Utils.h"
#include "amcmc.h"


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



/// Main program of uncertainty propagation of the ODE model parameters via non-intrusive spectral projection by Monte-Carlo (NISP_MC)
int main(int argc, char* argv[])
{
  // Number of MC samples
  int nSam;
  if (argc!=2)
    nSam=100;
  else
    nSam=atoi(argv[1]);


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

  // Instantiate a PC object for NISP computations with no quadrature
  PCSet myPCSet("NISPnoq",order,dim,pcType,0.0,1.0);

  // The number of PC terms
  const int nPCTerms = myPCSet.GetNumberPCTerms();
  cout << "The number of PC terms in an expansion is " << nPCTerms << endl;

  // Print the multiindices on screen
  myPCSet.PrintMultiIndex();

  // Generate Monte-Carlo samples
  Array2D<double> samPts(nSam,dim,0.e0);
  myPCSet.DrawSampleVar(samPts);

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

  // Evaluation the model at quadrature points
  Array2D<double> allmodelparams(nSam,fulldim,0.e0);
  for(int iq=0;iq<nSam;iq++){

    // Replace the model parameters corresponding to the sample points
      PCSet inpPC("NISPnoq",order,dim,pcType,0.0,1.0);
      Array1D<double> sampt(dim,0.e0);
      getRow(samPts,iq,sampt);
      for(int idim=0;idim<dim;idim++){
	Array1D<double> cur_cf;
	getCol(allPCcoefs,uncParamInd(idim),cur_cf);
	modelparams(uncParamInd(idim))=  inpPC.EvalPC(cur_cf, sampt);
      }

    printf("Model parameters for sample %d / %d: ", iq+1, nSam);
    for(int i=0;i<modelparams.XSize();i++){
      printf("%s=%lg ",modelparamnames(i).c_str(),modelparams(i));
    }
    printf("\n");
    for(int idim=0;idim<fulldim;idim++)
      allmodelparams(iq,idim)=modelparams(idim);
  }


    // Time and time step
    int step=0;
    double tym=t0;


 // Open files to write out
  FILE *f_dump,*modes_dump;
  // ... the mean and stdev
  if(!(f_dump = fopen(outPrint->dumpfile.c_str(),"w"))){
    printf("Could not open file '%s'\n",outPrint->dumpfile.c_str());
    exit(1);
  }

  // ... the full set of modes, filename hardwired
  string modes_dumpfile = "solution_NISP_MC_modes.dat";
  if(!(modes_dump = fopen(modes_dumpfile.c_str(),"w"))){
    printf("Could not open file '%s'\n",modes_dumpfile.c_str());
    exit(1);
  }



 // Model evaluations at the initial step
  Array1D<double> uu(nSam,0.e0),vv(nSam,0.e0),ww(nSam,0.e0);
  cout << "Stepping forward in time" << endl;

  while(tym < tf) {

    if (step % 1000 == 0)
      cout << (double) step*100.0/nStep << " % completed" << endl;


    for(int iq=0;iq<nSam;iq++){

      Array1D<double> outValues(3,0.e0);
      outValues(0)=uu(iq);
      outValues(1)=vv(iq);
      outValues(2)=ww(iq);
      Array1D<double> cur_modelparams;
      getRow(allmodelparams, iq, cur_modelparams);

      // Forward run, see model.h
      forwardFunctionDt(cur_modelparams.GetArrayPointer(),dTym,outValues.GetArrayPointer());

      // Update the output values
      uu(iq) = outValues(0);
      vv(iq) = outValues(1);
      ww(iq) = outValues(2);

    }


    // Non-intrusive spectral projection via Monte-Carlo integration
    myPCSet.GalerkProjectionMC(samPts,uu,u);
    myPCSet.GalerkProjectionMC(samPts,vv,v);
    myPCSet.GalerkProjectionMC(samPts,ww,w);


    // Get standard deviations
    double uStDv = myPCSet.StDv(u);
    double vStDv = myPCSet.StDv(v);
    double wStDv = myPCSet.StDv(w);

    if(step % outPrint->fdumpInt == 0){
      // write time, u, v, w (all modes) to file
      WriteModesToFilePtr(step*dTym, u.GetArrayPointer(), v.GetArrayPointer(), w.GetArrayPointer(), nPCTerms, modes_dump);
      // write u, v, w (mean and standard deviation) to file
      WriteMeanStdDevToFilePtr(step*dTym, u(0), v(0), w(0), uStDv, vStDv, wStDv, f_dump);
    }

    // write u, v, w (mean and standard deviation) to screen
    if(step % outPrint->dumpInt == 0 && step>0){
      WriteMeanStdDevToStdOut(step, step*dTym, u(0), v(0), w(0), uStDv, vStDv, wStDv);
    }

    // Step forward
      tym += dTym;
      step+=1;

  }// End of time loop






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
