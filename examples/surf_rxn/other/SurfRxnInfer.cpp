/* =====================================================================================
                      The UQ Toolkit (UQTk) version @UQTKVERSION@
                     Copyright (2013) Sandia Corporation
                      http://www.sandia.gov/UQToolkit/

     Copyright (2013) Sandia Corporation. Under the terms of Contract DE-AC04-94AL85000
     with Sandia Corporation, the U.S. Government retains certain rights in this software.

     This file is part of The UQ Toolkit (UQTk)

     UQTk is free software: you can redistribute it and/or modify
     it under the terms of the GNU Lesser General Public License as published by
     the Free Software Foundation, either version 3 of the License, or
     (at your option) any later version.

     UQTk is distributed in the hope that it will be useful,
     but WITHOUT ANY WARRANTY; without even the implied warranty of
     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
     GNU Lesser General Public License for more details.
 
     You should have received a copy of the GNU Lesser General Public License
     along with UQTk.  If not, see <http://www.gnu.org/licenses/>.
 
     Questions? Contact Bert Debusschere <bjdebus@sandia.gov>
     Sandia National Laboratories, Livermore, CA, USA
===================================================================================== */

#include <math.h>
#include "XMLUtils.h"
#include "uqtktools.h"
#include "uqtkmcmc.h"

#include "model.h"
#include "posterior.h"
#include "XMLreader.h"

using namespace std;

/// \brief Example that infers parameters, given data on species concentrations, 
/// in a 3-equation model for heterogeneous surface
/// reaction involving a monomer, dimer, and inert species adsorbing
/// onto a surface out of gas phase. This model mimics some aspects
/// of CO oxidation.
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


/// Main program: MCMC inference of the ODE model parameters given data
int main(int argc, char *argv[])
{
  // Pointer to posterior information
  postAux* pinfo=new postAux;
  // Read the xml tree
  RefPtr<XMLElement> xmlTree=readXMLTree("surf_rxn.in.xml");
  // Read the model specification and send them to the posterior information structure
  readXMLModelInput(xmlTree,pinfo->modelparams, pinfo->modelparamnames, pinfo->modelauxparams);
  // Read specific information needed by inference, e.g. data and the parameters to be inferred
  readXMLDataInput(xmlTree,pinfo->data, pinfo->postparams,&(pinfo->noisetype));


  // Number of MCMC steps
  int nsteps;
  // Array to hold the starting values of the chain
  Array1D<double> chstart;
  // Define the MCMC object
  MCMC mchain(LogPosterior,(void*) pinfo);
  // Read the xml file for MCMC-specific information
  readXMLChainInput(xmlTree,&mchain, chstart, &nsteps,pinfo->chainParamInd,pinfo->priortype,pinfo->priorparam1,pinfo->priorparam2);

  // Prepend the parameter names to the output file
  FILE* f_out;  
  string filename=mchain.getFilename();
  int chdim=chstart.XSize();
  
  f_out = fopen(filename.c_str(),"w"); 

  fprintf(f_out, "%s ","Step");
  for(int i=0;i<chdim;i++)
    fprintf(f_out, "%21s ",pinfo->modelparamnames(pinfo->chainParamInd(i)).c_str());
  fprintf(f_out, "%24s %24s \n","Accept_prob","Log_posterior");
  fclose(f_out);

  // Set the flag
  mchain.namesPrepended();
  // Run the chain
  mchain.runChain(nsteps, chstart);
   
  return 0;
}


