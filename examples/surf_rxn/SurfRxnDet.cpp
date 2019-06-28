/* =====================================================================================
                      The UQ Toolkit (UQTk) version 3.0.4
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
#include "PCSet.h"
#include "tools.h"
#include "mcmc.h"
#include "arrayio.h"
#include "arraytools.h"

#include "model.h"
#include "XMLreader.h"

using namespace std;


/// \brief Example that solves deterministic 3-equation model for heterogeneous surface
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


/// Main program: a single run of the ODE model given parameters in an xml file
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
  // Read the output preferences
  dumpInfo* outprint=new dumpInfo;
  readXMLDumpInfo( xmlTree, &(outprint->dumpInt), &(outprint->fdumpInt), &(outprint->dumpfile) );

  // Initial time
  double t0 = 0.0;
  // Final time
  double tf = modelauxparams(0);
  // Time step
  double dTym = modelauxparams(1);

  // Initial conditions of zero coverage (based on Makeev:2002)
  double outValues[]={0.e0,0.e0,0.e0};
  
  // Forward ODE function, see model.h
  forwardFunction(modelparams.GetArrayPointer(), t0, tf, dTym, outValues,outprint);


  return 0;
}

