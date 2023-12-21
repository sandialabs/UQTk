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

#include <math.h>
#include "tools.h"
#include "mcmc.h"
#include "amcmc.h"
#include "arrayio.h"
#include "XMLUtils.h"
#include "ss.h"
#include "model.h"
#include "posterior.h"
#include "XMLreader.h"

using namespace std;

/// \brief Example that infers parameters of a model given data measureements on its outputs


/// Main program: MCMC inference of the model parameters given data
int main(int argc, char *argv[])
{
  // Pointer to posterior information
  postAux* pinfo=new postAux;
  // Read the xml tree
  RefPtr<XMLElement> xmlTree=readXMLTree("line_infer.xml");
  // Read the model specification and send them to the posterior information structure
  readXMLModelInput(xmlTree,pinfo->modelparams, pinfo->modelparamnames, pinfo->modelauxparams);
  // Read specific information needed by inference, e.g. data and the parameters to be inferred
  readXMLDataInput(xmlTree,pinfo->data, pinfo->postparams,&(pinfo->noisetype));


  // Number of MCMC steps
  int nsteps;
  // Array to hold the starting values of the chain
  Array1D<double> chstart;
  // Define the MCMC object
    ///\note It appears based on the .xml file that the MCMC passed in is of the AMCMC variety
  AMCMC mchain(LogPosterior,(void*) pinfo);
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
