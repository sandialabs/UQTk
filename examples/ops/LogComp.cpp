/* =====================================================================================
                      The UQ Toolkit (UQTk) version 3.0.4
                     Copyright (2017) Sandia Corporation
                      http://www.sandia.gov/UQToolkit/

     Copyright (2017) Sandia Corporation. Under the terms of Contract DE-AC04-94AL85000
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

#include <iostream>
#include "PCSet.h"
#include "Utils.h"

/// \brief Short program to test and demonstrate UQTk operations on PC expansions.
///
/// For a summary of the theory behind the implemented algorithms,
/// see the paper
/// Numerical challenges in the use of polynomial chaos representations for
/// stochastic processes, B.J Debusschere, H.N Najm, P.P Pébay, O.M Knio, 
/// R.G Ghanem, O.P.L Maître. Siam J Sci Comput, 2005, vol. 26(2) pp.698-719
/// http://dx.doi.org/10.1137/S1064827503427741

int main()
{
  // Separator in output
  cout << endl << std::string(70,'=') << endl;
  cout << "Initialization" << endl;
  cout << std::string(70,'=') << endl << endl;

  // Initialize Polynomial Chaos (PC) Set
  int ord = 5;
  int dim = 1;
  PCSet myPCSet("ISP",ord,dim,"HG");

  // Number of terms in the PC expansion
  const int npc = myPCSet.GetNumberPCTerms();
  cout << "The number of PC terms in an expansion is " << npc << endl;

  // Print the indices of all basis terms and their norms
  myPCSet.PrintMultiIndexNormSquared();

  // Set accuracy for Taylor series approaches to function evaluations
  myPCSet.SetTaylorTolerance(1.e-15);

  // mean and std dev of a
  double ma = 4.0; // Mean
  double sa = 0.5; // Std Dev

  // Separator in output
  cout << endl << std::string(70,'=') << endl;
  cout << "Computing Log(a)" << endl;
  cout << std::string(70,'=') << endl << endl;

  // Test of routines with Array arrays to hold the PC coefficients
  // of the random variables
  Array1D<double> aa(npc,0.e0);

  myPCSet.InitMeanStDv(ma,sa,aa);

  // Check initialization
  cout << "PC coeffs for aa are :" << endl;
  for(int ip=0; ip < npc; ip++){
    cout <<"  " << ip <<": "<< aa(ip) << endl;
  }

  // Check initialization
  double s_dum = myPCSet.StDv(aa);
  cout << endl;
  cout << "Mean of aa is    : " << aa(0) << endl;
  cout << "Std Dev of aa is : " << s_dum << endl;


  // Logarithm
  cout << endl << "Computing the natural log with Taylor Series approach:" << endl;
  Array1D<double> ac_t(npc,0.e0);
  myPCSet.SetLogCompMethod(TaylorSeries);
  myPCSet.Log(aa,ac_t);

  cout << "PC coeffs for c= ln(a) are (Array1D routine, Taylor Series):" << endl;
  for(int ip=0; ip < npc; ip++){
    cout <<"  " << ip <<": "<< ac_t(ip) << endl;
  }

  cout << endl << "Computing the natural log with Integration approach:" << endl;
  Array1D<double> ac_i(npc,0.e0);
  myPCSet.SetLogCompMethod(Integration);
  myPCSet.Log(aa,ac_i);

  cout << "PC coeffs for c= ln(a) are (Array1D routine, Integration):" << endl;
  for(int ip=0; ip < npc; ip++){
    cout <<"  " << ip <<": "<< ac_i(ip) << endl;
  }

  // Sample PCEs to generate PDFs
  cout << endl << std::string(70,'=') << endl;
  cout << "Sample PCEs to generate PDFs" << endl;
  cout << std::string(70,'=') << endl << endl;

  const int nsamples = 50000;

  Array1D<double> aa_samp(nsamples,0.e0);
  myPCSet.DrawSampleSet(aa, aa_samp);
  WriteToFile(aa_samp,(char *)"samples.a.dat");

  Array1D<double> act_samp(nsamples,0.e0);
  myPCSet.DrawSampleSet(ac_t, act_samp);
  WriteToFile(act_samp,(char *)"samples.loga_tay.dat");

  Array1D<double> aci_samp(nsamples,0.e0);
  myPCSet.DrawSampleSet(ac_i, aci_samp);
  WriteToFile(aci_samp,(char *)"samples.loga_int.dat");

  return 0;

}
