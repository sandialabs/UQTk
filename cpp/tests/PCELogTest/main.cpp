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

/// \brief Short program to test and compare the Taylor Series and
/// Integration approach for computing the Log of a PCE
///
/// For a summary of the theory behind the implemented algorithms,
/// see the paper
/// Numerical challenges in the use of polynomial chaos representations for
/// stochastic processes, B.J Debusschere, H.N Najm, P.P Pébay, O.M Knio,
/// R.G Ghanem, O.P.L Maître. Siam J Sci Comput, 2005, vol. 26(2) pp.698-719
/// http://dx.doi.org/10.1137/S1064827503427741

#include <iostream>
#include "math.h"
#include "Array1D.h"
#include "Array2D.h"
#include "dsfmt_add.h"
#include "arrayio.h"
#include "arraytools.h"
#include "PCBasis.h"
#include "PCSet.h"
#include "quad.h"
#include "assert.h"


using namespace std;

/*************************************************
Begin main code
*************************************************/

int main()
{
  // Separator in output
  cout << endl << std::string(70,'=') << endl;
  cout << "Initialization" << endl;
  cout << std::string(70,'=') << endl << endl;

  // Initialize Polynomial Chaos (PC) set with Wiener-Hermite basis functions
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

	// Compute the rms difference between the PC coefficients obtained with
	// both methods
  cout << endl << std::string(70,'=') << endl;
	cout << "Assessing the difference between both approaches:" << endl;
	cout << std::string(70,'=') << endl << endl;

	double error = 0.0;
	for(int ip=0; ip < npc; ip++){
    error += pow(ac_i(ip)-ac_t(ip),2);
  }

	cout << "RMS difference between Integration and Taylor Series coefficients for Log(p): " << sqrt(error/float(npc)) << endl;
	assert( sqrt(error) <= 1e-10);

	// Conclusion
  cout << endl << std::string(70,'=') << endl;
	cout << "The test was successful!" << endl;
	cout << std::string(70,'=') << endl << endl;
	return 0;

}
