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

     Questions? Contact the UQTk Developers at https://github.com/sandialabs/UQTk/discussions
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
  int ord = 2;
  int dim = 1;
  PCSet myPCSet("ISP",ord,dim,"HG");

  // Number of terms in the PC expansion
  const int npc = myPCSet.GetNumberPCTerms();
  cout << "The number of PC terms in an expansion is " << npc << endl;

  // Print the indices of all basis terms and their norms
  myPCSet.PrintMultiIndexNormSquared();

  // Set accuracy for Taylor series approaches to function evaluations
  myPCSet.SetTaylorTolerance(1.e-15);

  // Separator in output
  cout << endl << std::string(70,'=') << endl;
  cout << "Arithmetic with routines that take double* arguments" << endl;
  cout << std::string(70,'=') << endl << endl;

  // Test of routines with double* arrays to hold
  // the PC coefficients of the random variables
  double* a = new double[npc];
  double* b = new double[npc];
  double* c = new double[npc];

  for(int ip=0; ip < npc; ip++){
    a[ip] = 0.e0;
    b[ip] = 0.e0;
    c[ip] = 0.e0;
  }

  // Initialize PC coefficients for the 
  // random variables a and b
  double ma = 2.0; // Mean
  double sa = 0.1; // Std Dev
  myPCSet.InitMeanStDv(ma,sa,a);

  // Check initialization
  cout << "PC coeffs for a are :" << endl;
  for(int ip=0; ip < npc; ip++){
    cout <<"  " << ip <<": "<< a[ip] << endl;
  }

  // Check initialization
  double s_dum = myPCSet.StDv(a);
  cout << endl;
  cout << "Mean of a is    : " << a[0] << endl;
  cout << "Std Dev of a is : " << s_dum << endl;

  // Initialize b
  b[0] = 2.0;
  b[1] = 0.2;
  b[2] = 0.01;

  // Perform some arithmetic on the random variables
  cout << endl << "Subtraction:" << endl;
  myPCSet.Subtract(a,b,c);

  cout << "PC coeffs for c= a-b are (double* routine):" << endl;
  for(int ip=0; ip < npc; ip++){
    cout <<"  " << ip <<": "<< c[ip] << endl;
  }

  cout << endl << "Product:" << endl;
  myPCSet.Prod(a,b,c);

  cout << "PC coeffs for c= a*b are (double* routine):" << endl;
  for(int ip=0; ip < npc; ip++){
    cout <<"  " << ip <<": "<< c[ip] << endl;
  }

  cout << endl << "Exponential:" << endl;
  myPCSet.Exp(a,c);

  cout << "PC coeffs for c= exp(a) are (double* routine):" << endl;
  for(int ip=0; ip < npc; ip++){
    cout <<"  " << ip <<": "<< c[ip] << endl;
  }

  cout << endl << "Division:" << endl;
  myPCSet.Div(a,b,c);

  cout << "PC coeffs for c= a/b are (double* routine):" << endl;
  for(int ip=0; ip < npc; ip++){
    cout <<"  " << ip <<": "<< c[ip] << endl;
  }

  cout << endl << "Standard Deviation: " << endl;
  double c_stdv = myPCSet.StDv(c);
  cout << "Std Dev of c=a/b is (double* routine):" << c_stdv << endl;

  cout << endl << "Natural Log: " << endl;
  myPCSet.Log(a,c);

  cout << "PC coeffs for c= ln(a) are (double* routine):" << endl;
  for(int ip=0; ip < npc; ip++){
    cout <<"  " << ip <<": "<< c[ip] << endl;
  }

  delete[] a;
  delete[] b;
  delete[] c;

  // Separator in output
  cout << endl << std::string(70,'=') << endl;
  cout << "Arithmetic with routines that take Array arguments" << endl;
  cout << std::string(70,'=') << endl << endl;

  // Test of routines with Array arrays to hold the PC coefficients
  // of the random variables
  Array1D<double> aa(npc,0.e0);
  Array1D<double> ab(npc,0.e0);
  Array1D<double> ac(npc,0.e0);

  // Initialize some of the PC coefficients
  aa(0) = 2.0;
  aa(1) = 0.1;
  ab(0) = 2.0;
  ab(1) = 0.2;
  ab(2) = 0.01;

  // Perform arithmetic operations on the random variables
  cout << endl << "Subtraction: " << endl;
  myPCSet.Subtract(aa,ab,ac);

  cout << "PC coeffs for c= a-b are (Array1D routine):" << endl;
  for(int ip=0; ip < npc; ip++){
    cout <<"  " << ip <<": "<< ac(ip) << endl;
  }

  cout << endl << "Product: " << endl;
  myPCSet.Prod(aa,ab,ac);

  cout << "PC coeffs for c= a*b are (Array1D routine):" << endl;
  for(int ip=0; ip < npc; ip++){
    cout <<"  " << ip <<": "<< ac(ip) << endl;
  }

  cout << endl << "Exponential: " << endl;
  myPCSet.Exp(aa,ac);

  cout << "PC coeffs for c= exp(a) are (Array1D routine):" << endl;
  for(int ip=0; ip < npc; ip++){
    cout <<"  " << ip <<": "<< ac(ip) << endl;
  }

  cout << endl << "Division: " << endl;
  myPCSet.Div(aa,ab,ac);

  cout << "PC coeffs for c= a/b are (Array1D routine):" << endl;
  for(int ip=0; ip < npc; ip++){
    cout <<"  " << ip <<": "<< ac(ip) << endl;
  }

  cout << endl << "Standard Deviation: " << endl;
  double ac_stdv = myPCSet.StDv(ac);
  cout << "Std Dev of c=a/b is (Array1D routine):" << ac_stdv << endl;

  // Logarithm
  cout << endl << "Computing the natural log with Taylor Series approach:" << endl;
  myPCSet.SetLogCompMethod(TaylorSeries);
  myPCSet.Log(aa,ac);

  cout << "PC coeffs for c= ln(a) are (Array1D routine, Taylor Series):" << endl;
  for(int ip=0; ip < npc; ip++){
    cout <<"  " << ip <<": "<< ac(ip) << endl;
  }

  cout << endl << "Computing the natural log with Integration approach:" << endl;
  myPCSet.SetLogCompMethod(Integration);
  myPCSet.Log(aa,ac);

  cout << "PC coeffs for c= ln(a) are (Array1D routine, Integration):" << endl;
  for(int ip=0; ip < npc; ip++){
    cout <<"  " << ip <<": "<< ac(ip) << endl;
  }

  // Sample PCEs to generate PDFs
  cout << endl << std::string(70,'=') << endl;
  cout << "Sample PCEs to generate PDFs" << endl;
  cout << std::string(70,'=') << endl << endl;

  const int nsamples = 50000;

  Array1D<double> aa_samp(nsamples,0.e0);
  myPCSet.DrawSampleSet(aa, aa_samp);
  WriteToFile(aa_samp,(char *)"samples.a.dat");

  Array1D<double> ac_samp(nsamples,0.e0);
  myPCSet.DrawSampleSet(ac, ac_samp);
  WriteToFile(ac_samp,(char *)"samples.loga.dat");

  return 0;

}
