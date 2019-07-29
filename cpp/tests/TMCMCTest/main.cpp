/* =====================================================================================
                     The UQ Toolkit (UQTk) version @UQTKVERSION@
                     Copyright (@UQTKYEAR@) Sandia Corporation
                     http://www.sandia.gov/UQToolkit/

     Copyright (@UQTKYEAR@) Sandia Corporation. Under the terms of Contract DE-AC04-94AL85000
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
#include "math.h"
#include "Array1D.h"
#include "Array2D.h"
#include "quad.h"
#include "dsfmt_add.h"
#include "arrayio.h"
#include "arraytools.h"
#include "assert.h"
#include "tmcmc2.h"

using namespace std; 

typedef std::vector<double>  RealVector;


/*************************************************
Begin main code
*************************************************/
int main(int argc, char ** argv){

	// Default Arguments
	int iseed  = 1;
	int ndim   = 5; // update rngs.dat with initial proposal (uniform) bounds
	int nprocs = 4;
	int nspl   = 10000;
	double gamma = -1.0;
	int maxniter = 1;
	double deltaRate = 0.8;
	double cv = 0.1;
	double a = 1.0/9.0;
	double b = 8.0/9.0;
	int MFactor = 1;
	bool basis = false;
	int CATSteps = 1;
	double logevid;


	RealVector rngs;
	std::ifstream myfile;
	double dtmp;
	myfile.open("rngs.dat", std::ios_base::in);
	while (myfile >> dtmp) rngs.push_back(dtmp);
	myfile.close();


RealVector dts; // Obsolete, Dbeta changes are adaptive now.
  std::ofstream evidFile("Evidence.dat");

  // Run TMCMC (Repeat for Manifolds)
  for (size_t i = 0; i < maxniter; ++i) {



    // Track Stage
    std::ofstream TMCMCiter("TMCMCiter.dat");
    TMCMCiter.clear();
    TMCMCiter << i;
    TMCMCiter.close();


    // Run TMCMC, get evidence
    logevid = tmcmc2(rngs, gamma, nspl, iseed, nprocs, ndim,
                        cv, a, b, MFactor, basis, CATSteps);

    evidFile << std::setprecision(18) << logevid << std::endl;
    // Move residual files out of the way if necessary.
    std::string filename;
    for (size_t f = 1000; f > 0; --f) {
      filename = "samples.dat." + std::to_string(f);

      std::ifstream sampleFile(filename);

      if (sampleFile.is_open()) {
        break;
      }
    }

    std::ifstream moveFile("move.sh");
    if (moveFile.is_open()) {
      std::string moveStr = "./move.sh Stage" + std::to_string(i) +
      " " + filename;
      system(moveStr.c_str());
    }
  }

  evidFile.close();

  // true log-evidence is -4.79991 = log(2*(1/3)^5)
  double err_threshold = 1.0; // threshold on error in log-evidence
  assert(fabs(logevid - (-4.79991)) < err_threshold);

  return 0; 

}
