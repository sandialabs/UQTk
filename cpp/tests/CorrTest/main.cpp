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
#include "arrayio.h"
#include "arraytools.h"
#include "probability.h"
#include "assert.h"


using namespace std; 

/*************************************************
Begin main code
*************************************************/
int main(int argc, char ** argv) {
  
  cout << "Testing implementation of distance correlation " << endl;
  double cMin[6]={3.6e-01,1.3e-01,1.4e-01,4.8e-01,3.0e-01,0.0};
  double cMax[6]={3.8e-01,1.5e-01,1.6e-01,5.0e-01,3.2e-01,0.01};

  Array2D<double> spls, dCor;

  for ( int i = 0; i < 6; i++ ) {
    char fname[20];
    sprintf(fname,"set_%d_w1.dat",i+1);
    read_datafileVS(spls, (char *) fname);
    distCorr( spls, dCor);
    printf("Set #%d: %e\n",i+1, dCor(1,0));
    assert((dCor(1,0)>cMin[i])&&(dCor(1,0)<cMax[i]));
  }

  return 0; 

}
