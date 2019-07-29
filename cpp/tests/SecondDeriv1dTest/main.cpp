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
int main(int argc, char ** argv){

		/********************************************
	********************************************/

	// get 1d legendre polynomials
	int ndim = 3; 
	int norder = 4;
	PCSet polymodel("NISPnoq",norder,ndim,"LU"); 
	polymodel.PrintMultiIndex();

	// define alpha, which is a specific multiindex
	Array1D<int> alpha(ndim);
	alpha(0) = 1;
	alpha(1) = 3;
	alpha(2) = 2;

	// define evaluation point and gradient
	Array1D<double> x(3,1); 
	Array1D<double> grad; 

	// get derivative of 3d legendre polynomial 
	// defined by the multiindex above
	polymodel.dPhi_alpha(x, alpha, grad);
	printarray(grad);
 	assert(grad(0) == 1.0);
 	assert(grad(1) == 6.0);
 	assert(grad(2) == 3.0);

 	


	return 0; 

}
