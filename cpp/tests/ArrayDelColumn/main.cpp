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
#include "math.h"
#include "Array1D.h"
#include "Array2D.h"
#include "mcmc.h"
#include "quad.h"
#include "dsfmt_add.h"
#include "arrayio.h"
#include "arraytools.h"
#include "assert.h"

using namespace std; 

/*************************************************
Begin main code
*************************************************/
int main(int argc, char ** argv){


	/*************************************************
	Array2D Test
	This test deletes the last column from a 2d array
	*************************************************/
	// create empty 1D array
	int m = 100; 
	int n = 3; 
	Array2D<double> x2d(m,n,0);

	// fill in with uniform random numbers
	dsfmt_t RandomState2; // define state
	int seed2 = 10; 
	dsfmt_init_gen_rand(&RandomState2,seed2);
	for (int i = 0; i < m; i++){
		for (int j = 0; j < n; j++){
			x2d(i,j) = dsfmt_genrand_urv(&RandomState2);
		}
	}

	// write array to file
	write_datafile(x2d,"array2d.txt");

	// delete last column of 2d array
	Array2D<double> y2d = mtxdel(x2d,2,1);
	assert(y2d.YSize() == n-1);

	return 0; 

}
