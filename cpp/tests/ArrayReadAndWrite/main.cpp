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
	Array1D Test
	*************************************************/
	// create empty 1D array
	int N = 100; 
	Array1D<double> x(N,0);

	// fill in with integer variables
	for (int i = 0; i < N; i++){
		x(i) = i;
	}

	// write array to file
	write_datafile_1d(x,"array1d.txt");

	// read back data file
	Array1D<double> xnew(N,0);
	read_datafile_1d(xnew,"array1d.txt");

	// check to make sure xnew is right
	for (int i = 0; i < N; i++){
		assert(xnew(i) == i);
	}
	

	/*************************************************
	Array2D Test
	*************************************************/
	// create empty 1D array
	int m = 100; 
	int n = 3; 
	Array2D<double> x2d(m,n,0);

	// fill in with uniform random numbers
	for (int i = 0; i < m; i++){
		for (int j = 0; j < n; j++){
			x2d(i,j) = i*j;
		}
	}

	// write array to file
	write_datafile(x2d,"array2d.txt");

	Array2D<double> y2d(m,n);
	read_datafile(y2d,"array2d.txt");

	// check to make sure y2d is right
	for (int i = 0; i < m; i++){
		for (int j = 0; j < n; j++){
			assert(y2d(i,j) == i*j);
		}
	}

	return 0; 

}
