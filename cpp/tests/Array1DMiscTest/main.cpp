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

using namespace std; 

/*************************************************
Begin main code
*************************************************/
int main(int argc, char ** argv){

	/**********************************
	Read and write 1D Array
	*********************************/
	
	int dim = 3; 

	// Create dim-D array with all ones
	Array1D<double> x(dim,1);

	// Write array to file
	write_datafile_1d(x,"x.dat");

	// create dim-D array of zeros
	Array1D<double> y(dim,0);

	// read in data file to y
	read_datafile_1d(y,"x.dat");

	/**********************************
	Fill in normal r.v.'s to 1D Array
	*********************************/
	
	// Feed in normal random numbers to array
	dsfmt_t RandomState; 
	int seed = 1; 
	dsfmt_init_gen_rand(&RandomState,seed);
	for (int i = 0; i < dim; i++){
		x(i) = dsfmt_genrand_nrv(&RandomState);
	}

	// Write array to file
	write_datafile_1d(x,"x_nrv.dat");

	/**********************************
	Elementary operations on arrays
	*********************************/

	// add two arrays and print output
	Array1D<double> z = add(x,y); 
	printarray(z); 

	// multiply array by scalar and print output
	z = scale(z,3.14159); 
	printarray(z); 

	// dot product of array
	double a = dot(z,z); 
	cout << a << endl; 

	// delete ith element of array
	z.erase(1);
	printarray(z);

	// add element to end of array
	z.PushBack(1);
	printarray(z);

	return 0; 

}
