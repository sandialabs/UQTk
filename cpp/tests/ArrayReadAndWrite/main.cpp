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
