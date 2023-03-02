/* =====================================================================================

                      The UQ Toolkit (UQTk) version 3.1.3
                          Copyright (2023) NTESS
                        https://www.sandia.gov/UQToolkit/
                        https://github.com/sandialabs/UQTk

     Copyright 2023 National Technology & Engineering Solutions of Sandia, LLC (NTESS).
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

     Questions? Contact the UQTk Developers at <uqtk-developers@software.sandia.gov>
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
