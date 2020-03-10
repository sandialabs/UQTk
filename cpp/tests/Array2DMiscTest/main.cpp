/* =====================================================================================

                      The UQ Toolkit (UQTk) version 3.1.0
                          Copyright (2020) NTESS
                        https://www.sandia.gov/UQToolkit/
                        https://github.com/sandialabs/UQTk

     Copyright 2020 National Technology & Engineering Solutions of Sandia, LLC (NTESS).
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

	int m = 3;
	int n = 3;

	// Create nxm array with all ones
	Array2D<double> A(m,n,1);

	// Write array to file
	write_datafile(A,"A.dat");

	// create dim-D array of zeros
	Array2D<double> B(m,n,0);

	// read in data file to B
	read_datafile(B,"A.dat");

	/**********************************
	Fill in normal r.v.'s to 1D Array
	*********************************/

	// Feed in uniform random numbers to array
	dsfmt_t RandomState;
	int seed = 1;
	dsfmt_init_gen_rand(&RandomState,seed);
	for (int i = 0; i < m; i++){
		for (int j = 0; j < n; j++){
			A(i,j) = dsfmt_genrand_urv(&RandomState);
		}
	}

	// Write array to file
	write_datafile(A,"A_nrv.dat");
	cout << "\nA : " << endl;
	printarray(A);

	/*********************************
	Linaer alg. operations on arrays
	********************************/

	// matrix vector product
	Array1D<double> x(3,1);
	Array1D<double> b = dot(A,x);
	cout << "\nb = " << endl;
	printarray(b);

	// print transpose
	Array2D<double> AT = Trans(A);
	cout << "\nTranspose:" << endl;
	printarray(AT);

	// get inverse of square matrix A
	if ( n == m){
		Array2D<double> Ainv = INV(A);
		cout << "\nA inverse:" << endl;
		printarray(Ainv);
	}

	// get least squares solution
	x.Resize(3,0);
	LSTSQ(A,b,x);
	cout << "\nleast squares solution:" << endl;
	printarray(x);

	// get QR factorization
	Array2D<double> Q,R;
	QR(A,Q,R);
	cout << "\nQ:" << endl;
	printarray(Q);
	cout << "\nR:" << endl;
	printarray(R);

	// get SVD factorization
	Array2D<double> U,VT;
	Array1D<double> S;
	SVD(A,U,S,VT);
	cout << "\nU:" << endl;
	printarray(U);
	cout << "\nS:" << endl;
	printarray(S);
	cout << "\nVT:" << endl;
	printarray(VT);


}
