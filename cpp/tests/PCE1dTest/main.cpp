/* =====================================================================================

                      The UQ Toolkit (UQTk) version 3.1.1
                          Copyright (2021) NTESS
                        https://www.sandia.gov/UQToolkit/
                        https://github.com/sandialabs/UQTk

     Copyright 2021 National Technology & Engineering Solutions of Sandia, LLC (NTESS).
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
#include "PCBasis.h"
#include "PCSet.h"
#include "quad.h"
#include "assert.h"


using namespace std;

/*************************************************
Begin main code
*************************************************/
int main(int argc, char ** argv){

	// Set up quadrature rule
	Array2D<double> x;
	Array1D<double> w;
	Array2D<int> index;

	int ndim = 1;
	int level = 16;

	Quad q("LU","full",ndim,level);
	q.SetRule();
	q.GetRule(x,w,index);

	// define function to be evaluated
	Array1D<double> y(x.XSize(),0);
	for (int i = 0; i < x.XSize(); i++){
		y(i) = 1./(1 + x(i,0)*x(i,0));
		// cout << x(i,0) << ", " << y(i) << endl;
	}

	// Define PCSet object with quadrature rule above
	int nord = level;
	PCSet pcmodel("NISPnoq",nord,ndim,"LU");
	pcmodel.SetQuadRule(q);

	// get multiindex
	Array2D<int> mindex;
	pcmodel.GetMultiIndex(mindex);

	// get coefficients
	Array1D<double> ck;
	pcmodel.GalerkProjection(y,ck);

	// evaluate PC at quadrature points
	Array1D<double> ytest(x.XSize(),0);
	pcmodel.EvalPCAtCustPoints(ytest,x,ck);
	double error = 0.0;
	for (int i = 0; i < x.XSize(); i++){
		error += pow(y(i) - ytest(i),2);
	}
	cout << "error at the quadrature points: " << sqrt(error) << endl;
	assert( sqrt(error) <= 1e-12);



	return 0;

}
