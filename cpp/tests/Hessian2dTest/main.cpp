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
	Get Legendre PCE for 2d quadratic
	********************************************/
	// Set up quadrature rule
	Array2D<double> xq;
	Array1D<double> wq;
	Array2D<int> index; 

	int ndim1 = 2; 
	int level1 = 5;

	Quad q("LU","sparse",ndim1,level1);
	q.SetRule();
	q.GetRule(xq,wq,index);

	// define function to be evaluated
	Array1D<double> y0(xq.XSize(),0);
	for (int i = 0; i < xq.XSize(); i++){
		y0(i) = ( -.5*( sqrt(2)*xq(i,0)*xq(i,0) + sqrt(3)*xq(i,1)*xq(i,1)) + xq(i,0)*xq(i,0)*xq(i,1) );
	} 

	// Define PCSet object with quadrature rule above
	int nord = 3; 
	PCSet pcmodel("NISPnoq",nord,ndim1,"LU");
	pcmodel.SetQuadRule(q);

	// get multiindex
	Array2D<int> mindex0; 
	pcmodel.GetMultiIndex(mindex0);
	pcmodel.PrintMultiIndex();
	// write_datafile(mindex0, "mindex.dat");

	// get coefficients
	Array1D<double> ck0; 
	pcmodel.GalerkProjection(y0,ck0);
	// write_datafile_1d(ck0, "ck.dat");

	/********************************************
	Test Hessian of Legendre PCE
	********************************************/

	// get 1d legendre polynomials
	int ndim = 2;

	// The first step, after reading in mindex, is to initialize the PCSet object. 
	PCSet polymodel("NISPnoq",mindex0,"LU");
	polymodel.PrintMultiIndex();

	//compute the hessian at a single x point
	Array1D<double> xpnt(ndim,0.0);
	Array2D<double> hessian; 
	polymodel.ddPhi(xpnt,mindex0,hessian,ck0);
	cout << fabs(hessian(0,0) - (-sqrt(2))) << endl;
	cout << fabs(hessian(1,1) - (-sqrt(3))) << endl;
	cout << fabs(hessian(0,1) - (2*xpnt(0))) << endl;
	printarray(hessian);
	assert(fabs(hessian(0,0) - (-sqrt(2))) <= 1e-12);
	assert(fabs(hessian(1,1) - (-sqrt(3))) <= 1e-12);
	assert(fabs(hessian(0,1) - (2*xpnt(0))) <= 1e-12);


	return 0; 

}
