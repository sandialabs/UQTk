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
#include "arrayio.h"
#include "PCBasis.h"
#include "PCSet.h"
#include "arraytools.h"
#include "dsfmt_add.h"
#include "bcs.h"
#include "assert.h"


using namespace std; 


int main(){

	// get pc object
	Array1D<double> ck(10,0.0);
	Array2D<int> mindex(10,2); 

	// set ck values
	ck(0) = 0.666666666666664;
	ck(1) = 1.600000000000499; 
	ck(2) = 1.000000000000289;
	ck(5) = -0.6666666666668039;
	ck(6) = 0.4000000000008473; 

	// set mindex values
	mindex(0,0) = 0; mindex(0,1) = 0;
	mindex(1,0) = 1; mindex(1,1) = 0;
	mindex(2,0) = 0; mindex(2,1) = 1;
	mindex(3,0) = 2; mindex(3,1) = 0;
	mindex(4,0) = 1; mindex(4,1) = 1;
	mindex(5,0) = 0; mindex(5,1) = 2;
	mindex(6,0) = 3; mindex(6,1) = 0;
	mindex(7,0) = 2; mindex(7,1) = 1;
	mindex(8,0) = 1; mindex(8,1) = 2;  
	mindex(9,0) = 0; mindex(9,1) = 3;

	// set pc model given multiindex
	PCSet pcmodel("NISPnoq",mindex,"LEG");

	// get 2d quadrature points
	Quad q("LU","sparse",2,5);
	Array2D<double> x; 
	Array1D<double> w; 
	q.SetRule();
	q.GetRule(x,w);
	// printarray(x);

	// evaluate PCE at quadrature points
	Array1D<double> y(x.XSize(),0.0);
	pcmodel.EvalPCAtCustPoints(y,x,ck);
	// printarray(y);

	// get projection matrix
	Array2D<double> Phi; 
	pcmodel.EvalBasisAtCustPts(x,Phi);
	// printarray(Phi);

	// Main inputs are Phi, ydata, sigma
	double sigma = 1e-8; 

	// params
	double eta = 1e-12; 
	Array1D<double> lambda_init; 
	double scale = .1; 

	// outputs
	Array1D<double> weights, errbars, basis, alpha;
	Array1D<int> used; 
	double lambda=0.0;


	int adaptive=1;
	int optimal=1;
	int verbose=0;

	// run bcs
	//bcs(Phi,y,sigma,eta,scale,weights,used,errbars);
	BCS(Phi,y,sigma,eta,lambda_init,adaptive,optimal,scale,verbose,weights,used,errbars,basis,alpha,lambda);


	printarray(weights);
	printarray(used);
	printarray(mindex);
	printarray(ck);

	assert(used(0) == 1);
	assert(used(1) == 2);
	assert(used(2) == 0);

	assert(fabs(weights(0) - 1.6) < 1e-8);
	assert(fabs(weights(1) - 1) < 1e-8);
	assert(fabs(weights(2) - 2./3) < 1e-8);




	return 0; 

}
