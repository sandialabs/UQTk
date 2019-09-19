/* =====================================================================================

                      The UQ Toolkit (UQTk) version @UQTKVERSION@
                          Copyright (@UQTKYEAR@) NTESS
                        https://www.sandia.gov/UQToolkit/
                        https://github.com/sandialabs/UQTk

     Copyright @UQTKYEAR@ National Technology & Engineering Solutions of Sandia, LLC (NTESS).
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
#include "lowrank.h"
#include "tools.h"
#include "arrayio.h"
#include "arraytools.h"
#include "assert.h"


using namespace std;


int main(){


  /// Create xdata as training and xcheck as test input data
  int nsam    = 113;
  int ncheck  = 11;
  int ndim    = 2;
  int seed    = 111;
  int rank    = 2;
  int order   = 3;
  int maxiter = 50;
  Array2D<double> xdata(nsam,ndim);
  generate_uniform(xdata,seed);
  Array2D<double> xcheck(ncheck,ndim);
  generate_uniform(xcheck,seed*seed);


  /// Generate true function evaluations
  Array1D<double> ydata(nsam,0.e0);
  for(int i=0;i<nsam;i++)
    ydata(i)=xdata(i,0)*xdata(i,1)+pow(xdata(i,0),3)*pow(xdata(i,1),2);

  Array1D<double> ycheck(ncheck,0.e0);
  for(int i=0;i<ncheck;i++)
    ycheck(i)=xcheck(i,0)*xcheck(i,1)+pow(xcheck(i,0),3)*pow(xcheck(i,1),2);

  //write_datafile_1d(ycheck,"ycheck.txt");

  /// Work arrays
  Array1D<double> ycheck_lr,error;//(nsam,0.e0);

  /// Declare tensor object
  CanonicalTensor f;

  /// Initialize an array of functional bases
  Array1D<FunctionalBases* > testFArr(1);

  /// Set monomial type of basis in each dimension
  FunctionalBases *newBase;
  Array1D<int> baseOrder(ndim,order);
  newBase = new PLBases(baseOrder);
  testFArr(0) = newBase;


  /// Initialize algorithm
  CanonicalTensorALSLeastSquares s;
  s.setBases(testFArr);
  s.setIteration(maxiter);

  /// Loop over one to rank
  for (int r=1;r<=rank;r++){
    s.setRank(r);
    f = s.solveDirect(xdata,ydata,false);

    FunctionalTensor fun(f,testFArr);
    ycheck_lr = fun.tensorEval(xcheck);

    //char filename[10];
    //sprintf(filename, "ycheck_%d.txt",r);
    //write_datafile_1d(ycheck_lr,filename);

    error = subtract(ycheck,ycheck_lr);
    printf("Norm of the error for rank %d approximation = %e\n",r,(norm(error)/norm(ydata)));
  }
  delete newBase;

  // The actual function is rank=2, so we should have a negligible error
  assert((norm(error)/norm(ydata))<1e-8);



	return 0;

}
