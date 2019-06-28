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
#include <cstdio>
#include <stddef.h>
#include <fstream>
#include <string>
#include <math.h>
#include <iostream>
#include <string>
#include <time.h>
#include <getopt.h>

#include "Array1D.h"
#include "Array2D.h"
#include "Array3D.h"
#include "error_handlers.h"
#include "PCSet.h"
#include "bcs.h"
#include "pce_bcs.h"
#include "arraytools.h"
#include "arrayio.h"
#include "tools.h"

using namespace std;





/// \brief Program to find sparse PC coefficients of a given set of data
int main(int argc, char *argv[])
{

  // Input settings
    int ntot=100;                    // Number of input samples
    int ndim=3;                      // The dimensionality of the input
    int nord=9;                      // The order of the initial, overcomplete basis
    int seed=13;                     // Seed for random sample generator
    string which_chaos="LU";          // PC type
    double eta=1.e-9;                // BCS stopping criterion


    // Generate input points
    Array2D<double> xdata(ntot,ndim,0.e0);
    generate_uniform(xdata,seed);
    
    // Evaluate a black-box forward function
    Array1D<double> ydata(ntot,0.e0);
    func(xdata,ydata);  

    // Compute the initial multiindex, i.e. the basis set up to given order for a given dimensionality
    Array2D<int> mindex;
    int npc=computeMultiIndex(ndim,nord,mindex);
  
    // Work variables
    Array1D<double> lambda_init ;                    // Parameter of the Laplace distribution, controlling the sparsity; 
                                                     // if empty, then lambda will be computed, otherwise lambda will be fixed to the given value
    double lambda ;                                  // Sparsity-controlling parameter on output
    Array1D<double> weights, errbars, basis, alpha ; // 
    Array1D<int> used ;                              // The indices of the selected basis terms
    Array2D<int> newmindex;                          // The new multiindex

    int    adaptive = 0      ;                       // Flag for adaptive implementation
    int    optimal  = 1      ;                       // Flag for optimal implementation
    double scale    = 0.1    ;                       // Diagonal loading parameter, relevant only in adaptive, non-optimal implementation
    int    verbose  = 0      ;                        // Verbosity
    double sigma2= (double) pow(get_std(ydata),2)/1.0e6;  // 'Data' noise variance, will be re-estimated on output


    // Declare a PC object with Non-intrusive, non-quadrature implementation
    PCSet PCCurrModel("NISPnoq",mindex,which_chaos,0.0,1.0);

    
    // Compute the projection matrix
    Array2D<double> Phi;
    PCCurrModel.EvalBasisAtCustPts(xdata, Phi); 
  
            
    // Clear the work variables
    lambda_init.Clear();
    weights.Clear();
    errbars.Clear();
    basis.Clear();
    alpha.Clear();
    used.Clear();
    newmindex.Clear();
 
    // Run the BCS algorithm
    BCS(Phi,ydata,sigma2,eta,lambda_init,adaptive, optimal,scale,
		verbose,weights,used,errbars,basis,alpha,lambda);
    
    // Pick the selected subset of the multiindices (bases)
    subMatrix_row(mindex,used,newmindex);

    cout << "BCS algorithm selected " << used.Length() << " out of " << npc << " terms" << endl;
    

    // Declare a new PC object with the sparse basis
    PCSet PCModel("NISPnoq",newmindex,which_chaos,0.0,1.0);
    // Evaluate the PC at the sampled points
    Array1D<double> ypc;
    PCModel.EvalPCAtCustPoints(ypc,xdata,weights);

    
    // Compute the relative L2 distance between values and the approximation
    double num=0.0, den=0.0;
    for(int i=0;i<ntot;i++){
      num+=pow(ydata(i)-ypc(i),2.0);
      den+=pow(ydata(i),2.0);
    }
    cout << "L2 relative error at the training points is " << sqrt(num/den) << endl;
    

    // Write out the sampled points
    write_datafile(xdata,"xdata.dat");
    // Write out the function evaluations
    write_datafile_1d(ydata,"ydata.dat");
    // Write out the initial multiindex
    write_datafile(mindex,"mindex.dat");
    // Write out PC coefficients 
    write_datafile_1d(weights,"PCcoeff.dat");
    // Write out the sparse multiindex
    write_datafile(newmindex,"mindex_new.dat");
    // Write out the indices of the bases selected from the initial multiindex
    write_datafile_1d(used,"used.dat");
    // Write out the sparse PC expansion values at the sampled points
    write_datafile_1d(ypc,"ypc.dat");
	  

   
    
  
    
  return 0;
}


// Forward function, f(x,y,z)=z*cos(x+y)
void func(Array2D<double>& xdata, Array1D<double>& ydata)
{


  int ntot=xdata.XSize();
  ydata.Resize(ntot,0.e0);
  for (int i=0;i<ntot;i++){
    ydata(i)=cos(xdata(i,0)+xdata(i,1))*xdata(i,2);
  }




  return;
}
