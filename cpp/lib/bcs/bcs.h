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
/// \file bcs.h 
/// \brief Header for the implemenations of Bayesian compressive sensing algorithm. 

#ifndef BCS_H
#define BCS_H

#include "Array1D.h"
#include "Array2D.h"


#define MAX_IT 1000


/// \brief Given the projection matrix PHI, the measurement vector y, initial noise variance sigma2,
/// stopping criterion eta, hierarchical prior parameter lambda_init, adaptivity,optimality, scale and verbosity flags,
/// produces weights of the retained bases, their corresponding number (in the array 'used'), errorbars,
/// next projection basis (if adaptive), and estimates for prior hyperparameters alpha and lambda.
void WBCS(Array2D<double> &PHI, Array1D<double> &y, double &sigma2,
                 double eta, Array1D<double> &lambda_init, 
		         int adaptive, int optimal, double scale, int verbose,
                 Array1D<double> &weights, Array1D<int> &used,
                 Array1D<double> &errbars, Array1D<double> &basis, 
                 Array1D<double> &alpha, double &lambda, Array2D<double> &Sig);

/// \brief Given the projection matrix PHI, the measurement vector y, initial noise variance sigma2,
/// stopping criterion eta, hierarchical prior parameter lambda_init, adaptivity,optimality, scale and verbosity flags,
/// produces weights of the retained bases, their corresponding number (in the array 'used'), errorbars,
/// next projection basis (if adaptive), and estimates for prior hyperparameters alpha and lambda.
void BCS(Array2D<double> &PHI, Array1D<double> &y, double &sigma2,
                 double eta, Array1D<double> &lambda_init, 
                 int adaptive, int optimal, double scale, int verbose,
                 Array1D<double> &weights, Array1D<int> &used,
                 Array1D<double> &errbars, Array1D<double> &basis, 
                 Array1D<double> &alpha, double &lambda) ;


#endif // BCS_H
