/* =====================================================================================

                      The UQ Toolkit (UQTk) version 3.1.4
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

     Questions? Contact the UQTk Developers at https://github.com/sandialabs/UQTk/discussions
     Sandia National Laboratories, Livermore, CA, USA
===================================================================================== */
/// \file mrv.h
/// \author K. Sargsyan  2016 -
/// \brief Header for multivariate random variable class

#ifndef MRV_H_SEEN
#define MRV_H_SEEN

#include "Array1D.h"
#include "Array2D.h"
#include "PCSet.h"

#include <iostream>
#include <string.h>
#include <stdio.h>
#include <sstream>

using namespace std; // needed for python string conversion

/// \class	Mrv: class for multivariate RV parameterized by PC expansions

class Mrv {
public:

  /// \brief Constructor with dimensionality, pdftype, randomized parameter indices, order, and pctype
  Mrv(int ndim,string pdfType, Array1D<int> rndInd, int order,string pctype);
  /// \brief Destructor
  ~Mrv() {}
  

  /// \brief Parameterization bookkeeping (i.e. alpha corresponds to certain parameter lambda and certain PC term)
  int Parametrize();

  /// \brief Get bounds on parameters
  /// \note Useful when some parameters forced to be positive to make use of invariance
  void getBounds(Array1D<double>& lower, Array1D<double>& upper);

  /// \brief Get dimensionailty of parameterization
  int getPDim(){ return this->pDim_;}
   
  /// \brief Given parameters of representation, fold them in a 2d-array of PC coefficients for convenience
  Array2D<double> getMultiPCcf(Array1D<double>& rvParams);
  /// \brief Evaluate at multivariate PC at given germ samples for given coefficient matrix
  Array2D<double> evalMultiPC(Array2D<double>& xiSam, Array2D<double>& multiPCcf);
  /// \brief Random-sample all parameters given coefficient matrix
  Array2D<double> mcParam(Array2D<double>& multiPCcf, int nsam);
  /// \brief Quadrature-sample all parameters given coefficient matrix
  Array2D<double> quadParam(Array2D<double>& multiPCcf);
  /// \brief Propagate the multivariate RV with given coefficeints through a given function at given values x
  Array2D<double> propNISP(Array2D<double> (*forwardFcn)(Array2D<double>&, Array2D<double>&, Array2D<double>&, void*), Array2D<double>& fixindnom,void* funcinfo, Array2D<double>& multiPCcf, Array2D<double>& x);
  /// \brief Sample values of a given function given input coefficeint matrix
  Array2D<double> propMC(Array2D<double> (*forwardFcn)(Array2D<double>&, Array2D<double>&, Array2D<double>&, void*), Array2D<double>& fixindnom,void* funcinfo,Array2D<double>& multiPCcf, Array2D<double>& x,int nsam);
  /// \brief Compute moments given coefficent matrix
  void computeMoments(Array2D<double>& funcCf, Array1D<double>& fcnMean,Array1D<double>& fcnStd,bool covFlag, Array2D<double>& fcnCov);

  /// \brief Get PC term ID
  void getPCTermId(Array1D<int>& pctermid){pctermid=pctermId_; return;}

private:
  /// \brief Randomized parameters indices
  Array1D<int> rndInd_;
  /// \brief For a given parameterization, id the corresponding physical parameter lambda
  Array1D<int> paramId_;
  /// \brief For a given parameterization, id the PC term/order for the corresponding parameter representation
  Array1D<int> pctermId_;

  /// \brief PDF type ('pct', 'pci' or 'full')
  string pdfType_;
  /// \brief PC type (see pce library for options)
  string pcType_;
  /// \brief Number of parameters in alpha parameterization
  int pDim_;
  /// \brief Number of randomized parameters
  int rDim_;
  /// \brief Number of physical parameters lambda
  int nDim_;
  /// \brief Order of function PC representation
  int order_;
  /// \brief Number of PC parameters for each independent component
  int nPC_;
  /// \brief Pointer to the corresponding PC object
  PCSet* pcModel_;


  
}; 

/*******************************************************************/
/*******************************************************************/
/*******************************************************************/


#endif /* MRV_H_SEEN */
