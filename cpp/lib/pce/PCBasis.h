/* =====================================================================================

                      The UQ Toolkit (UQTk) version 3.1.2
                          Copyright (2022) NTESS
                        https://www.sandia.gov/UQToolkit/
                        https://github.com/sandialabs/UQTk

     Copyright 2022 National Technology & Engineering Solutions of Sandia, LLC (NTESS).
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
/// \file PCBasis.h
/// \author B. Debusschere, C. Safta, K. Sargsyan, K. Chowdhary 2007 -
/// \brief Header file for the univariate PC class

#ifndef PCBASIS_H_SEEN
#define PCBASIS_H_SEEN

#include <iostream>
#include <string.h>
#include "Array1D.h"
#include "Array2D.h"
// #include "Array3D.h"
#include "ftndefs.h"
#include "dsfmt_add.h"


/// \class  PCBasis
/// \brief  Contains all basis type specific definitions and operations
///         needed to generate a PCSet
class PCBasis {
public:
  /// \brief Constructor: initializes the univariate basis type and order
  ///
  /// Currently, the only valid types are Hermite-Gaussian, denoted with "HG",
  /// Legendre-Uniform, denoted with "LU", or Laguerre-Gamma, denoted with "LG".
  /// (Where the shape parameter for the Gamma distribution is alpha + 1 = 2)
  /// \todo At some point, the basis selection should probably be implemented
  /// in a more elegant way using base and inherited classes. For the time being,
  /// Hermite-Gaussian or Legendre-Uniform will probably be the most commonly used
  /// cases.
  /// The parameters alpha and betta are relevant only for LG, SW and JB chaoses
  /// \note Maxord specifies the maximal order up to which the computations are performed
  PCBasis(const string type="LU", const double alpha=0.0, const double betta=1.0, const int maxord=10);


  /// \brief Destructor
  ~PCBasis() {};

  /// \brief Initialize the quadrature points and weights and store the information
  /// in arrays quadPoints_, quadWeights_,quadIndices_
  /// \note Uses an arbitrary number of quad. points.
  /// \note The default implementation relies on N_q=2*p+1 quadrature points,
  /// where p is the maximal order and N_q is the number of quadrature points
  /// \todo Come up with a smarter way to pick the number of quadrature points
  /// \note Quadrature points are set according to the basis function type
  /// \note quadPoints is a 2D array but its second dimension is equal to 1.
  void Init1dQuadPoints(int qdpts);

  /// \brief Evaluate polynomial 1d basis functions at quadrature points
  /// and store in the private variable psi1d_
  void Eval1dBasisAtQuadPoints();

  /// \brief Evaluate polynomial 1d basis functions up to the order kord at custom points
  /// given by an array custPoints
  /// Returns the evaluations in the first argument psi, where the number of rows are
  /// the number of points, and columns correspond to successive orders
  void Eval1dBasisAtCustPoints(Array2D<double>& psi,int kord, const Array1D<double>& custPoints);

  /// \brief Evaluate 1d basis functions for the given value of random
  /// variable xi. Return the value of the basis functions for all orders
  /// in the passed Array1D array (indexed by their order),
  /// also returns the highest-order value.
  /// \note For custom 'pdf' option, a file containing the polynomial recursion
  /// coefficients, called 'ab.dat', is required.
  /// \todo Import the recursion coefficients in a more friendly fashion.
  double EvalBasis(const double &xi, Array1D<double> &basisEvals) const;
  /// \brief Evaluate 1d basis functions for the given value of random
  /// variable xi. Return the value of the basis functions for all orders
  /// in the passed double * array (indexed by their order),
  /// also returns the highest-order value.
  double EvalBasis(const double &xi, const int kord, double *basisEvals) const;

  /// \brief Evaluate the norms (squared) of the basis functions exactly
  /// and stores in the private array psi1dSqExact_
  void Eval1dNormSq_Exact(int kord);


  /***************************************************
  New derivative functionality
  ***************************************************/
  /// \brief Evaluate derivative of 1d non-normalized Legendre basis.
  void EvalDerivBasis(const double& xi, Array1D<double>& basisDEvals);
  void Eval1dDerivBasisAtCustPoints(Array2D<double>& dpsi,int kord, const Array1D<double>& custPoints);

  void Eval2ndDerivBasis(const double& xi,Array1D<double>& ddP);
  void Eval2ndDerivCustPoints(Array2D<double>& psi, int kord, Array1D<double>& custPoints);
  /***************************************************
  ***************************************************/

  /// \brief Get the norms-squared of the basis functions.
  /// Returns the values for each basis function in the passed Array1D array
  void Get1dNormsSq(Array1D<double>& psi1dSq) const {psi1dSq=psi1dSq_; return;}

  /// \brief Get the analytic norms-squared of the basis functions.
  /// Returns the values for each basis function in the passed Array1D array
  void Get1dNormsSqExact(Array1D<double>& psi1dSqExact) const {psi1dSqExact=psi1dSqExact_; return;}

  /// \brief Get samples of the random variables associated
  /// with the current PC basis functions and return them in the 1D array randSamples.
  /// Take as many samples as the length of the array randSamples
  /// \note This function does NOT reset the random number seed before sampling
  void GetRandSample(Array1D<double>& randSamples);

  /// \brief Get nSamp samples of the random variables associated
  /// with the current PC basis functions and return them in the double* randSamples.
  /// \note This function does NOT reset the random number seed before sampling
  void GetRandSample(double* randSamples, const int& nSamp);

  /// \brief Get the random number generator seed
  int GetSeed() const {return rSeed_;}

  /// \brief Function to (re)seed the random number generator
  /// used to sample the Basis functions
  void SeedRandNumGen(const int& seed);

  /// \brief Get the quadrature integration information
  void GetQuadRule(Array2D<double>& qPoints, Array1D<double>& qWeights, Array2D<int>& qIndices);

  /// \brief Get the quadrature points in the passed Array2D array
  /// \note Although quadPoints is a 2D array, its second dimension is equal to 1
  void GetQuadPoints(Array2D<double>& quadPoints) const { quadPoints=quadPoints_; return;}

  /// \brief Get the quadrature weights in the passed Array1D array
  void GetQuadWeights(Array1D<double>& quadWeights) const { quadWeights=quadWeights_; return;}

  /// \brief Get the quadrature points' indices in the passed Array1D array
  void GetQuadIndices(Array2D<int>& quadIndices) const { quadIndices=quadIndices_; return;}

  /// \brief Get the basis values at quadrature points in the passed Array2D array
  void GetBasisAtQuadPoints(Array2D<double>& psi1d) const { psi1d=psi1d_; return;}

  /// \brief Get the PC type
  string GetPCType() const {return type_;}

  /// \brief Get the value of the parameter alpha
  double GetAlpha() const {return alpha_;}

  /// \brief Get the value of the parameter beta
  double GetBeta() const {return beta_;}

private:
  /// \brief Dummy default constructor, which should not be used as it is not well defined
  /// Therefore we make it private so it is not accessible
  /// \note All parameters are intialized to dummy values.
  // PCBasis(): type_("NA") {};

  /// \brief Dummy copy constructor, which should not be used as it is currently
  /// not well defined. Therefore we make it private so it is not accessible.
  /// \note I am not sure actually whether the initialization performed below
  /// is legal as it requires access to private data members of the class that
  /// is passed in.
  PCBasis(const PCBasis &obj):type_(obj.type_) {};


  /// \brief Evaluate the norms (squared) of the basis functions
  /// and stores in the private array psi1dSq_
  void Eval1dNormSq(int kord);


    /// \brief Evaluate 1d norm of order kord exactly
  double NormSq_Exact(int kord);

  /// \brief String indicator of type of basis functions used
  string type_;

  /// \brief Array to store quadrature points
  Array2D<double> quadPoints_;

  /// \brief Array to store quadrature weights
  Array1D<double> quadWeights_;

 /// \brief Array to store quadrature point indexing; useful only for nested rules
  Array2D<int> quadIndices_;


  /// \brief Array to store basis functions evaluated at quadrature points
  /// for each order: psi1d_(iqp,iord) contains the value of the polynomial
  /// chaos basis of order iord at the location of quadrature point iqp.
  Array2D<double> psi1d_;

  /// \brief Array with the norms squared of the 1D basis functions for each order
  Array1D<double> psi1dSq_;

  /// \brief Array with the exact norms squared of the 1D basis functions for each order
  Array1D<double> psi1dSqExact_;

  /// \brief Maximal order of any dimension
  int maxord_;

  /// \brief Number of parameters to specify the basis
  int narg_;

  /// \brief Parameter alpha for PCs that require a parameter (LG,SW,JB)
  double alpha_;

  /// \brief Parameter beta for PCs that require two parameters (SW,JB)
  double beta_;

  /// \brief Random sequence state for dsfmt
  /// \todo need more functionalities to get/set this variable from user
  dsfmt_t  rnstate_ ;

  /// \brief The seed used for the random number generators that
  /// sample the xi's in the basis functions
  ///
  /// This seed is set to 1 during the class construction and
  /// can be reset with the SeedRandNumGen function
  /// \sa SeedRandNumGen
  int rSeed_;

};

#endif /* PCBASIS_H_SEEN */
