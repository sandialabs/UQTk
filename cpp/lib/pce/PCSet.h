/* =====================================================================================
                     The UQ Toolkit (UQTk) version @UQTKVERSION@
                     Copyright (@UQTKYEAR@) Sandia Corporation
                     http://www.sandia.gov/UQToolkit/

     Copyright (@UQTKYEAR@) Sandia Corporation. Under the terms of Contract DE-AC04-94AL85000
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
/// \file PCSet.h
/// \author B. Debusschere, C. Safta, K. Sargsyan, K. Chowdhary 2007 -
/// \brief Header file for the Multivariate PC class

#ifndef PCSET_H_SEEN
#define PCSET_H_SEEN

#include <iostream>
#include <sstream>
#include <string.h>
#include <stdio.h>
#include <map>
#include "Array1D.h"
#include "Array2D.h"
#include "error_handlers.h"
#include "ftndefs.h"
#include "quad.h"

/*  CVODE headers  */
#include <cvode/cvode.h>             /* prototypes for CVODE fcts., consts.  */
#include <nvector/nvector_serial.h>  /* serial N_Vector types, fcts., macros */
#include <cvode/cvode_dense.h>       /* prototype for CVDense                */
#include <sundials/sundials_dense.h> /* definitions DlsMat DENSE_ELEM        */
#include <sundials/sundials_types.h> /* definition of type realtype          */

#include <iostream>
#include <string.h>
#include <stdio.h>
#include <sstream>
using namespace std; // needed for python string conversion

class PCBasis;
class Quad;

typedef enum {TaylorSeries=0, Integration} LogCompMethod;

/// \class  PCSet
/// \brief  Defines and initializes PC basis function set and provides functions
///         to manipulate PC expansions defined on this basis set.

class PCSet {
public:

  /// \brief Constructor: initializes the PC basis set for the
  /// order, number of dimensions and type that are passed in.
  ///
  /// Implementation type sp_type has three options "ISP" (intrusive methods), "NISP" (non-intrusive),
  /// or "NISPnoq" (non-intrusive without quadrature initialization)
  /// \note alpha and betta are parameters only relevant for GLG, JB or SW chaoses
  PCSet(const string sp_type, const int order, const int n_dim, const string pc_type, 
        const double alpha=0.0, const double betta=1.0);


  /// \brief Constructor: initializes the PC basis set for the
  /// order, number of dimensions and type that are passed in. It also
  /// customizes the multiindex sequence ( lexicographical-lex,
  /// colexicographical-colex, reverse lexicographical-revlex, and
  /// reverse clexicographical-revcolex).
  ///
  /// Implementation type sp_type has three options "ISP" (intrusive methods), "NISP" (non-intrusive),
  /// or "NISPnoq" (non-intrusive without quadrature initialization)
  /// \note alpha and betta are parameters only relevant for GLG, JB or SW chaoses
  PCSet(const string sp_type, const int order, const int n_dim, const string pc_type, const string pc_seq,
        const double alpha=0.0, const double betta=1.0);


  /// \brief Constructor: initializes the PC basis set ordered in an HDMR fashion
  /// given order per each HDMR rank (univariate, bivariate, etc...)
  ///
  /// Implementation type sp_type has three options "ISP" (intrusive methods), "NISP" (non-intrusive),
  /// or "NISPnoq" (non-intrusive without quadrature initialization)
  /// \note alpha and betta are parameters only relevant for GLG, JB or SW chaoses
  PCSet(const string sp_type, const Array1D<int>& maxOrders, const int n_dim, const string pc_type, 
        const double alpha=0.0, const double betta=1.0);

  /// \brief Constructor: initializes the PC basis set for a given
  /// custom multiIndex
  ///
  /// Implementation type sp_type has three options
  /// "ISP" (intrusive methods), "NISP" (non-intrusive),
  /// or "NISPnoq" (non-intrusive without quadrature initialization)
  /// \note alpha and betta are parameters only relevant for GLG, JB or SW chaoses
  PCSet(const string sp_type, const Array2D<int>& customMultiIndex, const string pc_type, 
        const double alpha=0.0, const double betta=1.0);

  /// \brief Destructor: cleans up all memory and destroys object
  ~PCSet();

 /////////////////////////////////////////////////////////////////////////////////////////
  /// Set Gradient and Hessian operators
  /////////////////////////////////////////////////////////////////////////////////////////

  /// \brief Evaluate Gradient at a single d-dim point
  /// and d-dim basis polynomial
  void dPhi_alpha(Array1D<double>& x, Array1D<int>& alpha, Array1D<double>& grad);
   /// \brief Evaluate Gradient at a single d-dim point
  /// for a PCSet object
  void dPhi(Array1D<double>& x, Array2D<int>& mindex, Array1D<double>& grad, Array1D<double>& ck);
     /// \brief Evaluate Gradient at a multiple d-dim x points
  /// for a PCSet object
  void dPhi(Array2D<double>& x, Array2D<int>& mindex, Array2D<double>& grad, Array1D<double>& ck);

  /// \brief Evaluate Hessian at a single d-dim point
  /// and d-dim basis polynomial
  void ddPhi_alpha(Array1D<double>& x, Array1D<int>& alpha, Array2D<double>& grad);
   /// \brief Evaluate Gradient at a single d-dim point
  /// for a PCSet object
  void ddPhi(Array1D<double>& x, Array2D<int>& mindex, Array2D<double>& grad, Array1D<double>& ck);

  /////////////////////////////////////////////////////////////////////////////////////////
  /// Set the quadrature rule
  /////////////////////////////////////////////////////////////////////////////////////////

  /// \brief Obtain 1d quadrature points and weights 
  /// \note This is used in triple or quadruple product computation for which the default quadrature is not enough
  void SetQd1d(Array1D<double>& qdpts1d,Array1D<double>& wghts1d, int nqd);

  /// \brief Set the quadrature points by specifying
  /// a grid type, a full/sparse indicator, and an integer parameter
  ///
  /// Full/sparse switch fs_type can be either 'full' or 'sparse'
  /// The parameter param is the number of points per dimension for
  /// full quadrature, and the level for sparse quadrature
  /// Options for grid_type are, besides the standard PC types,
  /// 'CC' (Clenshaw-Curtis), 'CCO' (Clenshaw-Curtis open),
  /// 'NC' (Newton-Cotes), 'NCO' (Newton-Cotes open),
  /// where open means that endpoints are expluded
  /// \note 'NC', 'NCO' quadratures are the same as uniformly spaced grids
  void SetQuadRule(const string grid_type,const string fs_type,int param);

  /// \brief Set a custom quadrature rule by pointing to the corresponding object
  void SetQuadRule(Quad &quadRule);

  /////////////////////////////////////////////////////////////////////////////////////////
  /// Print information on the screen
  /////////////////////////////////////////////////////////////////////////////////////////

  /// \brief Print the multi-indices for all terms on the screen
  void PrintMultiIndex() const;

  /// \brief For all terms, print their multi-index and norm^2 on the screen
  void PrintMultiIndexNormSquared() const;


  /////////////////////////////////////////////////////////////////////////////////////////
  /// Get and set variables/arrays inline
  /////////////////////////////////////////////////////////////////////////////////////////

  /// \brief Get the PC type
  string GetPCType() const {return pcType_;}

  /// \brief Get the value of the parameter alpha
  double GetAlpha() const {return alpha_;}

  /// \brief Get the value of the parameter beta
  double GetBeta() const {return beta_;}

  /// \brief Get the multiindex (return Array2D)
  void GetMultiIndex(Array2D<int> &mindex) const {mindex=multiIndex_;}

  /// \brief Get the multiindex (return double *)
  void GetMultiIndex(int *mindex) const;

  /// \brief Get the norm-squared
  /// \todo this seems like a duplication, see below GetPsiSq()
  void GetNormSq(Array1D<double>& normsq) const {normsq=psiSq_;}

  /// \brief Get the number of terms in a PC expansion of this order and dimension
  int GetNumberPCTerms() const {return nPCTerms_;}

  /// \brief Get the PC dimensionality
  int GetNDim() const {return nDim_;}

  /// \brief Get the PC order
  int GetOrder() const {return order_;}

  /// \brief Get the number of quadrature points
  int GetNQuadPoints() const {return nQuadPoints_;}

  /// \brief Get the quadrature points
  void GetQuadPoints(Array2D<double>& quad) const {quad=quadPoints_;}

  /// \brief Get the quadrature points and weights
  void GetQuadPointsWeights(Array2D<double>& quad, Array1D<double>& wghts) const { quad=quadPoints_; wghts=quadWeights_;}

  /// \brief Get the quadrature points folded into a one-dimensional array quad
  void GetQuadPoints(double* quad) const {for(int i=0;i<nQuadPoints_;i++) for(int j=0;j<nDim_;j++) quad[i*nDim_+j]=quadPoints_(i,j);}

  /// \brief Get the quadrature weights
  void GetQuadWeights(Array1D<double>& wghts) const {wghts=quadWeights_;}

  /// \brief Get the quadrature weights folded into a one-dimensional array wghts
  void GetQuadWeights(double* wghts) const {for(int i=0;i<nQuadPoints_;i++) wghts[i]=quadWeights_(i);}

  /// \brief Get the values of the
  /// basis polynomials evaluated at the quadrature points
  void GetPsi(Array2D<double>& psi) const {psi=psi_;}

  /// \brief Get the polynomials evaluated at the quadrature points
  /// folded into a one-dimensional array psi
  void GetPsi(double* psi) const  {for(int i=0;i<nQuadPoints_;i++) for(int j=0;j<nPCTerms_;j++) psi[i*nPCTerms_+j]=psi_(i,j);}

  /// \brief Get the basis polynomial norms-squared in an array class object psisq
  void GetPsiSq(Array1D<double>& psisq) const {psisq=psiSq_;}

  /// \brief Get the basis polynomial norms-squared in a double* array psisq
  void GetPsiSq(double* psisq) const {for(int i=0;i<nPCTerms_;i++) psisq[i]=psiSq_(i);}

  /// \brief Get relative tolerance for Taylor series approximations
  double GetTaylorTolerance() const {return rTolTaylor_;}

  /// \brief Set relative tolerance for Taylor series approximations
  void SetTaylorTolerance(const double& rTol) {rTolTaylor_ = rTol;}

  /// \brief Get maximum number of terms in Taylor series approximations
  int GetTaylorTermsMax() const {return maxTermTaylor_;}

  /// \brief Set maximum number of terms in Taylor series approximations
  void SetTaylorTermsMax(const int& maxTerm) {maxTermTaylor_ = maxTerm;}

  /// \brief Set method of computing the log function.
  ///
  /// Use the argument TaylorSeries to select the Taylor series approach or
  /// Integration to select the integration method.
  void SetLogCompMethod(const LogCompMethod& logMethod) {logMethod_ = logMethod;}

  /// \brief Get relative tolerance for GMRES in Div routine
  double GetGMRESDivTolerance() const {return rTolGMRESDiv_;}

  /// \brief Set the relative tolerance for GMRES in Div routine
  void SetGMRESDivTolerance(const double& rTol) {rTolGMRESDiv_ = rTol;}


  /////////////////////////////////////////////////////////////////////////////////////////
  /// Intrusive arithmetics
  /////////////////////////////////////////////////////////////////////////////////////////

  /// \brief Initializes a PC expansion p in a double* format to have the same distribution as the underlying PC germ,
  /// but with a specified mean m and standard deviation s
  /// \note This assumes that the zeroth order term is the first one in the multi-index - this assumption does not hold in general
  /// \note This function only holds for expansions with one stochastic dimension
  /// \note All existing coefficient values in p will be overwritten
  /// \todo Make this function work for general multi-indices, and for any number of stochastic dimensions
  void InitMeanStDv(const double& m, const double& s, double* p) const;

  /// \brief Initializes a PC expansion p in Array1D<double> format to have the same distribution as the underlying PC germ,
  /// but with a specified mean m and standard deviation s
  /// \note This assumes that the zeroth order term is the first one in the multi-index - this assumption does not hold in general
  /// \note This function only holds for expansions with one stochastic dimension
  /// \note All existing coefficient values in p will be overwritten
  /// \todo Make this function work for general multi-indices, and for any number of stochastic dimensions
  void InitMeanStDv(const double& m, const double& s, Array1D<double>& p) const;

  /// \brief Copy PC expansion p2 into p1 (i.e. p1 = p2).
  ///
  /// All arguments in double* format.
  void Copy(double* p1, const double* p2) const;

  /// \brief Copy PC expansion p2 into p1 (i.e. p1 = p2).
  ///
  /// All arguments in Array format
  /// \note Requires the size of the arrays that are passed in to equal the number of PC terms
  void Copy(Array1D<double>& p1, const Array1D<double>& p2) const;

  /// \brief Add two PC expansions given by double* arguments p1 and p2, and
  /// return the result in p3.
  void Add(const double* p1, const double* p2, double* p3) const;

  /// \brief Add two PC expansions given by Array1D arguments p1 and p2, and
  /// return the result in p3.
  /// \note Requires the size of the arrays that are passed in to equal the number of PC terms
  void Add(const Array1D<double>& p1, const Array1D<double>& p2, Array1D<double>& p3) const;

  /// \brief Add PC expansions given by double* argument p2 to p1 and
  /// return the result in p1.
  void AddInPlace(double* p1, const double* p2) const;

  /// \brief Add PC expansions given by Array1D argument p2 to p1 and
  /// return the result in p1.
  /// \note Requires the size of the arrays that are passed in to equal the number of PC terms
  void AddInPlace(Array1D<double>& p1, const Array1D<double>& p2) const;

  /// \brief Multiply PC expansion p1 with scalar a and return the result in p2.
  /// All PCEs are in double* format
  void Multiply(const double* p1, const double& a, double* p2) const;

  /// \brief Multiply PC expansion p1 with scalar a and return the result in p2.
  /// All PCEs are in Array1D format
  /// \note Requires the size of the arrays that are passed in to equal the number of PC terms
  void Multiply(const Array1D<double>& p1, const double& a, Array1D<double>& p2) const;

  /// \brief Multiply PC expansions given by double* argument p1 with scalar a and
  /// return the result in p1.
  void MultiplyInPlace(double* p1, const double& a) const;

  /// \brief Multiply PC expansions given by Array1D argument p1 with scalar a and
  /// return the result in p1.
  /// \note Requires the size of the arrays that are passed in to equal the number of PC terms
  void MultiplyInPlace(Array1D<double>& p1, const double& a) const;

  /// \brief Subtract PC expansion p2 from p1, and return the result in p3, with
  /// all arguments given as double*
  void Subtract(const double* p1, const double* p2, double* p3) const;

  /// \brief Subtract PC expansion p2 from p1, and return the result in p3, with
  /// all arguments given as Array1D structures.
  /// \note Requires the size of the arrays that are passed in to equal the number of PC terms
  void Subtract(const Array1D<double>& p1, const Array1D<double>& p2, Array1D<double>& p3) const;

  /// \brief Subtract PC expansion p2 from p1, and return the result in p1, with
  /// all arguments given as double*
  void SubtractInPlace(double* p1, const double* p2) const;

  /// \brief Subtract PC expansion p2 from p1, and return the result in p1, with
  /// all arguments given as Array1D structures.
  /// \note Requires the size of the arrays that are passed in to equal the number of PC terms
  void SubtractInPlace(Array1D<double>& p1, const Array1D<double>& p2) const;

  /// \brief Multiply two PC expansions given by double* arguments p1 and p2, and
  /// return the result in p3.
  void Prod(const double* p1, const double* p2, double* p3) const;

  /// \brief Multipy two PC expansions given by Array1D arguments p1 and p2, and
  /// return the result in p3.
  /// \note Requires the size of the input arrays to equal the number of PC terms
  void Prod(const Array1D<double>& p1, const Array1D<double>& p2, Array1D<double>& p3) const;

  /// \brief Multiply three PC expansions given by double* arguments
  /// p1, p2, and p3, and return the result in p4.
  void Prod3(const double* p1, const double* p2, const double* p3, double* p4) const;

  /// \brief Multipy three PC expansions given by Array1D arguments p1,
  /// p2, and p3, and return the result in p4.
  /// \note Requires the size of the input arrays to equal the number of PC terms
  void Prod3(const Array1D<double>& p1, const Array1D<double>& p2, const Array1D<double>& p3,
	    Array1D<double>& p4) const;

  /// \brief Evaluates a polynomial of PC that is given in double* argument p1.
  /// Polynomial coefficients are given in double* argument polycf of size npoly.
  /// The output PC is contained in double* argument p2.
  /// \note Recursive algorithm is implemented.
  void Polyn(const double* polycf, int npoly, const double* p1, double* p2) const;

  /// \brief Evaluates a polynomial of PC that is given by Array1D argument p1.
  /// Polynomial coefficients are given in the Array1D argument polycf.
  /// The output PC is contained in Array1D argument p2.
  /// \note Requires the size of array p1 to equal the number of PC terms
  void Polyn(const Array1D<double>& polycf, const Array1D<double>& p1, Array1D<double>& p2) const;

  /// \brief Evaluates a multivariate polynomial of a set of PC inputs given by
  /// Array2D argument p1 (each column of p1 is a PC input).
  /// Polynomial coefficients are given in Array1D argument polycf.
  /// Multiindex set for the multivariate polynomial is given in Array2D argument mindex.
  /// The output PC is contained in Array1D argument p2.
  /// \note Requires the size of the array polycf to equal the first dimension of argument mindex
  /// \note Requires the size of the array p1 to equal (the number of PC terms) X (second dimension of argument mindex)
  /// \note Uses a recursive algorithm
  /// \note Out of convenience, this function so far is implemented for Array classes, not double* arrays.
  /// \todo A double* version should be added.
  void PolynMulti(const Array1D<double>& polycf, const Array2D<int>& mindex, const Array2D<double>& p1, Array1D<double>& p2) const;

  /// \brief Take the exp() of the PC expansion given by double* argument p1, and
  /// return the result in p2.
  ///
  /// Relies on Taylor series expansion: exp(x) = 1 + x + x^2/2! + x^3/3! + ...
  /// However, for efficiency and to avoid overflow, the terms are computed
  /// as d_i = d_{i-1}*x/i. Also, to reduce the number of terms needed in the
  /// series, we subtract the mean out of a random variable u as u = u_0 + (u-u_0)
  /// and exp(u) = exp(u_0)*exp(u-u_0), where exp(u_0) can be computed with the
  /// regular exp(double& ) function
  /// \note The Taylor series is truncated after a tolerance criterium is
  /// achieved on the relative error defined as the max absolute value of the PC
  /// coefficients in the last added term, divided by the mean of exp(p1).
  /// The tolerance is set to 1.e-6 by default and can be changed with SetTaylorTolerance().
  /// \note The maximum number of terms in the Taylor series is set by default to 500
  /// and can be changed with SetTaylorTermsMax()
  void Exp(const double* p1, double* p2) const;

  /// \brief Take the exp() of the PC expansion given by Array1D argument p1, and
  /// return the result in p2.
  /// \note Requires the size of the arrays that are passed in to equal the number of PC terms
  void Exp(const Array1D<double>& p1, Array1D<double>& p2) const;

  /// \brief Take the natural logarithm log() of the PC expansion given by double*
  /// argument p1, and return the result in p2. The logarithm is evaluated
  /// either via Taylor series or via integration depending on the value of
  /// parameter logMethod_
  void Log(const double* p1, double* p2) const;

  /// \brief Take the natural logarithm, log(), of the PC expansion given by
  /// Array1D argument p1, and return the result in Array1D argument p2.
  /// \note Requires the size of the arrays that are passed in to equal the number of PC terms
  void Log(const Array1D<double>& p1, Array1D<double>& p2) const;

  /// \brief Take the logarithm to base 10 of the PC expansion given by double*
  /// argument p1, and return the result in p2.
  ///
  /// First use Log() to compute the natural logarithm and then divide it by log(10)
  void Log10(const double* p1, double* p2) const;

  /// \brief Take the logarithm to base 10 of the PC expansion given by
  /// Array1D argument p1, and return the result in Array1D argument p2.
  /// \note Requires the size of the arrays that are passed in to equal the number of PC terms
  void Log10(const Array1D<double>& p1, Array1D<double>& p2) const;

  /// \brief Evaluate power a (a real number) of PC expansion given by double*
  /// argument p1, and return the result in p2.
  /// The power is computed as p1^a = exp(a*log(p1)), where log(p1) is evaluated
  /// either via Taylor series or via integration depending on the value of
  /// parameter logMethod_
  void RPow(const double* p1, double* p2, const double& a) const;

  /// \brief Evaluate power a (a real number) of PC expansion given by
  /// Array1D argument p1, and return the result in Array1D argument p2.
  /// \note Requires the size of the arrays that are passed in to equal the number of PC terms
  void RPow(const Array1D<double>& p1, Array1D<double>& p2, const double& a) const;

  /// \brief Evaluate power ia (an integer number) of PC expansion given by double*
  /// argument p1, and return the result in p2.
  void IPow(const double* p1, double* p2, const int& ia) const;

  /// \brief Evaluate power ia (an integer number) of PC expansion given by
  /// Array1D argument p1, and return the result in Array1D argument p2.
  /// \note Requires the size of the arrays that are passed in to equal the number of PC terms
  void IPow(const Array1D<double>& p1, Array1D<double>& p2, const int& ia) const;

  /// \brief Evaluate the inverse of PC expansion given by double*
  /// argument p1, and return the result in p2.
  /// The inverse is computed using the division function
  void Inv(const double* p1, double* p2) const;

  /// \brief Evaluate the inverse of PC expansion given by
  /// Array1D argument p1, and return the result in Array1D argument p2.
  /// \note Requires the size of the arrays that are passed in to equal the number of PC terms
  void Inv(const Array1D<double>& p1, Array1D<double>& p2) const;

  /// \brief Divide the PC expansion p1 by p2, and return the result in p3
  /// (All arguments in double* format)
  ///
  /// The "division" p3 = p1/p2 is performed by solving the system
  /// of equations p2*p3 = p1 for the unknown p3.
  /// \note When GMRES is used to solve this system of equations (based on a
  /// preprocessor flag in the source code for this routine), a relative tolerance
  /// criterium is used that is set by default to 1.e-8, and can be
  /// changed with SetGMRESDivTolerance().
  /// \todo Remove duplication of data and parameters that was required for
  /// enforcing imposed "const" constraints on some of the arguments and the class data
  /// members when they are being passed to fortran.
  void Div(const double* p1, const double* p2, double* p3) const;

  /// \brief Divide the PC expansion p1 by p2, and return the result in p3
  /// (All arguments in Array1D<double> format)
  /// \note Requires the size of the arrays that are passed in to equal the number of PC terms
  void Div(const Array1D<double>& p1, const Array1D<double>& p2, Array1D<double>& p3) const;

  /// \brief Returns the standard deviation of PC expansion p in a double* format
  /// \note This assumes that the zeroth order term is the first one in the multi-index - this assumption does not hold in general
  /// \todo Lift the assumption by looking for the constant term in the multiindex
  double StDv(const double* p) const;

  /// \brief Returns the standard deviation of PC expansion p
  ///(Argument in Array1D<double> format)
  /// \note This assumes that the zeroth order term is the first one in the multi-index - this assumption does not hold in general
  /// \todo Lift the assumption by looking for the constant term in the multiindex
  /// \note For a more general implementation, see ComputeVarFrac()
  double StDv(const Array1D<double>& p) const;

  /// \brief Compute the rms average of the PC coefficients (i.e. the square root
  /// of the average of the square of the PC coefficients, not taking into
  /// account any basis functions). (Arguments in double* format)
  double GetModesRMS(const double* p) const;

  /// \brief Compute the rms average of the PC coefficients (i.e. the square root
  /// of the average of the square of the PC coefficients, not taking into
  /// account any basis functions). (Arguments in Array1D<double> format)
  /// \note Requires the size of the array that is passed in to equal the number of PC terms
  double GetModesRMS(const Array1D<double>& p) const;

  /// \brief Computes derivatives of univariate PC given by coefficients p1
  /// returns coefficient vector of the derivative in p2
  /// \note Makes use of intrusive computations on recursive formulae for derivatives
  /// \todo Supports LU and HG bases only
  /// \todo Supports only for 1d PCs
  void Derivative(const double* p1, double* p2) const;

  /// \brief Computes derivatives of univariate PC given by coefficients p1
  /// returns coefficient vector of the derivative in p2
  /// \note Makes use of intrusive computations on recursive formulae for derivatives
  /// \todo Supports LU and HG bases only
  /// \todo Supports only for 1d PCs
  void Derivative(const Array1D<double>& p1, Array1D<double>& p2) const;

  /// \brief Returns number of triple products
  int GetNumTripleProd() const;
  /// \brief Returns triple products indices (int*/double* version)
  void GetTripleProd(int *nTriple, int *iProd, int *jProd, double *Cijk) const;
  /// \brief Returns triple products indices (Array version)
  void GetTripleProd(Array1D<int>& nTriple, Array1D<int>& iProd, Array1D<int>& jProd, Array1D<double>& Cijk) const;
  /// \brief Returns number of quad products
  int  GetNumQuadProd() const;
  /// \brief Returns quad products indices (int*/double* version)
  void GetQuadProd(int *nQuad, int *iProd, int *jProd, int *kProd, double *Cijkl) const;
  /// \brief Returns quad products indices (Array version)
  void GetQuadProd(Array1D<int> &nQuad, Array1D<int> &iProd, Array1D<int> &jProd, Array1D<int> &kProd, 
		   Array1D<double> &Cijkl) const;

  /////////////////////////////////////////////////////////////////////////////////////////
  /// Random sample generator functions
  /////////////////////////////////////////////////////////////////////////////////////////

  /// \brief Reseed the random number generator used for the sampling
  /// of the PC variables and expansions
  void SeedBasisRandNumGen(const int& seed) const;

  /// \brief Draw a set of samples from the PC expansion p,
  /// and return the result in the array samples.
  /// All arguments are in Array1D<double> format
  /// The number of samples requested is assumed to be the size of the samples array
  /// \note The size of the array p that is passed in needs to equal the number of PC terms
  void DrawSampleSet(const Array1D<double>& p, Array1D<double>& samples);

  /// \brief Draw a set of samples from the PC expansion given in double* argument p,
  /// and return the result in double* array samples.
  /// The number of samples requested is the argument nSamples
  void DrawSampleSet(const double* p, double* samples, const int& nSamples);

  /// \brief Draw a set of samples of the underlying germ random variable
  /// \todo There is no double* version of this function
  void DrawSampleVar(Array2D<double>& samples) const;
  void DrawSampleVar(double *samples, const int &nS, const int &nD) const;

  /////////////////////////////////////////////////////////////////////////////////////////
  /// PC evaluation functionalities
  /////////////////////////////////////////////////////////////////////////////////////////

  /// \brief Evaluate the given PC expansion p, at the specified
  /// values of the random variables, randVarSamples.
  /// All arguments in const Array1D<double> format
  /// \note The number of elements in p needs to match the number of terms
  /// in the PC expansions in this PCSet.
  /// \note The number of elements in randVarSamples needs to match the
  /// number of dimensions in the PC expansion.
  double EvalPC(const Array1D<double>& p, Array1D<double>& randVarSamples);

  /// \brief Evaluate the given PC expansion p, at the specified
  /// values of the random variables, randVarSamples.
  /// All arguments in const double* format
  /// \note The number of elements in p is assumed to match the number of terms
  /// in the PC expansions in this PCSet.
  /// \note The number of elements in randVarSamples is assumed to match the
  /// number of dimensions in the PC expansion.
  double EvalPC(const double* p, const double* randVarSamples);

  /// \brief Evaluate the given PC expansion at given set of points with given coefficient vector and
  /// return the values in an 1D Array in the first argument.
  /// \todo There is no double* version of this function
  void EvalPCAtCustPoints(Array1D<double>& xch, Array2D<double>& custPoints,Array1D<double>& p);

  /// \brief Evaluate Basis Functions at given points custPoints and return in the array psi
  /// \todo There is no double* version of this function
  void EvalBasisAtCustPts(const Array2D<double>& custPoints,Array2D<double>& psi);

  void EvalBasisAtCustPts(const double* custPoints, const int npts, double* psi);
  // void EvalBasisAtCustPts(const int npts, const int ndim, const int npc, const double *custPoints, double *psi);

  /////////////////////////////////////////////////////////////////////////////////////////
  /// Galerkin projection functionalities
  /////////////////////////////////////////////////////////////////////////////////////////

  /// \brief Performs (NISP) Galerkin projection, given function evaluations at quadrature points
  /// Returns in the coefficient vector in the second argument
  /// \note User should make sure that the function HAS BEEN evaluated at the correct quadrature points
  /// by first extracting the quadrature points and evaluating the function externally
  /// \todo Overload this with forward function pointers
  /// \todo There is no double* version of this function
  void GalerkProjection(const Array1D<double>& fcn, Array1D<double>& ck);

  /// \brief Galerkin Projection via Monte-Carlo integration
  /// \note User should make sure that the function HAS BEEN evaluated at the correct sampling points
  /// by first sampling the proper PC germ distribution and evaluating the function externally
  /// \todo Overload this with forward function pointers
  /// \todo There is no double* version of this function
  void GalerkProjectionMC(const Array2D<double>& x, const Array1D<double>& fcn, Array1D<double>& ck);

  /////////////////////////////////////////////////////////////////////////////////////////
  /// Multiindex parsing functionalities
  /////////////////////////////////////////////////////////////////////////////////////////

  /// \brief Computes the order of each basis term and return it in the array orders,
  /// also returns the maximal order
  /// \todo There is no double* version of this function
  int ComputeOrders(Array1D<int>& orders);

  /// \brief Computes the effective dimensionality of each basis term,
  /// i.e., the number of dimensions that enter with a non-zero degree.
  /// also returns the maximal dimensionality among all basis terms
  /// \note This is not the classical effective dimensionality,
  /// since all dimensions can still be involved.
  int ComputeEffDims(int *effdim);

  /// \brief Computes the effective dimensionality of each basis term,
  /// i.e., the number of dimensions that enter with a non-zero degree.
  /// also returns the maximal dimensionality among all basis terms
  /// \note This is not the classical effective dimensionality,
  /// since all dimensions can still be involved.
  int ComputeEffDims(Array1D<int> &effdim);

  /// \brief Encode multiIndex into a 'sparse' format where the bases are ordered by their effective dimensionality.
  /// The i-th element in sp_mindex stores all the bases that have effective dimensionality  equal to i.  Also, only non-zero components are stored.
  /// \todo There is no double* version of this function
  void EncodeMindex(Array1D< Array2D<int> >& sp_mindex);

  /////////////////////////////////////////////////////////////////////////////////////////
  /// Moment/sensitivity extraction given coefficients
  /////////////////////////////////////////////////////////////////////////////////////////

  /// \brief Compute the mean of the PC given coefficients in double *coef 
  /// (seeking the zero-th order multiindex)
  double ComputeMean(const double *coef);

  /// \brief Compute the mean of the PC given coefficient array coef(seeking the zero-th order multiindex)
  double ComputeMean(Array1D<double>& coef);

  /// \brief Compute the variance fractions of each basis term given
  /// coefficients in double *coef; returns the variance fractions in the double *varfrac
  /// \note Also returns the variance
  /// \note The value for the zeroth order term has a special meaning: it is equal to mean^2/variance or (mean/std)^2.
  double ComputeVarFrac(const double *coef, double *varfrac);

  /// \brief Compute the variance fractions of each basis term given
  /// coefficient array coef; returns the variance fractions in the array varfrac
  /// \note Also returns the variance
  /// \note The value for the zeroth order term has a special meaning: it is equal to mean^2/variance or (mean/std)^2.
  double ComputeVarFrac(Array1D<double>& coef, Array1D<double>& varfrac);

  /// \brief Compute main effect sensitivity (Sobol) indices given coefficient array coef; returns the indices in the array mainsens
  /// \todo There is no double* version of this function
  void ComputeMainSens(Array1D<double>& coef, Array1D<double>& mainsens);

  /// \brief Compute total effect sensitivity (Sobol) indices given coefficient array coeff; returns the indices in the array totsens
  /// \todo There is no double* version of this function
  void ComputeTotSens(Array1D<double>& coef, Array1D<double>& totsens);

  /// \brief Compute joint effect sensitivity (Sobol) indices given coefficient array coeff; returns the indices in the array jointsens
  /// \note jointsens will be populated as a strictly upper-diagonal matrix
  /// \todo There is no double* version of this function
  void ComputeJointSens(Array1D<double>& coef, Array2D<double>& jointsens);


  /////////////////////////////////////////////////////////////////////////////////////////
  /// Other
  /////////////////////////////////////////////////////////////////////////////////////////

  /// \brief Set the verbosity level
  /// \note Currently, the values of 0, 1 and 2 are implemented
  void SetVerbosity(int verbosity) { uqtkverbose_ = verbosity; }

  /// \brief Evaluate norms-squared of all bases and return in the array normsq
  /// \todo There is no double* version of this function
  void EvalNormSq(Array1D<double>& normsq);
  void EvalNormSq(double* normsq, const int npc);

  /// \brief Evaluate norms-squared analytically of all bases and return in the array normsq
  /// \todo There is no double* version of this function
  /// \note Custom PCs do not have this capability
  void EvalNormSqExact(Array1D<double>& normsq);

  /// \brief Check if the point x is in the PC domain
  bool IsInDomain(double x);



 private:
  /// \brief Dummy default constructor, which should not be used as it is not well defined
  /// Therefore we make it private so it is not accessible
  /// \note All parameters are intialized to dummy values.
  PCSet(): order_(0), nDim_(0) {};


  /// \brief Dummy copy constructor, which should not be used as it is currently
  /// not well defined. Therefore we make it private so it is not accessible.
  /// \note I am not sure actually whether the initialization performed below
  /// is legal as it requires access to private data members of the class that
  /// is passed in.
  PCSet(const PCSet &obj):order_(obj.order_), nDim_(obj.nDim_) {};

  /// \brief Compute maximal order per dimension and fill in the array maxOrdPerDim_
  void ComputeMaxOrdPerDim();

  /// \brief Initialization of the appropriate variables
  /// \note Intrusive implementation only works with TotalOrder multiindes
  /// \todo Test and allow intrusive implementation with customized multiindices
  void Initialize(const string ordertype);

  /// \brief Initialize quadrature for computing triple products(ISP) and orthogonal projection(NISP)
  //  void InitQuadrature();

  /// \brief Initialize variables that are needed only in intrusive computations
  void InitISP();
  /// \brief Initialize variables that are needed only in non-intrusive computations
  void InitNISP();

  /// \brief Evaluate the expectation of product of three basis functions
  void EvalBasisProd3();

  /// \brief Evaluate the expectation of product of four basis functions
  void EvalBasisProd4();

  /// \brief Wrapper for Matrix-vector multiplication routine to be called by GMRES.
  ///
  /// As GMRES is a Fortran77 routine, this routine is
  /// defined as a static function. One of the function arguments (obj)
  /// was originally isym, a flag for matrix symmetry, but has been
  /// repurposed to carry an integer handle to identify this object.
  /// \note The matrix vector product here comes down to a product
  /// between two PC expansions.
  static void GMRESMatrixVectorProdWrapper(int* n, double* x, double* y, int* nelt,
                                    int* ia, int* ja, double* a, int* obj) {
    // Look up *obj in the map that relates integer indices to pointers to PCSet
    OMap_t::iterator it = omap_->find(*obj);
    if(it == omap_->end()) {
      string err_message = (string) "GMRESMatrixVectorProdWrapper():"
                           + " the callback object is not a valid entry in the map";
      throw Tantrum(err_message);
    }
    // Perform callback to the member function of the proper PCSet instance
    it->second->GMRESMatrixVectorProd(x, a, y);

    return;
  }

  /// \brief Wrapper for preconditioner routine to be called by GMRES.
  ///
  /// As GMRES is a Fortran77 routine, this routine is
  /// defined as a static function. One of the function arguments (obj)
  /// was originally isym, a flag for matrix symmetry, but has been
  /// repurposed to carry an integer handle to identify this object.
  /// \note Since we currently do not use preconditioning, this routine
  /// does nothing. It is a place holder for future use.
  static void GMRESPreCondWrapper(int* n, double* r, double* z, int* nelt,
                                  int* ia, int* ja, double* a, int* obj,
                                  double* rwork, int* iwork) { };

  /// \brief Actual C++ implementation of the matric vector multiplication
  /// for GMRES for the division operation.
  ///
  /// Given the structure of the problem, this boils down to the product
  /// between two PC variables.
  void GMRESMatrixVectorProd(const double* x, const double*a, double* y) const;

  /// \brief Computes natural logarithm using Taylor expansion:
  ///                                     N
  ///        p2 = ln(p1) =  ln(p1Mean) + sum d
  ///                                    n=1  n
  ///
  ///                       (n+1)
  ///                   (-1)       n            p1
  ///        where  d = ----     *x , and x = ------ - 1
  ///                n   n                    p1Mean
  ///
  /// \note See Exp notes for info related to tolerance and maximum number of terms criteria
  /// for truncating the Taylor series
  void LogTaylor(const double* p1, double* p2) const;

  /// \brief Computes natural logarithm by numerical integration:
  /// calculate p2=ln(p1) by integrating du=dx/x to get ln(x)
  void LogInt(const double* p1, double* p2) const;

  /// \brief Wrapper for LogIntRhs. The first component of f_data pointer
  /// carries an integer handle identifying the appropriate PC object
  /// \todo Why is this function a static int instead of static void? Should
  /// there be a return statement at the end?
  static int LogIntRhsWrapper(realtype t, N_Vector y, N_Vector ydot, void *f_data)
  {
    double indxobj = ((double*) f_data)[0] ;

    OMap_t::iterator it = omap_->find((int) indxobj);

    if (it == omap_->end())
    {
      string err_message = (string) "LogIntRhsWrapper():"
                           + " the callback object is not a valid entry in the map";
      throw Tantrum(err_message);
    }

    // Perform callback to the member function of the proper PCSet instance
    it->second->LogIntRhs(t,y,ydot,f_data);

    return ( 0 ) ;

  }

  /// \brief Evaluates rhs necessary to compute natural logarithm via integration
  int LogIntRhs(realtype t, N_Vector y, N_Vector ydot, void *f_data) const;

  /// \brief Verbosity level
  /// \note Currently the values of 0, 1 or 2 are implemented.
  int uqtkverbose_;

  /// \brief String indicator of ISP or NISP implementation type
  string spType_;

  /// \brief String indicator of PC type
  string pcType_;

  /// \brief String indicator of multiindex ordering 
  string pcSeq_;

  /// \brief Pointer to the class that defines the basis type and functions
  PCBasis* p_basis_;

  /// \brief Order of the PC representation
  int order_;

  /// \brief Maximal order within all dimensions
  int maxorddim_;

  /// \brief Array of maximum orders requested if custom(HDMR) ordering is requested
  Array1D<int> maxOrders_;

 /// \brief Array of maximum orders per dimension
  Array1D<int> maxOrdPerDim_;

  /// \brief Number of stochastic dimensions (degrees of freedom) in the PC representation
  const int nDim_;

  /// \brief Number of quadrature points used
  int nQuadPoints_;

  /// \brief Total number of terms in the PC expansions
  int nPCTerms_;

  /// \brief Relative tolerance for Taylor series approximations
  double rTolTaylor_;

  /// \brief Max number of terms in Taylor series approximations
  int maxTermTaylor_;

  /// \brief Tolerance to avoid floating-point errors
  double SMALL_;

  /// \brief GMRES tolerance in Div()
  double rTolGMRESDiv_;

  /// \brief Array to store basis functions evaluated at quadrature points
  /// for each order: psi_(iqp,ipc) contains the value of the polynomial
  /// chaos ipc-th basis at the location of quadrature point iqp.
  Array2D<double> psi_;

  /// \brief Array with the norms squared of the basis functions,
  /// corresponding to each term in the PC expansion
  Array1D<double> psiSq_;

  /// \brief Array to store quadrature points
  Array2D<double> quadPoints_;

  /// \brief Array to store quadrature weights
  Array1D<double> quadWeights_;

  /// \brief Array to store quadrature point indexing; useful for nested rules
  Array2D<int> quadIndices_;

  /// \brief Array to store multi-index: multiIndex_(ipc,idim) contains the order
  /// of the basis function associated with dimension idim, for the ipc-th term in the PC
  /// expansion.
  Array2D<int> multiIndex_;

  /// \brief i-indices of <\\Psi_i \\Psi_j \\Psi_k> terms that are not zero, for all k
  /// \note Stored as a vector over k, with each element being a vector of i-indices itself
  Array1D<Array1D<int> > iProd2_;

  /// \brief j-indices of <\\Psi_i \\Psi_j \\Psi_k> terms that are not zero, for all k
  /// \note Stored as a vector over k, with each element being a vector of j-indices itself
  Array1D<Array1D<int> > jProd2_;

  /// \brief <\\Psi_i \\Psi_j \\Psi_k> terms that are not zero, for all k
  /// \note Stored as a vector over k, with each element being a vector of <\\Psi_i \\Psi_j \\Psi_k> values
  Array1D<Array1D<double> > psiIJKProd2_;

  /// \brief i-indices of <\\Psi_i \\Psi_j \\Psi_k \\Psi_l> terms that are not zero, for all l
  /// \note Stored as a vector over l, with each element being a vector of i-indices itself
  Array1D<Array1D<int> > iProd3_;

  /// \brief j-indices of <\\Psi_i \\Psi_j \\Psi_k \\Psi_l> terms that are not zero, for all l
  /// \note Stored as a vector over l, with each element being a vector of j-indices itself
  Array1D<Array1D<int> > jProd3_;

  /// \brief k-indices of <\\Psi_i \\Psi_j \\Psi_k \\Psi_l> terms that are not zero, for all l
  /// \note Stored as a vector over l, with each element being a vector of k-indices itself
  Array1D<Array1D<int> > kProd3_;

  /// \brief <\\Psi_i \\Psi_j \\Psi_k \\Psi_l> terms that are not zero, for all l
  /// \note Stored as a vector over l, with each element being a vector of <\\Psi_i \\Psi_j \\Psi_k \\Psi_l> values
  Array1D<Array1D<double> > psiIJKLProd3_;

  /// \brief Flag for method to compute log: TaylorSeries or Integration
  LogCompMethod logMethod_ ;

  /// \brief CVODE parameter: maximal order
  int CVmaxord_;

  /// \brief CVODE parameter: maximal number of steps
  int CVmaxnumsteps_  ;

  /// \brief CVODE parameter: initial step size
  double CVinitstep_;

  /// \brief CVODE parameter: maximal step size
  double CVmaxstep_;

  /// \brief CVODE parameter: relative tolerance
  double CVrelt_;

  /// \brief CVODE parameter: absolute tolerance
  double CVabst_  ;

  /// \brief Check cvode return for errors
  int Check_CVflag(void *flagvalue, const char *funcname, int opt) const;

  /// \brief Index of this class
  int my_index_;

  /// \brief Number of free parameters to specify the basis
  int narg_;

  /// \brief Parameter alpha for PCs that require a parameter (GLG,SW,JB)
  double alpha_;
  /// \brief Parameter beta for PCs that require two parameters (SW,JB)
  double beta_;

  /// \brief Definition of a map to connect integer indexes with pointers to this class
  typedef std::map<int, PCSet*> OMap_t;
  /// \brief index of next object in map
  static int next_index_;
  /// \brief Map to connect integer indexes with pointers to this class
  static OMap_t *omap_;

};

#endif /* !PCSET_H_SEEN */
