#ifndef REGRESSION_H_
#define REGRESSION_H_

#include <boost/property_tree/ptree.hpp>

#include <Eigen/QR>

#include "MUQ/Utilities/MultiIndices/MultiIndexFactory.h"

#include "MUQ/Modeling/LinearAlgebra/AnyAlgebra.h"
#include "MUQ/Modeling/WorkPiece.h"

#include "MUQ/Optimization/CostFunction.h"

#include "MUQ/Approximation/Polynomials/IndexedScalarBasis.h"

namespace muq {
  namespace Approximation {
    class Regression : public muq::Modeling::WorkPiece, public std::enable_shared_from_this<Regression> {
    public:

      /**
	 <ol>
	 <li> order The order of the polynomial regression (<EM>Order</EM>)
	 <li> The type of polynomial basis to use (defaults to Legendre) (<EM>PolynomialBasis</EM>)
	 </ol>
	 @param[in] pt Options for the regression
      */
      Regression(boost::property_tree::ptree const& pt);

      /// Compute the coeffiecents of the polynomial given data
      /**
	 @param[in] xs The input points
	 @param[in] ys The output points
	 @param[in] center The center of the inputs (used to recenter the inputs)
      */
      void Fit(std::vector<Eigen::VectorXd> xs, std::vector<Eigen::VectorXd> const& ys, Eigen::VectorXd const& center);

      /// Compute the coeffiecents of the polynomial given data
      /**
	 @param[in] xs The input points
	 @param[in] ys The output points
       */
      void Fit(std::vector<Eigen::VectorXd> const& xs, std::vector<Eigen::VectorXd> const& ys);

      int NumInterpolationPoints() const;

       /**
	      @param[in] xs The input points in the ball (not necessarily normalized to the unit ball)
	       @param[in] center The center of the ball
         @param[in] kn Maximize the Lagrange polynomials in the smallest ball containing this number of points.  If \f$k_n<0\f$ (default) use the smallest ball that contains all of the points.  It may not be zero.
	        \return A tuple: poisedness constant, radius, index of input point associated with the poisedness consant, the location that maximizes the Lagrange polynomial.
       */
      std::pair<Eigen::VectorXd, double> PoisednessConstant(std::vector<Eigen::VectorXd> xs, Eigen::VectorXd const& center, int kn = -1) const;

      /// The order of the regression
      const unsigned int order;

    private:

      virtual void EvaluateImpl(muq::Modeling::ref_vector<boost::any> const& inputs) override;

      void ComputeBasisDerivatives(Eigen::VectorXd const& point, std::vector<Eigen::VectorXd>& gradient) const;

      /// Compute the coefficients for the basis functions
      /**
	 Given points, data, and a center compute the coefficients on the basis (inner product with the basis evalautes the local polynomial).  The data can be a have multiple outputs (e.g., fitting more than one polynomial at once), which leads to more than one set of basis coefficents.
	 @param[in] xs The points
	 @param[in] ys The output at each point
	 \return The coefficients for the basis
      */
      Eigen::MatrixXd ComputeCoefficients(std::vector<Eigen::VectorXd> const& xs, std::vector<Eigen::VectorXd> const& ys) const;

      /// Compute the right hand side given data to compute the polynomial coefficients
      /**
	 @param[in] vand The Vandermonde matrix
      	 @param[in] ys_data The output at each point
      */
      Eigen::MatrixXd ComputeCoefficientsRHS(Eigen::MatrixXd const& vand, std::vector<Eigen::VectorXd> const& ys_data) const;

      /// Create the Vandermonde matrix
      /**
	 @param[in] xs The points
	 \return The Vandermonde matrix
      */
      Eigen::MatrixXd VandermondeMatrix(std::vector<Eigen::VectorXd> const& xs) const;

      /// Center the input points
      /**
      Center the input points around currentCenter
      @param[in] xs Input points
      */
      Eigen::ArrayXd CenterPoints(std::vector<Eigen::VectorXd>& xs);

      /// Center the input points
      /**
      @param[in] xs Input points
      @param[in] center The center point
      */
      Eigen::ArrayXd CenterPoints(std::vector<Eigen::VectorXd>& xs, Eigen::VectorXd const& center) const;

      /// Center the input points
      /**
      Center the input points such that the first \f$k_n\f$ of them are in a unit ball (the rest may be inside or outside, depending on the their magnitude.)
      @param[in] xs Input points
      @param[in] center The center point
      @param[in] kn The first \f$kn\f$ points are in the unit ball
      */
      Eigen::ArrayXd CenterPoints(std::vector<Eigen::VectorXd>& xs, Eigen::VectorXd const& center, unsigned int const kn) const;

      class PoisednessCost : public muq::Optimization::CostFunction {
      public:

	       PoisednessCost(std::shared_ptr<Regression const> parent, std::vector<Eigen::RowVectorXd> const& lagrangeCoeff, unsigned int const inDim);

	        virtual ~PoisednessCost() = default;

      private:

	       virtual double CostImpl(muq::Modeling::ref_vector<Eigen::VectorXd> const& input) override;

	        virtual void GradientImpl(unsigned int const inputDimWrt, muq::Modeling::ref_vector<Eigen::VectorXd> const& input, Eigen::VectorXd const& sensitivity) override;

	        std::shared_ptr<Regression const> parent;

	        const std::vector<Eigen::RowVectorXd>& lagrangeCoeff;
      };

      class PoisednessConstraint : public muq::Modeling::ModPiece {
      public:

	       PoisednessConstraint(unsigned int const inDim, double const alpha);

	        virtual ~PoisednessConstraint() = default;

      private:

	virtual void EvaluateImpl(muq::Modeling::ref_vector<Eigen::VectorXd> const& input) override;

	virtual void JacobianImpl(unsigned int const outputDimWrt,
                                  unsigned int const inputDimWrt,
                                  muq::Modeling::ref_vector<Eigen::VectorXd> const& input) override;

          const double alpha;
      };

      /// The input dimension
      const unsigned int inputDim;

      /// The radius of the poisedness constraint
      /**
      Defaults to 1.  Max is 1.
      */
      const double alpha;

      /// The multi-index to so we know the order of each term
      std::shared_ptr<muq::Utilities::MultiIndexSet> multi;

      /// The polynomial basis (in one variable) used to compute the Vandermonde matrix
      std::shared_ptr<IndexedScalarBasis> poly;

      /// Current center of the inputs
      Eigen::VectorXd currentCenter;

      /// Current radius of inputs
      Eigen::ArrayXd currentRadius;

      /// Coeffients for the polynomial basis
      Eigen::MatrixXd coeff;

      /// Parameters for the poisedness cosntant optimization
      boost::property_tree::ptree optPt;
    };
  } // namespace Approximation
} // namespace muq

#endif
