#ifndef MONOTONEEXPANSION_H
#define MONOTONEEXPANSION_H

#include "MUQ/Approximation/Polynomials/BasisExpansion.h"
#include <boost/property_tree/ptree.hpp>

namespace muq{
  namespace Approximation{

    /** @class MonotoneExpansion
        @ingroup Polynomials
        @brief Defines a parameterized family of multivariate lower triangular monotone functions.
        @details This class defines functions that take the form
        \f[ f_i(x) = g_i(x_1, \ldots, x_{i-1}) + \int_{-\infty}^{x_i} \left[ h_i(x_1,\ldots, x_{i-1}, y)\right]^2 + \epsilon dy, \f]
        where the functions \f$g_i\f$ and \f$h_i\f$ are both represented through the
        muq::Approximation::BasisExpansion class, and \f$\epsilon>0\f$ is a small nugget
        ensuring that \f$f_i\f$ is strictly increasing with \f$x_i\f$.  By squaring \f$h_i\f$,
        we guarantee that the integrand is positive, which ensures that the integral itself increases
        with \f$x_i\f$, i.e., \f$\partial f_i / \partial x_i > 0 \f$.  Also notice how
        the \f$i^{th}\f$ output of the function \f$f\f$ only depends on the first \f$i\$ inputs,
        so the function \f$f\f$ is lower triangular. Combined with the fact that \f$\partial f_i / \partial x_i > 0 \f$,
        this means that \f$f\f$ will be invertible.
    */
    class MonotoneExpansion : public muq::Modeling::ModPiece{

    public:
      //MonotoneExpansion(boost::property_tree::ptree & params);

      MonotoneExpansion(std::shared_ptr<BasisExpansion> monotonePartsIn,
                        bool                            coeffInput=false);

      MonotoneExpansion(std::vector<std::shared_ptr<BasisExpansion>> const& generalPartsIn,
                        std::vector<std::shared_ptr<BasisExpansion>> const& monotonePartsIn,
                        bool                                                coeffInput=false);

      /** Returns the number of coefficients across all expansions in all dimension, i.e.,
          the total number of degrees of freedom describing this expansion.
      */
      unsigned NumTerms() const;

      /** Given a monotone expansion with components T(x) = [T_1(x_1), T_2(x_1,x_2), ..., T_N(x_1,...,x_N)],
          this function returns a new monotone expansion with just the first components:
          T_{new}(x) = [T_1(x_1), T_2(x_1,x_2), ..., T_M(x_1,...,x_M)].  The number
          of terms to keep, M, is the sole input to this function.
      */
      std::shared_ptr<MonotoneExpansion> Head(int numRows) const;

      /** Evaluate the inverse of this map. */
      Eigen::VectorXd EvaluateInverse(Eigen::VectorXd const& refPt) const;
      Eigen::VectorXd EvaluateInverse(Eigen::VectorXd const& refPt,
                                      Eigen::VectorXd const& tgtPt0) const;

      Eigen::VectorXd EvaluateForward(Eigen::VectorXd const& x) const;

      /** Get the current expansion coefficients.  The coefficients are ordered
          with coefficients from the general parts preceding coefficients for
          the monotone parts.  Within the general and monotone parts, the coefficients
          are ordered according to the output dimension.
      */
      virtual Eigen::VectorXd GetCoeffs() const;

      virtual void SetCoeffs(Eigen::VectorXd const& allCoeffs);

      /** Returns the log determinant of the Jacobian matrix (wrt x) at a particular point. */
      virtual double LogDeterminant(Eigen::VectorXd const& evalPt);
      virtual double LogDeterminant(Eigen::VectorXd const& evalPt,
                                    Eigen::VectorXd const& coeffs);

      /** Returns the gradient of the Jacobian log determinant with respect to the coefficients. */
      virtual Eigen::VectorXd GradLogDeterminant(Eigen::VectorXd const& evalPt);
      virtual Eigen::VectorXd GradLogDeterminant(Eigen::VectorXd const& evalPt,
                                                 Eigen::VectorXd const& coeffs);

    protected:

      virtual void EvaluateImpl(muq::Modeling::ref_vector<Eigen::VectorXd> const& inputs) override;

      virtual void JacobianImpl(unsigned int const                                wrtOut,
                                unsigned int const                                wrtIn,
                                muq::Modeling::ref_vector<Eigen::VectorXd> const& inputs) override;

      Eigen::MatrixXd JacobianWrtX(Eigen::VectorXd const& x) const;

      std::vector<std::shared_ptr<BasisExpansion>> generalParts;
      std::vector<std::shared_ptr<BasisExpansion>> monotoneParts;

      Eigen::VectorXd quadWeights;
      Eigen::VectorXd quadPts; // quadrature points on [0,1]

    private:
      static Eigen::VectorXi GetInputSizes(std::vector<std::shared_ptr<BasisExpansion>> const& generalPartsIn,
                                           std::vector<std::shared_ptr<BasisExpansion>> const& monotonePartsIn,
                                           bool                                                coeffInput);

    }; // class MonotoneExpansion

  } // namespace Approximation
} // namespace muq


#endif
