#ifndef MODPIECECOSTFUNCTION_H_
#define MODPIECECOSTFUNCTION_H_

#include "MUQ/Optimization/CostFunction.h"

namespace muq {
  namespace Optimization {
    class ModPieceCostFunction : public CostFunction {
    public:

      ModPieceCostFunction(std::shared_ptr<muq::Modeling::ModPiece> cost);

      virtual ~ModPieceCostFunction() = default;
    private:

      /// The value of the cost function
      /**
      @param[in] args The inputs \f$x\f$, \f$\theta_{1:n}\f$
      \return The value of the cost function
       */
      virtual double CostImpl(muq::Modeling::ref_vector<Eigen::VectorXd> const& input) override;

      /// Compute the gradient of the cost function
      /**
      @param[in] inputDimWrt Which input are we taking the derivative with respect to?
      @param[in] args The inputs \f$x\f$, \f$\theta_{1:n}\f$ and the sensitivity vector
       */
      virtual void GradientImpl(unsigned int const inputDimWrt, muq::Modeling::ref_vector<Eigen::VectorXd> const& input, Eigen::VectorXd const& sensitivity) override;

      // The muq::Modeling::ModPiece that holds the cost
      std::shared_ptr<muq::Modeling::ModPiece> cost;
    };
  } // namespace Optimization
} // namespace muq

#endif
