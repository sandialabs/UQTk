#ifndef EXPONENTIALGROWTHQUADRATURE_H
#define EXPONENTIALGROWTHQUADRATURE_H

#include "MUQ/Approximation/Quadrature/Quadrature.h"

namespace muq {
namespace Approximation {

  /** @class ExponentialGrowthQuadrature
      @ingroup Quadrature
      @brief 1d Quadrature rule with exponential growth
      @details In many cases, such as pseudo-spectral constructions of polynomial
      chaos expansions, it can be advantageous for the quadrature order to grow
      faster than normal.  This rule facilitates that by wrapping around another
      one dimensional quadrature rule but transforming the index.  If the original
      quadrature rule with index \f$j\f$ is denoted by \f$Q_0(j)\f$, then this
      rule will return \f$q_1(k) = Q(2^k)\f$ for a specified index $k$.
  */
  class ExponentialGrowthQuadrature : public Quadrature {
  public:

    ExponentialGrowthQuadrature(std::shared_ptr<Quadrature> const& quadIn);

    virtual ~ExponentialGrowthQuadrature() = default;

    virtual void Compute(unsigned int index) override;

    virtual unsigned int Exactness(unsigned int quadOrder) const override;

    virtual Eigen::MatrixXd const& Points() const override;

    virtual Eigen::VectorXd const& Weights() const override;

  protected:
    std::shared_ptr<Quadrature> otherQuad;

  };

} // namespace muq
} // namespace Approximation


#endif
