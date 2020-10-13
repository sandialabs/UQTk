#ifndef FULLTENSORQUADRATURE_H_
#define FULLTENSORQUADRATURE_H_

#include <Eigen/Core>
#include <vector>

#include "MUQ/Approximation/Quadrature/Quadrature.h"


namespace muq {
  namespace Approximation {

/**
   @class FullTensorQuadrature
   @ingroup Quadrature
   @brief Multivariate quadrature rule defined by the tensor product of 1d rules.
   @details One of the most straightforward multivariate quadrature rules, a full tensor expansion.
 * Creates a the complete tensor product of the points of the given order
 * 1D quadrature rules for each dimension. The 1D quadrature order can
 * be isotropic or vary per dimension. Note that an isotropic order
 * does not mean isotropic number of points in each dimension, if different
 * 1D quadrature rules are used.
 *
 * Probably not the most useful class for actual analysis, because you have to know exactly
 * what order you want, but is an important building block for sparse
 * quadrature routines.
 */
class FullTensorQuadrature : public Quadrature {
public:

  /**
  Sets up a tensor product of quadrature rules in dim dimensions.  Does not
  compute the tensor product because no order is specified.
  */
  FullTensorQuadrature(unsigned int                       dim,
                       std::shared_ptr<Quadrature> const& rules);

  /** Sets up and solves a tensor product of 1d quadrature rules.
  */
  FullTensorQuadrature(unsigned int                       dim,
                       std::shared_ptr<Quadrature> const& rules,
                       unsigned int                       order);

  /** Sets up a tensor product of quadrature rules with potentially different rules
      in each dimension.  If the orders row vector is specified, then the tensor
      product rule is also computed, otherwise an additional call to "Compute" is
      required.
  */
  FullTensorQuadrature(std::vector<std::shared_ptr<Quadrature>> const& rules,
                       Eigen::RowVectorXi                              orders = Eigen::RowVectorXi());


  virtual ~FullTensorQuadrature() = default;

  virtual void Compute(unsigned int order) override;
  virtual void Compute(Eigen::RowVectorXi const& orders) override;

  virtual unsigned int Exactness(unsigned int quadOrder) const override;

private:

  std::vector<std::shared_ptr<Quadrature>> rules;
};
}
}

#endif /* FULLTENSORQUADRATURE_H_ */
