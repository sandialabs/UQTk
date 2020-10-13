#ifndef CLENSHAWCURTISQUADRATURE_H
#define CLENSHAWCURTISQUADRATURE_H

#include "MUQ/Approximation/Quadrature/Quadrature.h"

namespace muq {
namespace Approximation {

  /** @class ClenshawCurtisQuadrature
      @ingroup Quadrature
      @brief 1d Clenshaw Curtis rule
      @details

      As discussed for the "CCN_RULE" at http://people.sc.fsu.edu/~jburkardt/cpp_src/ccn_rule/ccn_rule.html,
    The Clenshaw Curtis quadrature rule is only nested for certain orders, specifically,
    orders [1,3,5,9,17,33,...].  This class allows users to define a nested
    quadrature rule, where the index passed to `Compute` is computed to one of
    these nested orders, or a non-nested quadrature rule, where the order is
    just equal to the index.  If a nested rule is desired (the defaul behavior),
    "true" should be passed to the constructor.

<b>Construction of nested rule.</b>
@code{.cpp}
ClenshawCurtisQuadrature quad; // or equivalently quad(true)

// Compute a nested rule with index 3, which will have 9 points
unsigned int index = 3;
quad.Compute(index);

Eigen::MatrixXd pts = quad.Points();
Eigen::VectorXd wts = quad.Weights();
@endcode

<b>Construction of non-nested rule.</b>
@code{.cpp}
ClenshawCurtisQuadrature quad(false);

// Compute a nested rule with index 3, which will have 3+1 points
unsigned int index = 3;
quad.Compute(index);

Eigen::MatrixXd pts = quad.Points();
Eigen::VectorXd wts = quad.Weights();
@endcode
  */
  class ClenshawCurtisQuadrature : public Quadrature {
  public:

    ClenshawCurtisQuadrature(bool nestedIn=true);

    virtual ~ClenshawCurtisQuadrature() = default;

    virtual void Compute(unsigned int index) override;

    virtual unsigned int Exactness(unsigned int quadOrder) const override;
  private:


    unsigned int IndexToNumPoints(unsigned int index) const;

    const bool nested;
    const double pi = 3.14159265358979323846;
  };

} // namespace muq
} // namespace Approximation


#endif
