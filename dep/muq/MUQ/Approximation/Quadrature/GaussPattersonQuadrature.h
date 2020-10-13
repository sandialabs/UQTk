#ifndef GAUSSPATTERSONQUADRATURE_H
#define GAUSSPATTERSONQUADRATURE_H

#include "MUQ/Approximation/Quadrature/Quadrature.h"

namespace muq {
namespace Approximation {

  /** @class GaussPattersonQuadrature
      @ingroup Quadrature
      @brief 1d Gauss Patterson nested quadrature rule
      @details
  */
  class GaussPattersonQuadrature : public Quadrature {
  public:

    GaussPattersonQuadrature();

    virtual ~GaussPattersonQuadrature() = default;

    virtual void Compute(unsigned int index) override;

    virtual unsigned int Exactness(unsigned int quadOrder) const override;

  };

} // namespace muq
} // namespace Approximation


#endif
