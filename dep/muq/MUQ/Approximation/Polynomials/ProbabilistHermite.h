#ifndef PROBABILISTHERMITE_H_
#define PROBABILISTHERMITE_H_

#include "MUQ/Approximation/Polynomials/OrthogonalPolynomial.h"

namespace muq {
  namespace Approximation {
    class ProbabilistHermite : public OrthogonalPolynomial{
    public:

      /**
      @ingroup Polynomials
      @class ProbabilistHermite
      @brief A Probabilist Hermite polynomial (\f$1\f$, \f$2x\f$, \f$4x^2-2.0\f$, ect. ...)
      @details Probabilist Hermite polynomials are orthogonal with respect to \f$\exp\left[-\frac{1}{2}x^2\right]\f$, which helps with some conditioning problems.   Here we implement the probabilists' Hermite polynomials.
               note: Hermite may be unstable for high orders.
       */
      ProbabilistHermite() = default;

      virtual ~ProbabilistHermite() = default;

      virtual double DerivativeEvaluate(int const polyOrder, int const derivOrder, double const x) const override;

      virtual double Normalization(unsigned int polyOrder) const override;

    private:

      virtual double ak(unsigned int k) const override;
      virtual double bk(unsigned int k) const override;
      virtual double ck(unsigned int k) const override;
      virtual double phi0(double x) const override;
      virtual double phi1(double x) const override;
    };
  } // namespace Approximation
} // namespace muq

#endif
