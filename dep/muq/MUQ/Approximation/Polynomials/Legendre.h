#ifndef LEGENDRE_H_
#define LEGENDRE_H_

#include "MUQ/Approximation/Polynomials/OrthogonalPolynomial.h"

namespace muq {
  namespace Approximation {

    /** @ingroup Polynomials
        @class Legendre
        @brief Family of Legendre orthogonal polynomials.
    */
    class Legendre : public OrthogonalPolynomial {
    public:

      /// A Legendre polynomial (\f$1\f$, \f$x\f$, \f$\frac{1}{2}(3x^2-1)\f$, ect. ...)
      /**
	 Legendre polynomials are orthogonal, which helps with some conditioning problems.
       */
      Legendre();

      virtual ~Legendre();

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
