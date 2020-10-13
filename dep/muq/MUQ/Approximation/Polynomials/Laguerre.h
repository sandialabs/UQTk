#ifndef LAGUERRE_H_
#define LAGUERRE_H_

#include "MUQ/Approximation/Polynomials/OrthogonalPolynomial.h"

namespace muq {
    namespace Approximation {

      /** @ingroup Polynomials
          @class Laguerre
          @brief Family of Laguerre orthogonal polynomials.
      */
        class Laguerre : public OrthogonalPolynomial {
        public:

            /**
               Laguerre polynomials are orthogonal over \f$[0,\infty)\f$ with respect to \f$\exp(-x)x^a\f$.
            */
            Laguerre(const double aIn=0.0) : a(aIn){};

            virtual ~Laguerre() = default;

            virtual double DerivativeEvaluate(int const polyOrder, int const derivOrder, double const x) const override;

            virtual double Normalization(unsigned int polyOrder) const override;

        private:

            const double a;

            virtual double ak(unsigned int k) const override;
            virtual double bk(unsigned int k) const override;
            virtual double ck(unsigned int k) const override;
            virtual double phi0(double x) const override;
            virtual double phi1(double x) const override;
        };
    } // namespace Approximation
} // namespace muq

#endif
