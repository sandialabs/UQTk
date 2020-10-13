#ifndef JACOBI_H_
#define JACOBI_H_

#include "MUQ/Approximation/Polynomials/OrthogonalPolynomial.h"

namespace muq {
    namespace Approximation {

        /** @ingroup Polynomials
            @class Jacobi
            @brief Family of Jacobi orthogonal polynomials.
        */
        class Jacobi : public OrthogonalPolynomial {
        public:

            /**
               Jacobi polynomials are orthogonal over \f$[-1,1]\f$ with respect to the weighting function \f$(1-x)^a(1+x)^b\f$.  Note that Jacobi polynomials are a generalization of Legendre and Chebyshev polynomials.
            */
            Jacobi(const double aIn=1.0, const double bIn=1.0) : a(aIn), b(bIn){};

            virtual ~Jacobi() = default;

            virtual double DerivativeEvaluate(int const polyOrder, int const derivOrder, double const x) const override;

            virtual double Normalization(unsigned int polyOrder) const override;

        private:

            const double a;
            const double b;

            virtual double ak(unsigned int k) const override;
            virtual double bk(unsigned int k) const override;
            virtual double ck(unsigned int k) const override;
            virtual double phi0(double x) const override;
            virtual double phi1(double x) const override;
        };
    } // namespace Approximation
} // namespace muq

#endif
