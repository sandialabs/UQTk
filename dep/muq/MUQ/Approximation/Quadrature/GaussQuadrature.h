#ifndef GAUSSQUADRATURE_H_
#define GAUSSQUADRATURE_H_

#include "MUQ/Approximation/Polynomials/OrthogonalPolynomial.h"

#include "MUQ/Approximation/Quadrature/Quadrature.h"

#include <Eigen/Core>
#include <Eigen/Eigenvalues>

namespace muq {

  namespace Approximation {

    /** @ingroup Quadrature
        @brief Class for computing Gauss Quadrature rules from an orthogonal polynomial family.
        @details Uses the Golub-Welsch algorithm to construct a Gauss Quadrature rule of a
        specified order.
    */
    class GaussQuadrature : public Quadrature {

    public:

      GaussQuadrature();

      virtual ~GaussQuadrature() = default;

      GaussQuadrature(std::shared_ptr<OrthogonalPolynomial> polyIn);

      GaussQuadrature(std::shared_ptr<OrthogonalPolynomial> polyIn,
                      int                                   polyOrderIn);

      virtual void Compute(unsigned int quadOrder) override;

      virtual unsigned int Exactness(unsigned int quadOrder) const override{return 2*quadOrder+1;};

    private:

      std::shared_ptr<OrthogonalPolynomial> poly;

      int polyOrder;

    };

  }

}

#endif
