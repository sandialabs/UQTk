#ifndef HERMITEFUNCTION_H_
#define HERMITEFUNCTION_H_


#include "MUQ/Approximation/Polynomials/IndexedScalarBasis.h"
#include "MUQ/Approximation/Polynomials/PhysicistHermite.h"

namespace muq {
  namespace Approximation {

    /** @ingroup Polynomials
         @class HermiteFunction
         @brief A 1D hermite function based on Physicist Hermite Polynomials
    */
    class HermiteFunction : public IndexedScalarBasis {
    public:


      /// Create a polynomial
      HermiteFunction() : IndexedScalarBasis(),
                          polyBase(std::make_shared<PhysicistHermite>()){};

      virtual ~HermiteFunction() = default;

      /// Evaluate the hermite function
      virtual double BasisEvaluate(int const order, double const x) const override;

      /** Evaluate the \f$n^{th}\f$ derivative of a \f$p^{th}\f$ order hermite function. */
      virtual double DerivativeEvaluate(int const polyOrder, int const derivOrder, double const x) const override;

    private:

      static unsigned nChoosek(unsigned n, unsigned k);

      std::shared_ptr<PhysicistHermite> polyBase;

    };
  } // namespace Approximation
} // namespace muq


#endif
