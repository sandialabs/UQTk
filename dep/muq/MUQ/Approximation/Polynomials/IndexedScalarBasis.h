#ifndef INDEXEDSCALARBASIS_H_
#define INDEXEDSCALARBASIS_H_

#include <functional>
#include <string>

#include "MUQ/Modeling/WorkPiece.h"
#include "MUQ/Utilities/RegisterClassName.h"

namespace muq {
  namespace Approximation {

    /**
       @ingroup Polynomials
       @class IndexScalarBasis
       @brief An abstract basis for scalar basis functions
       @details Many basis functions can be naturally indexed by an integer.  For example,
       polynomial families are naturally indexed by their order.  This class
       defines an abstract base class for polynomials and other functions with
       scalar arguments that can be naturally index by an integer.
     */
    class IndexedScalarBasis : public muq::Modeling::WorkPiece {
    public:


      typedef std::function<std::shared_ptr<IndexedScalarBasis>()> ScalarBasisConstructorType;
      typedef std::map<std::string, ScalarBasisConstructorType> ScalarBasisMapType;
      static std::shared_ptr<ScalarBasisMapType> GetScalarBasisMap();

      // Factory method for consttructing polynomial from family
      static std::shared_ptr<IndexedScalarBasis> Construct(std::string const& polyName);


      /// Create a basis function
      IndexedScalarBasis();

      virtual ~IndexedScalarBasis() = default;

      /// Evaluate the specific basis type (must be implemented by the child)
      /** @param[in] order The integer index of the polynomial (e.g., polynomial order)
	        @param[in] x The point where we are evaluating the polynomial
	        @return The polynomial value
       */
      virtual double BasisEvaluate(int    const order,
                                   double const x) const = 0;


      /// Evaluates all basis functions with order <= maxOrder
      virtual Eigen::VectorXd EvaluateAllTerms(int    const maxOrder,
                                               double const x) const
      {
        Eigen::VectorXd output(maxOrder+1);
        for(int i=0; i<=maxOrder; ++i)
          output(i) = BasisEvaluate(i,x);
        return output;
      };

      /** Evaluate the \f$n^{th}\f$ derivative of a \f$p^{th}\f$ basis.
          @param[in] polyOrder The integer index \f$p\f$ of the basis.
          @param[in] derivOrder The order \f$n\f$ of the derivative.
          @param[in] x The location to evaluate the derivative.
          @return The the \f$n^{th}\f$ derivative of a \f$p^{th}\f$ basis
      */
      virtual double DerivativeEvaluate(int    const polyOrder,
                                        int    const derivOrder,
                                        double const x) const = 0;


    private:

      /// Evaluate the polynomial at a given point and order
      virtual void EvaluateImpl(muq::Modeling::ref_vector<boost::any> const& inputs);


    };
  } // namespace Approximation
} // namespace muq


#define REGISTER_SCALARBASIS_FAMILY(NAME) static auto regScalarBasis ##NAME		\
    = muq::Approximation::IndexedScalarBasis::GetScalarBasisMap()->insert(std::make_pair(#NAME, muq::Utilities::shared_factory<NAME>()));



#endif
