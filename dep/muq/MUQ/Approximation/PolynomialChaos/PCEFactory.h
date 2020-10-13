#ifndef PCEFACTORY_H
#define PCEFACTORY_H

#include "MUQ/Approximation/PolynomialChaos/PolynomialChaosExpansion.h"
#include "MUQ/Approximation/Quadrature/Quadrature.h"
#include "MUQ/Approximation/Quadrature/FullTensorQuadrature.h"

#include "MUQ/Utilities/MultiIndices/MultiIndex.h"

#include "MUQ/Modeling/ModPiece.h"

#include <boost/property_tree/ptree.hpp>

namespace muq {
namespace Approximation {

  /** @class PCEFactory
      @ingroup PolynomialChaos
      @brief Factory class for constructing a pseudo-spectral polynomial chaos
             approximation using a fixed quadrature rule.
  */
  class PCEFactory {
  public:

    PCEFactory(std::vector<std::shared_ptr<Quadrature>>         const& quadTypesIn,
               std::vector<std::shared_ptr<IndexedScalarBasis>> const& polyTypesIn);

    /** Constructs a factory using a prespecified tensor product quadrature rule.
        The terms in the polynomial expansion are chosen from the quadrature rule.
        Polynomial terms that can be integrated exactly are used.
        @param[in] quadTypesIn A vector of one dimensional quadrature rules for each input dimension.
        @param[in] quadOrders A multiindex holding the quadrature order in each direction.
        @param[in] polyTypesIn A vector of one dimensional polynomial families
                   used to define the multivariate polynomials used in a PCE expansion.
    */
    PCEFactory(std::vector<std::shared_ptr<Quadrature>>         const& quadTypesIn,
               std::shared_ptr<muq::Utilities::MultiIndex>      const& quadOrders,
               std::vector<std::shared_ptr<IndexedScalarBasis>> const& polyTypesIn);

    /** Constructs a factory using a tensor product quadrature rule and independently
        specified polynomial terms.  This allows users to employ the "classic"
        pseudo-spectral approach for computing polynomial chaos expansion.  Coefficients
        for each terms defined in polyTypesIn will be computed using the quadrature
        rule defined by the tensor product quadrature rule defined by quadTypesIn
        and quadOrders.  Note that you should carefully choose polyMultisIn and
        quadTypesIn because significant approximation errors can arise if terms
        in the polynomial approximation cannot be integrated exactly by the quadrature
        rule.

        @param[in] quadTypesIn A vector of one dimensional quadrature rules for each input dimension.
        @param[in] quadOrders A multiindex holding the quadrature order in each direction.
        @param[in] polyTypesIn A vector of one dimensional polynomial families
                   used to define the multivariate polynomials used in a PCE expansion.
        @param[in] polyMultisIn A multindex set defining the polynomial expansion
                   this factory should generate.
    */
    PCEFactory(std::vector<std::shared_ptr<Quadrature>>         const& quadTypesIn,
               std::shared_ptr<muq::Utilities::MultiIndex>      const& quadOrders,
               std::vector<std::shared_ptr<IndexedScalarBasis>> const& polyTypesIn,
               std::shared_ptr<muq::Utilities::MultiIndexSet>   const& polyMultisIn);


    /**
      Computes a polynomial chaos approximation to the given model.
      @param[in] model A ModPiece implementing the model we wish to approximate.
      @return A polynomial chaos expansion constructed from the polynomials and
              qudrature rules given to the PCEFactory constructor.
    */
    std::shared_ptr<PolynomialChaosExpansion> Compute(std::shared_ptr<muq::Modeling::ModPiece> const& model);


    /**
    Computes a polynomial chaos approximation using evaluations of the model at
    the quadrature points returned by PCEFactory::QuadPts.  In most cases, it is
    preferrable to call the Compute function with ModPiece argument, which will
    help avoid potential issues matching the model evaluation points with the
    quadrature points.   However, in some cases (e.g., when parallel model evaluations
    used) it is advantageous to manually compute the model evaluations and compute
    the expansion with this function.   This function is typically used in a workflow
    like the following code snippet:
@code{cpp}
// Set up quadrature types, orders, etc...
PCEFactory factory(quadTypes, quadOrders, polyTypes);

std::vector<Eigen::VectorXd> const& quadPts = factory.QuadPts();

std::vector<Eigen::VectorXd> const& quadEvals(quadPts.size());
for(int i=0; i<quadPts.size(); ++i)
  quadEvals.at(i) = model->Evaluate(quadPts.at(i))(0);

std::shared_ptr<PolynomialChaosExpansion> pce = factory.Compute(quadEvals);

@endcode
    */
    std::shared_ptr<PolynomialChaosExpansion> Compute(std::vector<Eigen::VectorXd> const& quadEvals);

    std::shared_ptr<PolynomialChaosExpansion> Compute(std::vector<std::reference_wrapper<const Eigen::VectorXd>> const& quadEvals);

    std::shared_ptr<PolynomialChaosExpansion> Compute(std::vector<Eigen::VectorXd>                const& quadEvals,
                                                      std::shared_ptr<muq::Utilities::MultiIndex> const& quadOrders);

    std::shared_ptr<PolynomialChaosExpansion> Compute(std::vector<std::reference_wrapper<const Eigen::VectorXd>> const& quadEvals,
                                                      std::shared_ptr<muq::Utilities::MultiIndex> const& quadOrders);
    /**
      Returns the quadrature points in the tensor product quadrature rule.
    */
    std::vector<Eigen::VectorXd> const& QuadPts() const{return quadPts;};
    std::vector<Eigen::VectorXd> const& QuadPts(std::shared_ptr<muq::Utilities::MultiIndex> const& quadOrders);


  protected:


    /** Sets up the tensor product quadrature and polynomial expansions based on
      the specified quadrature order.
    */
    void Setup(std::shared_ptr<muq::Utilities::MultiIndex> const& quadOrders);

    std::shared_ptr<muq::Utilities::MultiIndex> quadOrdersCache;


    std::vector<std::shared_ptr<Quadrature>> quadTypes;
    std::vector<std::shared_ptr<IndexedScalarBasis>> polyTypes;
    FullTensorQuadrature tensQuad;

    std::shared_ptr<muq::Utilities::MultiIndexSet>   polyMultis;

    std::vector<Eigen::VectorXd> quadPts;
    Eigen::VectorXd quadWts;

  };

}
}




#endif // PCEFACTORY_H
