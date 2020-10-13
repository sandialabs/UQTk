#ifndef ADAPTIVESMOLYAKPCE_H
#define ADAPTIVESMOLYAKPCE_H

#include "MUQ/Approximation/PolynomialChaos/SmolyakEstimator.h"
#include "MUQ/Approximation/PolynomialChaos/PolynomialChaosExpansion.h"
#include "MUQ/Approximation/PolynomialChaos/PCEFactory.h"

#include "MUQ/Approximation/Quadrature/Quadrature.h"

namespace muq {
namespace Approximation {

  class AdaptiveSmolyakPCE : public SmolyakEstimator<std::shared_ptr<PolynomialChaosExpansion>> {

  public:
    AdaptiveSmolyakPCE(std::shared_ptr<muq::Modeling::ModPiece>         const& modelIn,
                       std::vector<std::shared_ptr<Quadrature>>         const& quad1dIn,
                       std::vector<std::shared_ptr<IndexedScalarBasis>> const& polys1dIn);

    virtual ~AdaptiveSmolyakPCE() = default;

  protected:
    virtual std::vector<Eigen::VectorXd> OneTermPoints(std::shared_ptr<muq::Utilities::MultiIndex> const& multi) override;

    virtual std::shared_ptr<PolynomialChaosExpansion> ComputeOneTerm(std::shared_ptr<muq::Utilities::MultiIndex>                const& multi,
                                                                     std::vector<std::reference_wrapper<const Eigen::VectorXd>> const& modEvals) override;

    virtual std::shared_ptr<PolynomialChaosExpansion> AddEstimates(double w1, std::shared_ptr<PolynomialChaosExpansion> const& part1,
                                                                   double w2, std::shared_ptr<PolynomialChaosExpansion> const& part2) const override;

    virtual double ComputeMagnitude(std::shared_ptr<PolynomialChaosExpansion> const& estimate) const override;

    PCEFactory tensFactory;

  }; // class AdaptiveSmolyakPCE

} // namespace Approximation
} // namespace muq


#endif
