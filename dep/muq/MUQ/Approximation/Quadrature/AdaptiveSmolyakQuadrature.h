#ifndef ADAPTIVESMOLYAKPCE_H
#define ADAPTIVESMOLYAKPCE_H

#include "MUQ/Approximation/PolynomialChaos/SmolyakEstimator.h"
#include "MUQ/Approximation/PolynomialChaos/PolynomialChaosExpansion.h"

#include "MUQ/Approximation/Quadrature/FullTensorQuadrature.h"

namespace muq {
namespace Approximation {

  class AdaptiveSmolyakQuadrature : public SmolyakEstimator<Eigen::VectorXd> {

  public:
    AdaptiveSmolyakQuadrature(std::shared_ptr<muq::Modeling::ModPiece> const& modelIn,
                              std::vector<std::shared_ptr<Quadrature>> const& quad1d);


    virtual ~AdaptiveSmolyakQuadrature() = default;

  protected:
    virtual std::vector<Eigen::VectorXd> OneTermPoints(std::shared_ptr<muq::Utilities::MultiIndex> const& multi) override;

    virtual Eigen::VectorXd ComputeOneTerm(std::shared_ptr<muq::Utilities::MultiIndex>                const& multi,
                                           std::vector<std::reference_wrapper<const Eigen::VectorXd>> const& modEvals) override;

    virtual Eigen::VectorXd AddEstimates(double w1, Eigen::VectorXd const& part1,
                                         double w2, Eigen::VectorXd const& part2) const override;

    virtual double ComputeMagnitude(Eigen::VectorXd const& estimate) const override;

    std::shared_ptr<muq::Utilities::MultiIndex> cachedMulti;
    FullTensorQuadrature tensQuad;


  }; // class AdaptiveSmolyakQuadrature

} // namespace Approximation
} // namespace muq


#endif
