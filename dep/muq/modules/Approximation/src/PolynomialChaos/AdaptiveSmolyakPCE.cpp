#include "MUQ/Approximation/PolynomialChaos/AdaptiveSmolyakPCE.h"


using namespace muq::Approximation;
using namespace muq::Utilities;





AdaptiveSmolyakPCE::AdaptiveSmolyakPCE(std::shared_ptr<muq::Modeling::ModPiece>         const& modelIn,
                                       std::vector<std::shared_ptr<Quadrature>>         const& quad1dIn,
                                       std::vector<std::shared_ptr<IndexedScalarBasis>> const& polys1dIn)
                        : SmolyakEstimator<std::shared_ptr<PolynomialChaosExpansion>>(modelIn),
                          tensFactory(quad1dIn, polys1dIn)
{}



std::vector<Eigen::VectorXd> AdaptiveSmolyakPCE::OneTermPoints(std::shared_ptr<MultiIndex> const& multi)
{
  return tensFactory.QuadPts(multi);
}

std::shared_ptr<PolynomialChaosExpansion> AdaptiveSmolyakPCE::ComputeOneTerm(std::shared_ptr<MultiIndex>         const& multi,
                                                                             std::vector<std::reference_wrapper<const Eigen::VectorXd>> const& modEvals)
{
  return tensFactory.Compute(modEvals,multi);
}

std::shared_ptr<PolynomialChaosExpansion> AdaptiveSmolyakPCE::AddEstimates(double w1,
                                                                           std::shared_ptr<PolynomialChaosExpansion> const& part1,
                                                                           double w2,
                                                                           std::shared_ptr<PolynomialChaosExpansion> const& part2) const
{
  Eigen::VectorXd wts(2);
  wts << w1, w2;
  return PolynomialChaosExpansion::ComputeWeightedSum({part1,part2},wts);
}


double AdaptiveSmolyakPCE::ComputeMagnitude(std::shared_ptr<PolynomialChaosExpansion> const& estimate) const
{
  return estimate->Magnitude().norm();
}
