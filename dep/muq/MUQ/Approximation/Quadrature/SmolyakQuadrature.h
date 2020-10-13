#ifndef SMOLYAKQUADRATURE_H
#define SMOLYAKQUADRATURE_H

#include "MUQ/Approximation/Quadrature/Quadrature.h"
#include "MUQ/Approximation/Quadrature/FullTensorQuadrature.h"

#include "MUQ/Utilities/MultiIndices/MultiIndexSet.h"

namespace muq {
namespace Approximation {

  /** @class SmolyakQuadrature
      @ingroup Quadrature
      @brief Computes static Smolyak quadrature rules for multivariate integration.
  */
  class SmolyakQuadrature : public Quadrature
  {
  public:

    SmolyakQuadrature(unsigned int dim, std::shared_ptr<Quadrature> const& scalarRule);

    SmolyakQuadrature(std::vector<std::shared_ptr<Quadrature>> const& scalarRulesIn);

    virtual ~SmolyakQuadrature() = default;

    virtual void Compute(unsigned int order) override;

    virtual void Compute(Eigen::RowVectorXi const& orders) override;

    virtual void Compute(std::shared_ptr<muq::Utilities::MultiIndexSet> const& multis);

    /** The Smolyak rule is defined as the tensor product of differences
        of 1d quadrature rules
        \f[
        \sum_{\mathbf{k}\in\mathcal{K}} \Delta_{k_1}^{(1)} \otimes \cdots \otimes \Delta^{(d)}_{k_d},
        \f]
        where \f$\Delta_{n}^{(i)} = \mathscr{L}_{n}^{(i)} - \mathscr{L}_{n-1}^{(i)}\f$.
          Computationally however, it is easier to
        work directly with tensor products of the quadrature rules themselves.
        Rearranging the above tensor product, we have an expansion of the form
        \f[
        \sum_{\mathbf{k}\in\mathcal{K}} c_{\mathbf{k}} \mathscr{L}_{k_1}^{(i)} \otimes \mathscr{L}_{k_d}^{(d)}.
        \f]
        This function returns the weights \f$c_{\mathbf{k}}\f$ and is called
        during the call the SmolyakQuadrature::Compute().
    */
    static Eigen::VectorXd ComputeWeights(std::shared_ptr<muq::Utilities::MultiIndexSet> const& multis);

    static void UpdateWeights(unsigned int activeInd,
                              std::shared_ptr<muq::Utilities::MultiIndexSet> const& multis,
                              Eigen::VectorXd                                     & multiWeights);

    std::shared_ptr<muq::Utilities::MultiIndexSet> BuildMultis(Eigen::RowVectorXi const& orders) const;

  private:

    std::vector<std::shared_ptr<Quadrature>> scalarRules;

  }; // class SmolyakQuadrature
}
}


#endif // #ifndef SMOLYAKQUADRATURE_H
