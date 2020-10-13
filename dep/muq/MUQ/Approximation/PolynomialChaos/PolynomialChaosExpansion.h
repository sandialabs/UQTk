#ifndef POLYNOMIALCHAOSEXPANSION_H_
#define POLYNOMIALCHAOSEXPANSION_H_

#include "MUQ/Approximation/Polynomials/BasisExpansion.h"
#include "MUQ/Approximation/Polynomials/OrthogonalPolynomial.h"

#include <vector>
#include <set>

namespace muq {
namespace Approximation {

  class PCEFactory;

  /**
   @class PolynomialChaosExpansion
   @ingroup PolynomialChaos
   @brief A class for representing and using expansions of orthogonal multivariate polynomials
   @details
   * A particular polynomial chaos expansion for a function from R^n->R^m. This class uses
   * some the MultiIndexSet class to define PCE terms that are in the expansion.
   * For each PCE term and output, there is a coefficient. The PCE is built from 1D polynomials
   * specified so each dimension may have independent polynomial
   * choices.
   *
  @see BasisExpansion, PCEFactory
   */
  class PolynomialChaosExpansion : public BasisExpansion {
  friend class PCEFactory;

  public:

    PolynomialChaosExpansion(std::shared_ptr<OrthogonalPolynomial>          const& basisCompsIn,
                             std::shared_ptr<muq::Utilities::MultiIndexSet> const& multisIn,
                             Eigen::MatrixXd                                const& coeffsIn);

    PolynomialChaosExpansion(std::shared_ptr<OrthogonalPolynomial>          const& basisCompsIn,
                             std::shared_ptr<muq::Utilities::MultiIndexSet> const& multisIn,
                             unsigned int                                          outputDim);

    PolynomialChaosExpansion(std::vector<std::shared_ptr<IndexedScalarBasis>> const& basisCompsIn,
                             std::shared_ptr<muq::Utilities::MultiIndexSet>   const& multisIn,
                             Eigen::MatrixXd                                  const& coeffsIn);

    PolynomialChaosExpansion(std::vector<std::shared_ptr<IndexedScalarBasis>> const& basisCompsIn,
                             std::shared_ptr<muq::Utilities::MultiIndexSet>   const& multisIn,
                             unsigned int                                            outputDim);


    virtual ~PolynomialChaosExpansion() = default;

    ///compute the variance of the current expansion
    Eigen::VectorXd Variance() const;
    Eigen::MatrixXd Covariance() const;
    Eigen::VectorXd Mean() const;

    ///Compute the L2 norm of each output.
    Eigen::VectorXd Magnitude() const;

    /**
     * Compute the weighted sum of polynomial expansions. Slow, because it assumes you haven't been
     * tracking which polynomials to keep, so it does that, then calls the private method. If you
     * know already, your class should be a friend and use the private method directly, as that
     * will be faster.
     * */
    static std::shared_ptr<PolynomialChaosExpansion> ComputeWeightedSum(std::vector<std::shared_ptr<PolynomialChaosExpansion>> expansions,
                                                                        Eigen::VectorXd                                 const& weights);


    ///Compute the Sobol total sensitivity index for the input dimension, for each output dimension
    Eigen::VectorXd TotalSensitivity(unsigned int targetDim) const;

    ///Compute all Sobol total sensitivities. Rows are outputs, each column is an input.
    Eigen::MatrixXd TotalSensitivity() const;

    ///Compute the main sensitivity index for the input dimension, for each output dimension
    Eigen::VectorXd SobolSensitivity(unsigned int targetDim) const;

    /** Computes the Sobol sensitivity for a group of input parameters.
    */
    Eigen::VectorXd SobolSensitivity(std::vector<unsigned int> const& targetDims) const;

    /** Compute all the main sensitivities (e.g., one term Sobol sensitivity indices).
        Rows are outputs, each column is an input.
    */
    Eigen::MatrixXd MainSensitivity() const;

  private:

    /**
     * This function returns the sqrt of the normalization, sqrt(<p*p>),
     * for each PCE basis function p.
     * NB: This function does not change depending on the coefficients.
     */
    Eigen::VectorXd GetNormalizationVec() const;

    /**
     * An internal function to compute the weighted sum of polynomial expansions if you already know which polynomials are
     * included, where polynomials
     * is a clean copy not used by other expansions. Actually does the addition.
     * */
    static std::shared_ptr<PolynomialChaosExpansion> ComputeWeightedSum(std::vector<std::shared_ptr<PolynomialChaosExpansion>> expansions,
                                                                        Eigen::VectorXd const&                                 weights,
                                                                        std::shared_ptr<muq::Utilities::MultiIndexSet> const&  polynomials,
                                                                        std::vector<std::vector<unsigned int>>         const&  locToGlob);

  };

} // namespace muq
} // namespace Approximation



#endif /* POLYNOMIALCHAOSEXPANSION_H_ */
