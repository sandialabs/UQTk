#ifndef IMPORTANCESAMPLING_H_
#define IMPORTANCESAMPLING_H_

#include <boost/property_tree/ptree.hpp>

#include "MUQ/Modeling/Distributions/Distribution.h"
#include "MUQ/Modeling/Distributions/Density.h"

#include "MUQ/SamplingAlgorithms/SamplingAlgorithm.h"

namespace muq {
  namespace SamplingAlgorithms {
    /** @ingroup MCMC
      @class ImportanceSampling
      @brief Defines an imporance sampline sampler
      @details
      <B>Configuration Parameters:</B>
      Parameter Key | Type | Default Value | Description |
      ------------- | ------------- | ------------- | ------------- |
      "NumSamples"  | Int | - | The total number of steps (including burnin) to take, i.e., the length of the Markov chain. |
    */
    class ImportanceSampling : public SamplingAlgorithm {
    public:

      /**
      Use the target distribution as the baising distribution.  This is equivalent to Monte Carlo
      @param[in] target The target distribution
      @param[in] pt Options for the importance sampler
      */
      ImportanceSampling(std::shared_ptr<muq::Modeling::Distribution> const& target, boost::property_tree::ptree const& pt);

      /**
      @param[in] target The target distribution
      @param[in] bias The biasing distribution
      @param[in] pt Options for the importance sampler
      */
      ImportanceSampling(std::shared_ptr<muq::Modeling::ModPiece> const& target, std::shared_ptr<muq::Modeling::Distribution> const& bias, boost::property_tree::ptree const& pt);

      /**
      @param[in] target The target distribution
      @param[in] bias The biasing distribution
      @param[in] hyperparameters Hyperparameters for the biasing distribution
      @param[in] pt Options for the importance sampler
      */
      ImportanceSampling(std::shared_ptr<muq::Modeling::ModPiece> const& target, std::shared_ptr<muq::Modeling::Distribution> const& bias, std::vector<Eigen::VectorXd> hyperparameters, boost::property_tree::ptree const& pt);

      /**
      @param[in] target The target distribution
      @param[in] hyperparameters Hyperparameters for the biasing distribution
      @param[in] pt Options for the importance sampler
      */
      ImportanceSampling(std::shared_ptr<muq::Modeling::Distribution> const& target, std::vector<Eigen::VectorXd> hyperparameters, boost::property_tree::ptree const& pt);

      virtual ~ImportanceSampling() = default;

    private:

      virtual std::shared_ptr<SampleCollection> RunImpl(std::vector<Eigen::VectorXd> const& x0) override;

      /// The number of samples
      const unsigned int numSamps;

      // The target density
      std::shared_ptr<muq::Modeling::ModPiece> target;

      /// The biasing distribution
      std::shared_ptr<muq::Modeling::Distribution> bias;

      /// Hyperparameters for the biasing distribution
      const std::vector<Eigen::VectorXd> hyperparameters = std::vector<Eigen::VectorXd>();
    };
  } // namespace SamplingAlgorithms
} // namespace muq

#endif
