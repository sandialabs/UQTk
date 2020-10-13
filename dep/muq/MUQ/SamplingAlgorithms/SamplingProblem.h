#ifndef SAMPLINGPROBLEM_H_
#define SAMPLINGPROBLEM_H_

// include Density and not ModPiece so that if a SamplingProblem is constructed with a Density the compiler knows it is a child of ModPiece
#include "MUQ/Modeling/Distributions/Density.h"
#include "MUQ/SamplingAlgorithms/AbstractSamplingProblem.h"

namespace muq {
  namespace SamplingAlgorithms {

    /**
    @ingroup SamplingAlgorithms
    @class SamplingProblem
    @brief Class for sampling problems based purely on a density function.
    */
    class SamplingProblem : public AbstractSamplingProblem{
    public:

      /**
	     @param[in] target The target distribution
       */
      SamplingProblem(std::shared_ptr<muq::Modeling::ModPiece> const& targetIn);

      virtual ~SamplingProblem() = default;


      virtual double LogDensity(unsigned int const t, std::shared_ptr<SamplingState> const& state, AbstractSamplingProblem::SampleType type) override;

      virtual Eigen::VectorXd GradLogDensity(std::shared_ptr<SamplingState> const& state,
                                             unsigned                       blockWrt);


      std::shared_ptr<muq::Modeling::ModPiece> GetDistribution(){return target;};

    protected:

      /// The target distribution (the prior in the inference case)
      std::shared_ptr<muq::Modeling::ModPiece> target;

    private:

      static unsigned GetNumBlocks(std::shared_ptr<muq::Modeling::ModPiece> const& target);
      static std::vector<int> GetBlockSizes(std::shared_ptr<muq::Modeling::ModPiece> const& target);

    };
  } // namespace SamplingAlgorithms
} // namespace muq

#endif
