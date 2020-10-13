#ifndef SAMPLINGALGORITHM_H_
#define SAMPLINGALGORITHM_H_

#include "MUQ/config.h"

#if MUQ_HAS_PARCER
#include <parcer/Communicator.h>
#endif

#include "MUQ/Modeling/WorkPiece.h"

#include "MUQ/SamplingAlgorithms/AbstractSamplingProblem.h"
#include "MUQ/SamplingAlgorithms/SampleCollection.h"

/**
@defgroup SamplingAlgorithms

*/

namespace muq {
  namespace SamplingAlgorithms {

    class SamplingAlgorithm {//} : public muq::Modeling::WorkPiece {
    public:

      SamplingAlgorithm(std::shared_ptr<SampleCollection> const& samples);

      SamplingAlgorithm(std::shared_ptr<SampleCollection> const& samplesIn, std::shared_ptr<SampleCollection> const& QOIsIn);

#if MUQ_HAS_PARCER
      SamplingAlgorithm(std::shared_ptr<SampleCollection> const& samplesIn, std::shared_ptr<parcer::Communicator> const& comm);
#endif

      virtual ~SamplingAlgorithm() = default;

      virtual std::shared_ptr<SampleCollection> GetSamples() const;

      virtual std::shared_ptr<SampleCollection> GetQOIs() const;

      virtual void SetState(std::vector<Eigen::VectorXd> const& x0);

      virtual std::shared_ptr<SampleCollection> Run(std::vector<Eigen::VectorXd> const& x0 = std::vector<Eigen::VectorXd>());

      template<typename... Args>
        inline std::shared_ptr<SampleCollection> Run(Args const&... args) {
          std::vector<Eigen::VectorXd> vec;
          return RunRecurse(vec, args...);
        }

#if MUQ_HAS_PARCER
      std::shared_ptr<parcer::Communicator> GetCommunicator() const;
#endif

    protected:

      virtual std::shared_ptr<SampleCollection> RunImpl(std::vector<Eigen::VectorXd> const& x0) = 0;

      /**
	 Inputs:
	 <ol>
	 <li> Parameters for the algorithm
	 <li> The muq::SamplingAlgorithms::SamplingProblem that evaluates/samples the target distribution
	 </ol>
	 @param[in] inputs Inputs to the algorithm
       */
      //virtual void EvaluateImpl(muq::Modeling::ref_vector<boost::any> const& inputs) override;

      std::shared_ptr<SampleCollection> samples;

      std::shared_ptr<SampleCollection> QOIs;

#if MUQ_HAS_PARCER
      std::shared_ptr<parcer::Communicator> comm = nullptr;
#endif

    private:

      template<typename... Args>
        inline std::shared_ptr<SampleCollection> RunRecurse(std::vector<Eigen::VectorXd>& vec, Eigen::VectorXd const& ith, Args const&... args) {
          vec.push_back(ith);
          return RunRecurse(vec, args...);
        }

        inline std::shared_ptr<SampleCollection> RunRecurse(std::vector<Eigen::VectorXd>& vec, Eigen::VectorXd const& last) {
          vec.push_back(last);
          return Run(vec);
        }
    };
  } // namespace SamplingAlgorithms
} // namespace muq

#endif
