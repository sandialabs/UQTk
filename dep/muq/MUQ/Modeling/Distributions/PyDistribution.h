#ifndef PYDISTRIBUTION_H_
#define PYDISTRIBUTION_H_

#include "MUQ/Modeling/Distributions/Distribution.h"
#include "MUQ/Modeling/Distributions/GaussianBase.h"


namespace muq {
  namespace Modeling {
    class PyDistribution : public Distribution {
    public:
      PyDistribution(unsigned int varSizeIn, Eigen::VectorXi const& hyperSizesIn = Eigen::VectorXi());

      static std::vector<Eigen::VectorXd> ToStdVec(ref_vector<Eigen::VectorXd> const& input);
      static ref_vector<Eigen::VectorXd> ToRefVec(std::vector<Eigen::VectorXd> const& input);

    protected:

      virtual Eigen::VectorXd SampleImpl(std::vector<Eigen::VectorXd> const& inputs) = 0;
      virtual Eigen::VectorXd SampleImpl(ref_vector<Eigen::VectorXd> const& inputs) override;


      virtual double LogDensityImpl(ref_vector<Eigen::VectorXd> const& inputs) override;

      virtual double LogDensityImpl(std::vector<Eigen::VectorXd> const& inputs) = 0;
    };


    class PyGaussianBase : public GaussianBase {
    public:

      using GaussianBase::GaussianBase;

      virtual ~PyGaussianBase() = default;

      // This is needed to make SampleImpl public in PyGaussianBase
      using GaussianBase::SampleImpl;

    };

  } // namespace Modeling
} // namespace muq

#endif
