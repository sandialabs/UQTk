#ifndef RANDOMVARIABLE_H
#define RANDOMVARIABLE_H

#include "MUQ/Modeling/Distributions/Distribution.h"
#include "MUQ/Modeling/ModPiece.h"

namespace muq{
  namespace Modeling{

    class RandomVariable : public Distribution, public ModPiece{

    public:
      RandomVariable(std::shared_ptr<Distribution> distIn);

      virtual ~RandomVariable() = default;

    protected:
      std::shared_ptr<Distribution> dist;

      virtual void EvaluateImpl(ref_vector<Eigen::VectorXd> const& inputs) override;

      virtual double LogDensityImpl(ref_vector<Eigen::VectorXd> const& inputs) override;
      virtual Eigen::VectorXd GradLogDensityImpl(unsigned int wrt, ref_vector<Eigen::VectorXd> const& inputs) override;
      virtual Eigen::VectorXd SampleImpl(ref_vector<Eigen::VectorXd> const& inputs) override;

    }; // class RandomVariable

  } // namespace Modeling
} // namespace muq



#endif // #ifndef RANDOMVARIABLE_H
