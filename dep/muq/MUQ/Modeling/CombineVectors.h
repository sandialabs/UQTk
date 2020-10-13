#ifndef COMBINEVECTORS_H_
#define COMBINEVECTORS_H_

#include "MUQ/Modeling/ModPiece.h"

namespace muq {
  namespace Modeling {
    class CombineVectors : public ModPiece {
    public:
      CombineVectors(Eigen::VectorXi const& inputSizes);

      virtual ~CombineVectors() = default;
    private:

      virtual void EvaluateImpl(ref_vector<Eigen::VectorXd> const& inputs) override;

      virtual void JacobianImpl(unsigned int const outwrt, unsigned int const inwrt, ref_vector<Eigen::VectorXd> const& inputs) override;

      virtual void GradientImpl(unsigned int const outwrt, unsigned int const inwrt, ref_vector<Eigen::VectorXd> const& inputs, Eigen::VectorXd const& sens) override;

      virtual void ApplyJacobianImpl(unsigned int const outwrt, unsigned int const inwrt, ref_vector<Eigen::VectorXd> const& inputs, Eigen::VectorXd const& targ) override;
    };
  } // namespace Modeling
} // namespace muq

#endif
