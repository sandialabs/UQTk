#ifndef SCALEVECTOR_H_
#define SCALEVECTOR_H_

#include "MUQ/Modeling/ModPiece.h"

namespace muq {
  namespace Modeling {
    class ScaleVector : public ModPiece {
    public:

      ScaleVector(double const scale, unsigned int const dim);

      virtual ~ScaleVector() = default;

    private:

      virtual void EvaluateImpl(ref_vector<Eigen::VectorXd> const& inputs) override;

      virtual void JacobianImpl(unsigned int const outwrt, unsigned int const inwrt, ref_vector<Eigen::VectorXd> const& inputs) override;

      virtual void GradientImpl(unsigned int const outwrt, unsigned int const inwrt, ref_vector<Eigen::VectorXd> const& inputs, Eigen::VectorXd const& sens) override;

      virtual void ApplyJacobianImpl(unsigned int const outwrt, unsigned int const inwrt, ref_vector<Eigen::VectorXd> const& inputs, Eigen::VectorXd const& targ) override;

      const double scale;
    };
  } // namespace Modeling
} // namespace muq

#endif
