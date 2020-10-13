#ifndef SPLITVECTOR_H_
#define SPLITVECTOR_H_

#include "MUQ/Modeling/ModPiece.h"

namespace muq {
  namespace Modeling {
    class SplitVector : public ModPiece {
    public:
      /**
        @param[in] ind The first index of the segment for each output
        @param[in] size The size of each segment
        @param[in] insize The size of the input vector
      */
      SplitVector(Eigen::VectorXi const& ind, Eigen::VectorXi const& size, unsigned int const insize);

      virtual ~SplitVector() = default;

    private:

      virtual void EvaluateImpl(ref_vector<Eigen::VectorXd> const& inputs) override;

      virtual void JacobianImpl(unsigned int const outwrt, unsigned int const inwrt, ref_vector<Eigen::VectorXd> const& inputs) override;

      virtual void GradientImpl(unsigned int const outwrt, unsigned int const inwrt, ref_vector<Eigen::VectorXd> const& inputs, Eigen::VectorXd const& sens) override;

      virtual void ApplyJacobianImpl(unsigned int const outwrt, unsigned int const inwrt, ref_vector<Eigen::VectorXd> const& inputs, Eigen::VectorXd const& targ) override;

      const Eigen::VectorXi ind;

      const Eigen::VectorXi size;

    };
  } // namespace Modeling
} // namespace muq

#endif
