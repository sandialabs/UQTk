#ifndef CONSTANTVECTOR_H_
#define CONSTANTVECTOR_H_

#include "MUQ/Modeling/ModPiece.h"

namespace muq {
  namespace Modeling {

    /// A muq::Modeling::ModPiece with no inputs and a single known vector-valued output
    class ConstantVector : public ModPiece {
    public:

      /// Create a muq::Modeling::ConstantPiece with the outputs given in a vector of vectors
      /**
      @param[in] outs The outputs
       */
      ConstantVector(std::vector<Eigen::VectorXd> const& outs);

      /// Create a muq::Modeling::ConstantPiece with the outputs given in a vector
      /**
	       @param[in] outs The outputs
       */
      ConstantVector(Eigen::VectorXd const& valIn);

      /// Set the outputs
      /**
	       @param[in] outs The new outputs
       */
      void SetValue(Eigen::VectorXd const& valIn);

    private:
      /// The outputs are already set and not cleared so don't do anything
      /**
	       @param[in] inputs An empty vector of inputs
      */
      virtual void EvaluateImpl(muq::Modeling::ref_vector<Eigen::VectorXd> const& inputs) override;

      static Eigen::VectorXi OutSizes(std::vector<Eigen::VectorXd> const& outs);

    }; // class ConstantVector

  } // namespace Modeling
} // namespace muq

#endif
