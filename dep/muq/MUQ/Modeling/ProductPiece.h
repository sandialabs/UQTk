#ifndef PRODUCTPIECE_H_
#define PRODUCTPIECE_H_

#include "MUQ/Modeling/ModPiece.h"

namespace muq {
  namespace Modeling {

    /**
      @brief Returns the components in a single vector input.
      @details return \f$y = \prod_{i=1}^N x_i\f$.
    */
    class ProductPiece: public ModPiece {
    public:

      ProductPiece(unsigned int vectorDim);

    private:

      virtual void EvaluateImpl(muq::Modeling::ref_vector<Eigen::VectorXd> const& inputs) override;

    }; // class ProductPiece

  } // namespace Modeling
} // namespace muq

#endif
