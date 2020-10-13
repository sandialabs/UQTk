#ifndef REPLICATEOPERATOR_H_
#define REPLICATEOPERATOR_H_

#include "MUQ/Modeling/ModPiece.h"

namespace muq {
  namespace Modeling {

    /**
      @brief A muq::Modeling::ModPiece that replicates a single input vector \f$N\f$ times.
      @details Given a single input vector of length \f$M\f$, this ModPiece returns a
      vector of length \f$NM\f$ that contains \f$N\f$ replicates of the input vector.

    */
    class ReplicateOperator: public ModPiece {
    public:

      /**
	       @param[in] vectorDim The dimension of the input vector
         @param[in] numRepeat The number of times to replicate the input vector.
       */
      ReplicateOperator(unsigned int vectorDim, unsigned int numRepeat);


    private:
      unsigned int numRepl;
      
      virtual void EvaluateImpl(muq::Modeling::ref_vector<Eigen::VectorXd> const& inputs) override;

    }; // class ConstantVector

  } // namespace Modeling
} // namespace muq

#endif
