#ifndef IDENTITYPIECE_H_
#define IDENTITYPIECE_H_

#include "MUQ/Modeling/WorkPiece.h"

namespace muq {
  namespace Modeling {
    /// A muq::Modeling::WorkPiece that returns its inputs (identity operator)
    class IdentityPiece : public WorkPiece {
    public:
      
      /// Create a muq::Modeling::IdentityPiece with no fixed number of inputs and outputs and variable input/output types.
      IdentityPiece();
      
      /// Create a muq::Modeling::IdentityPiece with fixed number of inputs/outputs and variable input/output types.
      /**
	 @param[in] num The number of inputs/outputs (both must be fixed)
      */
      IdentityPiece(int const num);
      
      /// Create a muq::Modeling::IdentityPiece with a fixed number of inputs/outputs with specified types 
      /**
	 @param[in] types A vector of strings, each element is the type of an input or output (the number of inputs or outputs is the size of this vector)
      */
      IdentityPiece(std::vector<std::string> const& types);

      /// Create a muq::Modeling::IdentityPiece where some of the inputs/outputs have specified types 
      /**
	 @param[in] types A map from the input/output number to the input/output type
      */
      IdentityPiece(std::map<unsigned int, std::string> const& types);

      /// Create a muq::Modeling::IdentityPiece where some of the inputs/outputs have specified types the number of inputs/outputs is fixed
      /**
	 @param[in] types A map from the input/output number to the input/output type
	 @param[in] num The number of inputs/outputs
       */
      IdentityPiece(std::map<unsigned int, std::string> const& types, unsigned int const num);

    private:

      /**
	 @param[in] inputs The inputs to this WorkPiece
       */
      virtual void EvaluateImpl(muq::Modeling::ref_vector<boost::any> const& inputs) override;
      
    };
  } // namespace Modeling
} // namespace muq

#endif
