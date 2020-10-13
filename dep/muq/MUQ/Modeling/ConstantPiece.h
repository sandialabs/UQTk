#ifndef CONSTANTPIECE_H_
#define CONSTANTPIECE_H_

#include "MUQ/Modeling/WorkPiece.h"

namespace muq {
  namespace Modeling {
    /// A muq::Modeling::WorkPiece with no inputs and known, constant outputs
    class ConstantPiece : public WorkPiece {
    public:

      /// Create a muq::Modeling::ConstantPiece with the outputs given in a vector
      /**
	 @param[in] outs The outputs
       */
      ConstantPiece(std::vector<boost::any> const& outs);

      /// Create a muq::Modeling::ConstantPiece with no given outputs
      /**
	 This creates a WorkPiece that takes no inputs and returns no outputs.  The user can set the outputs to any thing after construction.  The number of outputs is set to be variable.
       */
      ConstantPiece();

      /// Create a muq::Modeling::ConstantPiece with the outputs given indivually
      /**
	 The number of outputs is -1 (variable) even though we know the outputs.
	 @param[in] in The \f$i^{th}]f$ input
	 @param[in] args The \f$[i+1, N\f$] outputs
       */
      template<typename ith, typename... Args>
	ConstantPiece(ith const& in, Args... args) : ConstantPiece(args...) {
	// add this output to the begining of the output vector
	outputs.insert(outputs.begin(), in);

	// set the number of outputs
	numOutputs = outputs.size();

	// set the output types
	outputTypes = Types(Types(outputs));
      }

      /// Create a muq::Modeling::ConstantPiece with the outputs given indivually
      /**
	 The number of outputs is -1 (variable) even though we know the outputs.
	 @param[in] last the last output
      */
      template<typename last>
	ConstantPiece(last const& in) : WorkPiece(0, WorkPiece::Fix::Inputs) {
	// the outputs will never change so we should not clear them
	clearOutputs = false;

	// add this output to the begining of the output vector
	outputs.insert(outputs.begin(), in);

	// set the number of outputs
	numOutputs = outputs.size();

	// set the output types
	outputTypes = Types(Types(outputs));
      }

      /// Set the outputs to be empty
      void SetOutputs();

      /// Set the outputs
      /**
	 @param[in] outs The new outputs
       */
      void SetOutputs(std::vector<boost::any> const& outs);

      /// set the outputs with the outputs given indivually
      /**
	 @param[in] in The \f$i^{th}]f$ input
	 @param[in] args The \f$[i+1, N\f$] outputs
       */
      template<typename ith, typename... Args>
	void SetOutputs(ith const& in, Args... args) {
	// create an empty vector of outputs
	std::vector<boost::any> outs(0);

	// add this output to the begining of the output vector
	outs.insert(outs.begin(), in);

	// call recursively
	SetOutputs(outs, args...);
      }

    private:

      /// Set the outputs with the outputs given indivually
      /**
	 @param[in] outs The outputs so far
	 @param[in] in The \f$i^{th}]f$ input
	 @param[in] args The \f$[i+1, N\f$] outputs
       */
      template<typename ith, typename... Args>
	void SetOutputs(std::vector<boost::any>& outs, ith const& in, Args... args) {
	// add this output
	outs.push_back(in);

	// call recursively
	SetOutputs(outs, args...);
      }

      /// Set the outputs with the outputs given indivually
      /**
	 @param[in] outs The outputs so far
	 @param[in] in The last input
       */
      template<typename last>
	void SetOutputs(std::vector<boost::any>& outs, last const& in) {
	// add this output
	outs.push_back(in);

	// call recursively
	SetOutputs(outs);
      }

      /// The outputs are already set and not cleared so don't do anything
      /**
	 @param[in] inputs An empty vector of inputs
       */
      virtual void EvaluateImpl(muq::Modeling::ref_vector<boost::any> const& inputs) override;

    };
  } // namespace Modeling
} // namespace muq

#endif
