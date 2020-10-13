#ifndef WORKPIECE_H_
#define WORKPIECE_H_

#include<iostream>
#include<vector>
#include<map>
#include<string>
#include<sstream>
#include<cassert>
#include<memory>
#include<string>

#include "boost/any.hpp"
#include "boost/optional.hpp"

#include <Eigen/Dense>

#include "MUQ/config.h"

#include "MUQ/Modeling/WorkPiece.h"

#if MUQ_HAS_SUNDIALS==1
// Sundials includes
#include <nvector/nvector_serial.h>
#include <sundials/sundials_dense.h> // definitions DlsMat DENSE_ELEM
#endif

namespace muq {
  namespace Modeling {

    /**
      @defgroup Modeling
      @brief Tools for constructing coupled physical and statistical models.
      @details

    */


    // Forward declaration of WorkGraphPiece
    class WorkGraphPiece;
    class WorkGraph;

    /// A vector of references to something ...
    template <typename T>
      using ref_vector = std::vector<std::reference_wrapper<const T>>;

    /// Base class for MUQ's modelling envronment
    class WorkPiece {

      // Make WorkGraphPiece a friend so it can access WorkPiece's outputs directly
      friend class WorkGraphPiece;
      friend class WorkGraph;

    protected:

      /// Does the constructor fix the inputs or the outputs?
      enum Fix {
	/// The constructor fixes the input number and possibly the types
	Inputs,
	/// The constructor fixes the output number and possibly the types
	Outputs };

    public:

      /// Create a muq::Modeling::WorkPiece with no fixed number of inputs and outputs and variable input/output types.
      WorkPiece();

      /// Create a muq::Modeling::WorkPiece with either a fixed number of inputs or outputs and variable input/output types.
      /**
	 @param[in] num The number of inputs or outputs (which one depends on the second parameter)
	 @param[in] fix WorkPiece::Fix::Inputs (default): the first parameter is the number of inputs; WorkPiece::Fix::Outputs: the first parameter is the number of outputs
      */
      WorkPiece(int const num, WorkPiece::Fix const fix = WorkPiece::Fix::Inputs);

      /// Create a muq::Modeling::WorkPiece with a fixed number of inputs and outputs but variable input/output types
      /**
	 @param[in] numIns The number of inputs
	 @param[in] numOuts The number of outputs
      */
      WorkPiece(int const numIns, int const numOuts);

      /// Create a muq::Modeling::WorkPiece with either a fixed number of inputs with specified types or a fixed number of outputs with specified types
      /**
	 If the number and type of the inputs is specified then the number and type of the outputs is variable.  The opposite is true if the number and type of the outputs is specified.
	 @param[in] types A vector of strings, each element is the type of an input or output (the number of inputs or outputs is the size of this vector)
	 @param[in] fix WorkPiece::Fix::Inputs (default): the elements of the first parameter are the types of the inputs; WorkPiece::Fix::Outputs: the elements of the first parameter are the types of the outputs
      */
      WorkPiece(std::vector<std::string> const& types, WorkPiece::Fix const fix = WorkPiece::Fix::Inputs);

      /// Create a muq::Modeling::WorkPiece where either some of the inputs have specified types or some of the outputs have specified types
      /**
	 If the inputs are specified, then the outputs are not (and vice versa).  The number of in/outputs is variable but some of them have specified type.  For example, if the first input is a string and the third input is a double then
	 <ul>
	 <li> Evaluate("string", 1.0, 2.0);
	 <li> Evaluate("string");
	 <li> Evaluate();
	 <li> Evaluate("string", "another string", 5.0, std::shared_ptr<AnObject>);
	 </ul>
	 are valid muq::Modeling::WorkPiece::Evaluate calls but
	 <ul>
	 <li> Evaluate(1.0, 2.0, "string");
	 </ul>
	 is not valid.
	 @param[in] types A map from the input/output number to the input/output type
	 @param[in] fix WorkPiece::Fix::Inputs (default): the elements of the first parameter are the types of the inputs; WorkPiece::Fix::Outputs: the elements of the first parameter are the types of the outputs
      */
      WorkPiece(std::map<unsigned int, std::string> const& types, WorkPiece::Fix const fix = WorkPiece::Fix::Inputs);

      /// Create a muq::Modeling::WorkPiece where either some of the inputs have specified types or some of the outputs have specified types and either the number of inputs or the number of outputs is fixed
      /**
	 @param[in] types A map from the input/output number to the input/output type
	 @param[in] num The number of inputs/outputs
	 @param[in] fixTypes WorkPiece::Fix::Inputs (default): the elements of the first parameter are the types of the inputs; WorkPiece::Fix::Outputs: the elements of the first parameter are the types of the outputs
	 @param[in] fixNum WorkPiece::Fix::Inputs (default): the second parameter is the number of inputs; WorkPiece::Fix::Outputs: the second parameter is the number of outputs
       */
      WorkPiece(std::map<unsigned int, std::string> const& types, int const num, WorkPiece::Fix const fixTypes = WorkPiece::Fix::Inputs, WorkPiece::Fix const fixNum = WorkPiece::Fix::Inputs);

      /// Create a muq::Modeling::WorkPiece with a fixed number of inputs with specified types and a fixed number of outputs (of uknown type)
      /**
	 @param[in] types A vector of strings, each element is the type of an input (the number of inputs is the size of this vector)
	 @param[in] num The number of outputs
      */
      WorkPiece(std::vector<std::string> const& types, int const num);

      /// Create a muq::Modeling::WorkPiece with a fixed number of outputs with specified types and a fixed number of inputs (of uknown type)
      /**
	 @param[in] num The number of inputs
	 @param[in] types A vector of strings, each element is the type of an output (the number of outputs is the size of this vector)
      */
      WorkPiece(int const num, std::vector<std::string> const& types);

      /// Create a muq::Modeling::WorkPiece where some of the inputs are known and we know the input and output numbers
      /**
	 @param[in] inTypes A map from the input number to the input type
	 @param[in] numIns The number of inputs
	 @param[in] numOuts The number of outputs
      */
      WorkPiece(std::map<unsigned int, std::string> const& inTypes, int const numIns, int const numOuts);

      /// Create a muq::Modeling::WorkPiece where some of the outputs are known and we know the input and output numbers
      /**
	 @param[in] numIns The number of inputs
	 @param[in] outTypes A map from the input number to the input type
	 @param[in] numOuts The number of outputs
      */
      WorkPiece(int const numIns, std::map<unsigned int, std::string> const& outTypes, int const numOuts);

      /// Create a muq::Modeling::WorkPiece with a fixed number of inputs and outputs with specified types
      /**
	 @param[in] inTypes A vector of strings, each element is the type of an input (the number of inputs is the size of this vector)
	 @param[in] outTypes A vector of strings, each element is the type of an output (the number of outputs is the size of this vector)
      */
      WorkPiece(std::vector<std::string> const& inTypes, std::vector<std::string> const& outTypes);

      /// Create a muq::Modeling::WorkPiece where some of the inputs are known and all of the outputs have specified types
      /**
	 @param[in] inTypes A map from the input number to the input type
	 @param[in] outTypes A vector of strings, each element is the type of an output (the number of outputs is the size of this vector)
      */
      WorkPiece(std::map<unsigned int, std::string> const& inTypes, std::vector<std::string> const& outTypes);

      /// Create a muq::Modeling::WorkPiece where some of the inputs are known with a known number of inputs and all of the outputs have specified types
      /**
	 @param[in] inTypes A map from the input number to the input type
	 @param[in] num The number of inputs
	 @param[in] outTypes A vector of strings, each element is the type of an output (the number of outputs is the size of this vector)
      */
      WorkPiece(std::map<unsigned int, std::string> const& inTypes, int const num, std::vector<std::string> const& outTypes);

      /// Create a muq::Modeling::WorkPiece where some of the outputs and all of the inputs have specified types
      /**
	 @param[in] inTypes A vector of strings, each element is the type of an input (the number of inputs is the size of this vector)
	 @param[in] outTypes A map from the output number to the output type
      */
      WorkPiece(std::vector<std::string> const& inTypes, std::map<unsigned int, std::string> const& outTypes);

      /// Create a muq::Modeling::WorkPiece where some of the outputs with a known number of outputs and all of the inputs have specified types
      /**
	 @param[in] inTypes A vector of strings, each element is the type of an input (the number of inputs is the size of this vector)
	 @param[in] outTypes A map from the output number to the output type
	 @param[in] num The number of outputs
      */
      WorkPiece(std::vector<std::string> const& inTypes, std::map<unsigned int, std::string> const& outTypes, int const num);

      /// Create a muq::Mdoeling::WorkPiece where some of the inputs and some of the outputs have specified types
      /**
	 @param[in] inTypes A map from the input number to the input type
	 @param[in] outTypes A map from the output number to the output type
      */
      WorkPiece(std::map<unsigned int, std::string> const& inTypes, std::map<unsigned int, std::string> const& outTypes);

      /// Create a muq::Mdoeling::WorkPiece where some of the inputs and some of the outputs have specified types with a fixed number of inputs
      /**
	 @param[in] inTypes A map from the input number to the input type
	 @param[in] numIn The number of inputs
	 @param[in] outTypes A map from the output number to the output type
      */
      WorkPiece(std::map<unsigned int, std::string> const& inTypes, int const numIn, std::map<unsigned int, std::string> const& outTypes);

      /// Create a muq::Mdoeling::WorkPiece where some of the inputs and some of the outputs have specified types with a fixed number of outputs
      /**
	 @param[in] inTypes A map from the input number to the input type
	 @param[in] outTypes A map from the output number to the output type
	 @param[in] outNum The number of outputs
      */
      WorkPiece(std::map<unsigned int, std::string> const& inTypes, std::map<unsigned int, std::string> const& outTypes, int const numOut);

      /// Create a muq::Mdoeling::WorkPiece where some of the inputs and some of the outputs have specified types with a fixed number of inputs and outputs
      /**
	 @param[in] inTypes A map from the input number to the input type
	 @param[in] numIn The number of inputs
	 @param[in] outTypes A map from the output number to the output type
	 @param[in] outNum The number of outputs
      */
      WorkPiece(std::map<unsigned int, std::string> const& inTypes, int const numIn, std::map<unsigned int, std::string> const& outTypes, int const numOut);

      /// Default destructor
      virtual ~WorkPiece() {}

      /// Evaluate this muq::Modeling::WorkPiece
      /**
	 This function takes the inputs to the muq::Modeling::WorkPiece, which must match WorkPiece::numInputs and WorkPiece::inputTypes if they are specified.  It then calls WorkPiece::EvaluateImpl(), which populates WorkPiece::outputs using the input arguments (passed to WorkPiece::EvaluateImpl()).  This function then checks WorkPiece::outputs, which much match WorkPiece::numOutputs and WorkPiece::outputTypes if they are specified.

	 This function builds a reference vector to the inputs and calls WorkPiece::Evaluate(ref_vector<boost::any> const& ins)
	 @param[in] ins A vector of inputs
	 \return The outputs of this muq::Modeling::WorkPiece
      */
      std::vector<boost::any> const& Evaluate(std::vector<boost::any> const& ins);

      /// Evaluate this muq::Modeling::WorkPiece using references to the inputs
      /**
	 This function takes the references to the inputs to the muq::Modeling::WorkPiece, which must match WorkPiece::numInputs and WorkPiece::inputTypes if they are specified.  It then calls WorkPiece::EvaluateImpl(), which populates WorkPiece::outputs using the input arguments (passed to WorkPiece::EvaluateImpl()).  This function then checks WorkPiece::outputs, which much match WorkPiece::numOutputs and WorkPiece::outputTypes if they are specified.

	 References are used for efficiency in the muq::Modeling::WorkGraph.
	 @param[in] ins A vector of references to the inputs
	 \return The outputs of this muq::Modeling::WorkPiece
      */
      std::vector<boost::any> const& Evaluate(ref_vector<boost::any> const& ins);

      /// Evaluate this muq::Modeling::WorkPiece in the case that there are no inputs
      /**
	 \return The outputs of this muq::Modeling::WorkPiece
      */
      std::vector<boost::any> const& Evaluate();

      /// Evalaute this muq::Modeling::WorkPiece using multiple arguments
      /**
	 This function allows the user to call WorkPiece::Evaluate without first creating a vector of inputs.  Instead, the user calls WorkPiece::Evaluate with multiple arguments (if specified, the number of arguments must match the number of inputs) and this function creates the input vector.
	 @param[in] args The inputs (may be more than one)
	 \return The outputs of this muq::Modeling::WorkPiece
      */
      template<typename... Args>
	std::vector<boost::any> const& Evaluate(Args... args) {
	// clear the outputs
	Clear();

	// create the reference input vector
	ref_vector<boost::any> inputs;
	inputs.reserve(numInputs<0? 0 : numInputs);

	// begin calling the EvaluateRecursive with the first input
	return EvaluateRecursive(inputs, args...);
      }

  //     /// Evaluate the Jacobian of this muq::Modeling::WorkPiece using references to the inputs
  //     /**
	//  This function takes a vector of inputs to the muq::Modeling::WorkPiece, which must match WorkPiece::numInputs and WorkPiece::inputTypes if they are specified.  It then calls WorkPiece::JacobianImpl(), which computes the Jacobian.  The Jacobian must be implemented by a child class, unless both the input and the output are Eigen::VectorXd's.  In this case, we default to finite difference.
	//  @param[in] wrtIn The input number we are taking the Jacobian with respect to
	//  @param[in] wrtOut The output number we are taking the Jacobian with respect to
	//  @param[in] ins A vector of inputs
	//  \return The Jacobian of this muq::Modeling::WorkPiece
  //     */
  //     boost::any Jacobian(unsigned int const wrtIn, unsigned int const wrtOut, std::vector<boost::any> const& ins);
  //
  //     /// Evaluate the Jacobian of this muq::Modeling::WorkPiece using references to the inputs
  //     /**
	//  This function takes the references to the inputs to the muq::Modeling::WorkPiece, which must match WorkPiece::numInputs and WorkPiece::inputTypes if they are specified.  It then calls WorkPiece::JacobianImpl(), which computes the Jacobian.  The Jacobian must be implemented by a child class, unless both the input and the output are Eigen::VectorXd's.  In this case, we default to finite difference.
  //
	//  References are used for efficiency in the muq::Modeling::WorkGraph.
	//  @param[in] wrtIn The input number we are taking the Jacobian with respect to
	//  @param[in] wrtOut The output number we are taking the Jacobian with respect to
	//  @param[in] ins A vector of references to the inputs
	//  \return The Jacobian of this muq::Modeling::WorkPiece
  //     */
  //     boost::any Jacobian(unsigned int const wrtIn, unsigned int const wrtOut, ref_vector<boost::any> const& ins);
  //
  //     /// Evaluate the Jacobian of this muq::Modeling::WorkPiece using multiple arguments
  //     /**
	//  This function allows the user to compute the Jacobian of an output with respect to one of the inputs.  If both the input and the output type is an Eigen::VectorXd and the user has not implemented the Jacobian, we default to finite difference.
	//  @param[in] wrtIn The input number we are taking the Jacobian with respect to
	//  @param[in] wrtOut The output number we are taking the Jacobian with respect to
	//  @param[in] args The inputs (may be more than one)
	//  \return The Jacobian of this muq::Modeling::WorkPiece
  //      */
  //     template<typename... Args>
	// boost::any Jacobian(unsigned int const wrtIn, unsigned int const wrtOut, Args... args) {
	// // make sure the input and output number are valid
	// assert(numInputs<0 || wrtIn<numInputs);
	// assert(numOutputs<0 || wrtOut<numOutputs);
  //
	// // clear the outputs and derivative information
	// ClearDerivatives();
  //
	// // create the reference input vector
	// ref_vector<boost::any> inputs;
	// inputs.reserve(numInputs<0? 0 : numInputs);
  //
	// // begin calling the JacobianRecursive with the first input
	// return JacobianRecursive(wrtIn, wrtOut, inputs, args...);
  //     }
  //
  //     /// Compute the Jacobian using finite differences
  //     /**
	//  Assume the input and output type are Eigen::VectorXd.
	//  @param[in] wrtIn The input number we are taking the Jacobian with respect to
	//  @param[in] wrtOut The output number we are taking the Jacobian with respect to
	//  @param[in] args The inputs (may be more than one)
	//  @param[in] refTol Scaled value for the finite difference step size (defaults to 1e-4)
	//  @param[in] minTol Minimum value for the finite difference step size (defaults to 1e-6)
  //      */
  //     void JacobianByFD(unsigned int const wrtIn, unsigned int const wrtOut, ref_vector<boost::any> const& inputs, double const relTol = 1.0e-4, const double minTol = 1.0e-6);
  //
  //     /// Evaluate the action of the Jacobian of this muq::Modeling::WorkPiece using references to the inputs
  //     /**
	//  This function takes a vector of inputs to the muq::Modeling::WorkPiece, which must match WorkPiece::numInputs and WorkPiece::inputTypes if they are specified.  It then calls WorkPiece::JacobianActionImpl(), which computes the action of the Jacobian.  The Jacobian must be implemented by a child class, unless the input and the output and the vector the Jacobian is action on are Eigen::VectorXd's.  In this case, we call WorkPiece::Jacobian() and apply it to the vector.
	//  @param[in] wrtIn The input number we are taking the Jacobian with respect to
	//  @param[in] wrtOut The output number we are taking the Jacobian with respect to
	//  @param[in] vec We want to apply the Jacobian to this vector
	//  @param[in] ins A vector of inputs
	//  \return The Jacobian of this muq::Modeling::WorkPiece
  //     */
  //     boost::any JacobianAction(unsigned int const wrtIn, unsigned int const wrtOut, boost::any const& vec, std::vector<boost::any> const& ins);
  //
  //     /// Evaluate the action of the Jacobian of this muq::Modeling::WorkPiece using references to the inputs
  //     /**
	//  References are used for efficiency in the muq::Modeling::WorkGraph.
	//  @param[in] wrtIn The input number we are taking the Jacobian with respect to
	//  @param[in] wrtOut The output number we are taking the Jacobian with respect to
	//  @param[in] vec We want to apply the Jacobian to this vector
	//  @param[in] ins A vector of references to the inputs
	//  \return The Jacobian of this muq::Modeling::WorkPiece
  //     */
  //     boost::any JacobianAction(unsigned int const wrtIn, unsigned int const wrtOut, boost::any const& vec, ref_vector<boost::any> const& ins);
  //
  //     /// Evaluate the action of the Jacobian of this muq::Modeling::WorkPiece using multiple arguments
  //     /**
	//  This function allows the user to compute the action of the Jacobian of an output with respect to one of the inputs.  If the vector given to this function, the input type and the output type is an Eigen::VectorXd and the user has not implemented the Jacobian action, we call muq::Modeling::Jacobian and apply it to the given vector
	//  @param[in] wrtIn The input number we are taking the Jacobian with respect to
	//  @param[in] wrtOut The output number we are taking the Jacobian with respect to
	//  @param[in] vec We want to apply the Jacobian to this vector
	//  @param[in] args The inputs (may be more than one)
	//  \return The action of the Jacobian of this muq::Modeling::WorkPiece on vec
  //      */
  //     template<typename... Args>
	// boost::any JacobianAction(unsigned int const wrtIn, unsigned int const wrtOut, boost::any const& vec, Args... args) {
	// // make sure the input and output number are valid
	// assert(numInputs<0 || wrtIn<numInputs);
	// assert(numOutputs<0 || wrtOut<numOutputs);
  //
	// // clear the outputs and derivative information
	// ClearDerivatives();
  //
	// // create the reference input vector
	// ref_vector<boost::any> inputs;
	// inputs.reserve(numInputs<0? 0 : numInputs);
  //
	// // begin calling the JacobianActionRecursive with the first input
	// return JacobianActionRecursive(wrtIn, wrtOut, vec, inputs, args...);
  //     }
  //
  //     /// Compute the action of the Jacobian using finite differences
  //     /**
	//  Assume the input and output type are Eigen::VectorXd.
	//  @param[in] wrtIn The input number we are taking the Jacobian with respect to
	//  @param[in] wrtOut The output number we are taking the Jacobian with respect to
	//  @param[in] vec We want to apply the Jacobian to this vector
	//  @param[in] args The inputs (may be more than one)
  //      */
  //     void JacobianActionByFD(unsigned int const wrtIn, unsigned int const wrtOut, boost::any const& vec, ref_vector<boost::any> const& input);
  //
  //     /// Evaluate the action of the Jacobian transpose of this muq::Modeling::WorkPiece using references to the inputs
  //     /**
	//  This function takes a vector of inputs to the muq::Modeling::WorkPiece, which must match WorkPiece::numInputs and WorkPiece::inputTypes if they are specified.  It then calls WorkPiece::JacobianTransposeActionImpl(), which computes the action of the Jacobian transpose.  The Jacobian must be implemented by a child class, unless the input and the output and the vector the Jacobian transpose is action on are Eigen::VectorXd's.  In this case, we call WorkPiece::Jacobian() and apply it to the vector.
	//  @param[in] wrtIn The input number we are taking the Jacobian with respect to
	//  @param[in] wrtOut The output number we are taking the Jacobian with respect to
	//  @param[in] vec We want to apply the Jacobian to this vector
	//  @param[in] ins A vector of inputs
	//  \return The Jacobian of this muq::Modeling::WorkPiece
  //     */
  //     boost::any JacobianTransposeAction(unsigned int const wrtIn, unsigned int const wrtOut, boost::any const& vec, std::vector<boost::any> const& ins);
  //
  //     /// Evaluate the action of the Jacobian transpose of this muq::Modeling::WorkPiece using references to the inputs
  //     /**
	//  References are used for efficiency in the muq::Modeling::WorkGraph.
	//  @param[in] wrtIn The input number we are taking the Jacobian with respect to
	//  @param[in] wrtOut The output number we are taking the Jacobian with respect to
	//  @param[in] vec We want to apply the Jacobian to this vector
	//  @param[in] ins A vector of references to the inputs
	//  \return The Jacobian of this muq::Modeling::WorkPiece
  //     */
  //     boost::any JacobianTransposeAction(unsigned int const wrtIn, unsigned int const wrtOut, boost::any const& vec, ref_vector<boost::any> const& ins);
  //
  //     /// Evaluate the action of the Jacobian transpose of this muq::Modeling::WorkPiece using multiple arguments
  //     /**
	//  This function allows the user to compute the action of the Jacobian transpose of an output with respect to one of the inputs.  If the vector given to this function, the input type and the output type is an Eigen::VectorXd and the user has not implemented the Jacobian transpose action, we call muq::Modeling::Jacobian and apply it to the given vector
	//  @param[in] wrtIn The input number we are taking the Jacobian with respect to
	//  @param[in] wrtOut The output number we are taking the Jacobian with respect to
	//  @param[in] vec We want to apply the Jacobian transpose to this vector
	//  @param[in] args The inputs (may be more than one)
	//  \return The action of the Jacobian transpose of this muq::Modeling::WorkPiece on vec
  //      */
  //     template<typename... Args>
	// boost::any JacobianTransposeAction(unsigned int const wrtIn, unsigned int const wrtOut, boost::any const& vec, Args... args) {
	// // make sure the input and output number are valid
	// assert(numInputs<0 || wrtIn<numInputs);
	// assert(numOutputs<0 || wrtOut<numOutputs);
  //
	// // clear the outputs and derivative information
	// ClearDerivatives();
  //
	// // create the reference input vector
	// ref_vector<boost::any> inputs;
	// inputs.reserve(numInputs<0? 0 : numInputs);
  //
	// // begin calling the JacobianTransposeActionRecursive with the first input
	// return JacobianTransposeActionRecursive(wrtIn, wrtOut, vec, inputs, args...);
  //     }
  //
  //     /// Compute the action of the Jacobian transpose using finite differences
  //     /**
	//  Assume the input and output type are Eigen::VectorXd.
	//  @param[in] wrtIn The input number we are taking the Jacobian with respect to
	//  @param[in] wrtOut The output number we are taking the Jacobian with respect to
	//  @param[in] vec We want to apply the Jacobian transpose to this vector
	//  @param[in] args The inputs (may be more than one)
  //      */
  //     void JacobianTransposeActionByFD(unsigned int const wrtIn, unsigned int const wrtOut, boost::any const& vec, ref_vector<boost::any> const& inputs);

      /// Get the (unique) name of this work piece
      /**
	 \return The name of the muq::Modeling::WorkPiece
       */
      std::string const& Name();

      /// Set the name of this work piece
      void SetName(std::string const& newName);

      /// Get the input type (if we know it) for a specific input
      /**
	 The return input type name is "demangled" so it is more human readable.
	 @param[in] inputNum The input we want the name of
	 @param[in] demangle true (default): demangle the input so it is human-readable, false: do not demangle the input
	 \return If we know the input types, the input type name is returned.  If we do not know the input types, return ""
       */
      std::string InputType(unsigned int inputNum, bool const demangle = true) const;

      /// Get the length of a vector valued input with fixed size
      int InputSize(unsigned int inputNum) const;

      void SetInputSize(unsigned int inputNum, int newSize);

      /// Get the output type (if we know it) for a specific output
      /**
	 The return output type name is "demangled" so it is more human readable.
	 @param[in] outputNum The output we want the name of
	 @param[in] demangle true (default): demangle the input so it is human-readable, false: do not demangle the input
	 \return If we know the output types, the output type name is returned.  If we do not know the output types, return ""
       */
      std::string OutputType(unsigned int outputNum, bool const demangle = true) const;

      /// Get the output types
      std::map<unsigned int, std::string> OutputTypes() const;

      /// Get the input types
      std::map<unsigned int, std::string> InputTypes() const;

      /// Get the unique ID number
      /**
	 \return The ID number
       */
      unsigned int ID() const;

      /** @brief Get the average run time for one of the implemented methods.
       *   @details This function returns the average wall clock time (in milliseconds) for the EvaluateImpl, GradientImpl,
       * JacobianImpl, JacobianActionImpl, or HessianImp functions.
       *            If the function was never called, -1 is returned.
       *   @param[in] method The implemented function of interest.  Possible options are "Evaluate", "Gradient", "Jacobian",
       * "JacobianAction", or "Hessian"
       *   @return The average wall clock time in milli-seconds.
       */
      virtual double GetRunTime(const std::string& method="Evaluate") const;

      /** @brief get the number of times one of the implemented methods has been called.
       *   @details This function returns the number of times the EvaluateImpl, GradientImpl, JacobianImpl,
       * JacobianActionImpl, or HessianImp functions have been called.
       *   @param[in] method The implemented function of interest.  Possible options are "Evaluate", "Gradient", "Jacobian",
       * "JacobianAction", or "Hessian"
       *   @return An integer with the number of calls.
       */
      virtual unsigned long int GetNumCalls(const std::string& method = "Evaluate") const;

      /** @brief Resets the number of call and times.
       *   @details Sets the number of calls to each implemented function to zero and the recorded wall clock times to zero.
       */
      virtual void ResetCallTime();


      /// The number of inputs
      int numInputs;

      /// The number of outputs
      int numOutputs;

      /// Create vector of references from a vector of boost::any's
      static ref_vector<const boost::any> ToRefVector(std::vector<boost::any> const& anyVec);
      static ref_vector<const Eigen::VectorXd> ToRefVector(std::vector<Eigen::VectorXd> const& anyVec);

    protected:

      /// Check the input type
      /**
	 @param[in] inputNum The input number --- we are check that the type has the same type as this input
	 @param[in] type The type of the input
	 \return true: The input matches the specified type or no input type is specified, false: the input type does not match the specified type
       */
      bool CheckInputType(unsigned int const inputNum, std::string const& type) const;

      /// Check the output type
      /**
	 @param[in] outputNum The output number --- we are check that the computed type has the same type as this output
	 @param[in] type The type of the output
	 \return true: The output matches the specified type or no output type is specified, false: the output type does not match the specified type
       */
      bool CheckOutputType(unsigned int const outputNum, std::string const& type) const;


      /// Get the types from a vector of boost::any's
      /**
	 @param[in] vec A vector of boost::any's
	 \return A vector of the types of the boost::any's
       */
      std::vector<std::string> Types(std::vector<boost::any> const& vec) const;

      /// Clear outputs every time Evaluate is called
      /**
	 If true, muq::Modeling::WorkPiece::outputs is cleared everytime muq::Modeling::WorkPiece::Evaluate is called.  If false, the outputs are not cleared.  Defaults to true.
       */
      bool clearOutputs = true;

      /// The outputs
      /**
	 The outputs of this muq::Modeling::WorkPiece are filled by WorkPiece::EvaluateImpl().  If the number of outputs is specified (i.e., WorkPiece::numOutputs is not -1) then WorkPiece::Evaluate() checks to make sure the size of this vector is equal to WorkPiece::numOutputs after calling WorkPiece::EvaluateImpl().  If the output types are specified (i.e., WorkPiece::outputTypes is not an empty vector) then WorkPiece::Evaluate() checks that the output types match WorkPiece::outputTypes after calling WorkPiece::EvaluateImpl().
      */
      std::vector<boost::any> outputs = std::vector<boost::any>(0);

      /// The input types
      /**
	 Each element specifies the type of the corresponding input.  This vector must have the same number of elements as WorkPiece::numInputs or it is empty (default), which indicates that the input types are variable.
      */
      std::map<unsigned int, std::string> inputTypes;

      /// The output types
      /**
	 Each element specifies the type of the corresponding output.  This vector must have the same number of elements as WorkPiece::numOutputs or it is empty (default), which indicates that the output types are variable.
      */
      std::map<unsigned int, std::string> outputTypes;

      std::map<unsigned int, int> inputSizes;

      /// Convert a vector of input types to a map
      /**
	 The key in the map is the index of the vector.
	 @param[in] typesVec A vector of input types
	 \return A map of input types
       */
      std::map<unsigned int, std::string> Types(std::vector<std::string> const& typesVec) const;


      // The following variables keep track of how many times the Implemented functions, i.e. EvaluateImpl, GradientImpl,
      // etc... are called.
      unsigned long int numEvalCalls   = 0;

      // these variables keep track of the total wall-clock time spent in each of the Implemented functions.  They are in
      // units of milliseconds
      double evalTime   = 0;

      /// A unique ID number assigned by the constructor
      const unsigned int id;

      /// A unique name for this WorkPiece.  Defaults to <ClassName>_<id>
      std::string name;

    private:

      virtual std::string CreateName() const;


      /// User-implemented function that determines the behavior of this muq::Modeling::WorkPiece
      /**
	 This function determines how the WorkPiece::inputs determine WorkPiece::outputs.  Must be implemented by a child.

	 WorkPiece::Evaluate() calls this function after checking the inputs and storing them in WorkPiece::inputs.  This function populates WorkPiece::outputs, the outputs of this muq::Modeling::WorkPiece.  WorkPiece::Evaluate() checks the outputs after calling this function.
      */
      virtual void EvaluateImpl(ref_vector<boost::any> const& inputs) = 0;

      /// Creates WorkPiece::inputs when the WorkPiece::Evaluate is called with multiple arguments
      /**
	 @param[in] inputs The growing vector of inputs
	 @param[in] in The input corresponding to the \f$i^{th}\f$ input
	 @param[in] args The remaining (greater than \f$i\f$) inputs
	 \return The outputs of this muq::Modeling::WorkPiece
      */
      template<typename ith, typename... Args>
	std::vector<boost::any> const& EvaluateRecursive(ref_vector<boost::any> &inputs, ith const& in, Args... args) {
	const int inputNum = inputs.size();

	// we have not yet put all of the inputs into the map, the ith should be less than the total number
	assert(numInputs<0 || inputNum+1<numInputs);

	// check the input type
	assert(CheckInputType(inputNum, typeid(in).name()));

	// add the last input to the input vector
	const boost::any in_any(in);
	inputs.push_back(std::cref(in_any));

	// call with EvaluateRecursive with the remaining inputs
	return EvaluateRecursive(inputs, args...);
      }

      /// Creates WorkPiece::inputs when the WorkPiece::Evaluate is called with multiple arguments
      /**
	 @param[in] inputs The growing vector of inputs
	 @param[in] in The input corresponding to the last input
	 \return The outputs of this muq::Modeling::WorkPiece
      */
      template<typename last>
	std::vector<boost::any> const& EvaluateRecursive(ref_vector<boost::any> &inputs, last const& in) {

	const int inputNum = inputs.size();

	// this is the last input, the last one should equal the total number of inputs
	assert(numInputs<0 || inputNum+1==numInputs);

	// check the input type
	assert(CheckInputType(inputNum, typeid(in).name()));

	// add the last input to the input vector
	const boost::any in_any(in);
	inputs.push_back(std::cref(in_any));

	return Evaluate(inputs);
      }

  //     /// User-implemented function that to compute the Jacobian
  //     /**
	//  If th user does not implement this function, the Jacobian cannot be computed.  However, if both the input and the output type are Eigen::VectorXd, we default to finite difference.
	//  @param[in] wrtIn The input number we are taking the Jacobian with respect to
	//  @param[in] wrtOut The output number we are taking the Jacobian with respect to
	//  @param[in] inputs The vector of references to the inputs
  //     */
  //     virtual void JacobianImpl(unsigned int const wrtIn, unsigned int const wrtOut, ref_vector<boost::any> const& inputs);
  //
  //     /// Creates WorkPiece::inputs when the WorkPiece::Jacobian is called with multiple arguments
  //     /**
	//  @param[in] wrtIn The input number we are taking the Jacobian with respect to
	//  @param[in] wrtOut The output number we are taking the Jacobian with respect to
	//  @param[in] inputs The growing vector of inputs
	//  @param[in] in The input corresponding to the \f$i^{th}\f$ input
	//  @param[in] args The remaining (greater than \f$i\f$) inputs
	//  \return The Jacobian of this muq::Modeling::WorkPiece
  //     */
  //     template<typename ith, typename... Args>
	// boost::any JacobianRecursive(unsigned int const wrtIn, unsigned int const wrtOut, ref_vector<boost::any> &inputs, ith const& in, Args... args) {
	// const int inputNum = inputs.size();
  //
	// // we have not yet put all of the inputs into the map, the ith should be less than the total number
	// assert(numInputs<0 || inputNum+1<numInputs);
  //
	// // check the input type
	// assert(CheckInputType(inputNum, typeid(in).name()));
  //
	// // add the last input to the input vector
	// const boost::any in_any(in);
	// inputs.push_back(std::cref(in_any));
  //
	// // call with JacobianRecursive with the remaining inputs
	// return JacobianRecursive(wrtIn, wrtOut, inputs, args...);
  //     }
  //
  //     /// Creates WorkPiece::inputs when the WorkPiece::Jacobian is called with multiple arguments
  //     /**
	//  @param[in] wrtIn The input number we are taking the Jacobian with respect to
	//  @param[in] wrtOut The output number we are taking the Jacobian with respect to
	//  @param[in] inputs The growing vector of inputs
	//  @param[in] in The input corresponding to the last input
	//  \return The Jacobian of this muq::Modeling::WorkPiece
  //     */
  //     template<typename last>
	// boost::any JacobianRecursive(unsigned int const wrtIn, unsigned int const wrtOut, ref_vector<boost::any> &inputs, last const& in) {
	// const int inputNum = inputs.size();
  //
	// // this is the last input, the last one should equal the total number of inputs
	// assert(numInputs<0 || inputNum+1==numInputs);
  //
	// // check the input type
	// assert(CheckInputType(inputNum, typeid(in).name()));
  //
	// // add the last input to the input vector
	// const boost::any in_any(in);
	// inputs.push_back(std::cref(in_any));
  //
	// return Jacobian(wrtIn, wrtOut, inputs);
  //     }
  //
  //     /// User-implemented function that to compute the action of the Jacobian
  //     /**
	//  @param[in] wrtIn The input number we are taking the Jacobian with respect to
	//  @param[in] wrtOut The output number we are taking the Jacobian with respect to
	//  @param[in] vec We want to apply the Jacobian to this vector
	//  @param[in] inputs The vector of references to the inputs
  //     */
  //     virtual void JacobianActionImpl(unsigned int const wrtIn, unsigned int const wrtOut, boost::any const& vec, ref_vector<boost::any> const& inputs);
  //
  //     /// Creates WorkPiece::inputs when the WorkPiece::JacobianAction is called with multiple arguments
  //     /**
	//  @param[in] wrtIn The input number we are taking the Jacobian with respect to
	//  @param[in] wrtOut The output number we are taking the Jacobian with respect to
	//  @param[in] vec We want to apply the Jacobian to this vector
	//  @param[in] inputs The growing vector of inputs
	//  @param[in] in The input corresponding to the \f$i^{th}\f$ input
	//  @param[in] args The remaining (greater than \f$i\f$) inputs
	//  \return The action of the Jacobian of this muq::Modeling::WorkPiece on vec
  //     */
  //     template<typename ith, typename... Args>
	// boost::any JacobianActionRecursive(unsigned int const wrtIn, unsigned int const wrtOut, boost::any const& vec, ref_vector<boost::any>& inputs, ith const& in, Args... args) {
	// const int inputNum = inputs.size();
  //
	// // we have not yet put all of the inputs into the map, the ith should be less than the total number
	// assert(numInputs<0 || inputNum+1<numInputs);
  //
	// // check the input type
	// assert(CheckInputType(inputNum, typeid(in).name()));
  //
	// // add the last input to the input vector
	// const boost::any in_any(in);
	// inputs.push_back(std::cref(in_any));
  //
	// // call with JacobianActionRecursive with the remaining inputs
	// return JacobianActionRecursive(wrtIn, wrtOut, vec, inputs, args...);
  //     }
  //
  //     /// Creates WorkPiece::inputs when the WorkPiece::JacobianAction is called with multiple arguments
  //     /**
	//  @param[in] wrtIn The input number we are taking the Jacobian with respect to
	//  @param[in] wrtOut The output number we are taking the Jacobian with respect to
	//  @param[in] vec We want to apply the Jacobian to this vector
	//  @param[in] inputs The growing vector of inputs
	//  @param[in] in The input corresponding to the last input
	//  \return The action of the Jacobian of this muq::Modeling::WorkPiece on vec
  //     */
  //     template<typename last>
	// boost::any JacobianActionRecursive(unsigned int const wrtIn, unsigned int const wrtOut, boost::any const& vec, ref_vector<boost::any>& inputs, last const& in) {
	// const int inputNum = inputs.size();
  //
	// // this is the last input, the last one should equal the total number of inputs
	// assert(numInputs<0 || inputNum+1==numInputs);
  //
	// // check the input type
	// assert(CheckInputType(inputNum, typeid(in).name()));
  //
	// // add the last input to the input vector
	// const boost::any in_any(in);
	// inputs.push_back(std::cref(in_any));
  //
	// return JacobianAction(wrtIn, wrtOut, vec, inputs);
  //     }
  //
  //     /// User-implemented function that to compute the action of the Jacobian transpose
  //     /**
	//  @param[in] wrtIn The input number we are taking the Jacobian with respect to
	//  @param[in] wrtOut The output number we are taking the Jacobian with respect to
	//  @param[in] vec We want to apply the Jacobian transpose to this vector
	//  @param[in] inputs The vector of references to the inputs
  //     */
  //     virtual void JacobianTransposeActionImpl(unsigned int const wrtIn, unsigned int const wrtOut, boost::any const& vec, ref_vector<boost::any> const& inputs);
  //
  //     /// Creates WorkPiece::inputs when the WorkPiece::JacobianTransposeAction is called with multiple arguments
  //     /**
	//  @param[in] wrtIn The input number we are taking the Jacobian with respect to
	//  @param[in] wrtOut The output number we are taking the Jacobian with respect to
	//  @param[in] vec We want to apply the Jacobian transpose on this vector
	//  @param[in] inputs The growing vector of inputs
	//  @param[in] in The input corresponding to the \f$i^{th}\f$ input
	//  @param[in] args The remaining (greater than \f$i\f$) inputs
	//  \return The action of the Jacobian of this muq::Modeling::WorkPiece on vec
  //     */
  //     template<typename ith, typename... Args>
	// boost::any JacobianTransposeActionRecursive(unsigned int const wrtIn, unsigned int const wrtOut, boost::any const& vec, ref_vector<boost::any>& inputs, ith const& in, Args... args) {
	// const int inputNum = inputs.size();
  //
	// // we have not yet put all of the inputs into the map, the ith should be less than the total number
	// assert(numInputs<0 || inputNum+1<numInputs);
  //
	// // check the input type
	// assert(CheckInputType(inputNum, typeid(in).name()));
  //
	// // add the last input to the input vector
	// const boost::any in_any(in);
	// inputs.push_back(std::cref(in_any));
  //
	// // call with JacobianTransposeActionRecursive with the remaining inputs
	// return JacobianTransposeActionRecursive(wrtIn, wrtOut, vec, inputs, args...);
  //     }
  //
  //     /// Creates WorkPiece::inputs when the WorkPiece::JacobianTransposeAction is called with multiple arguments
  //     /**
	//  @param[in] wrtIn The input number we are taking the Jacobian with respect to
	//  @param[in] wrtOut The output number we are taking the Jacobian with respect to
	//  @param[in] vec We want to apply the Jacobian transpose to this vector
	//  @param[in] inputs The growing vector of inputs
	//  @param[in] in The input corresponding to the last input
	//  \return The action of the Jacobian of this muq::Modeling::WorkPiece on vec
  //     */
  //     template<typename last>
	// boost::any JacobianTransposeActionRecursive(unsigned int const wrtIn, unsigned int const wrtOut, boost::any const& vec, ref_vector<boost::any>& inputs, last const& in) {
	// const int inputNum = inputs.size();
  //
	// // this is the last input, the last one should equal the total number of inputs
	// assert(numInputs<0 || inputNum+1==numInputs);
  //
	// // check the input type
	// assert(CheckInputType(inputNum, typeid(in).name()));
  //
	// // add the last input to the input vector
	// const boost::any in_any(in);
	// inputs.push_back(std::cref(in_any));
  //
	// return JacobianTransposeAction(wrtIn, wrtOut, vec, inputs);
  //     }

      /// Clear muq::Modeling::WorkPiece::outputs when muq::Modeling::Evaluate is called
      void Clear();

      /// Destroy a boost any
      /**
	 Destroys the object associated with a boost::any.  If the boost::any is not a smart pointer, for example, the function must destory it.

	 This function knows how to destory commonly used objects (those in muq::Modeling::WorkPiece::types).  Overloading muq::Modeling::WorkPiece::DestroyAnyImpl() destroys other types of objects
       */
      void DestroyAny(boost::any& obj) const;

      /// Destroy a boost any
      /**
	 By default, this function does nothing.  It can be overloaded by the user to destroy objects contained within boost::any's
       */
      virtual void DestroyAnyImpl(boost::any& obj) const;

      /// Creates a unique ID number, must be called by the constructor
      static unsigned int CreateID();

      /// A map of common types
      /**
	 A list of common types and their corresponding names:
	 <ol>
	 <li> "N_Vector" for <TT>N_Vector</TT> type (requires Sundials)
	 <li> "N_Vector vector" for <TT>std::vector<N_Vector></TT> type (requires Sundials)
	 <li> "DlsMat" for <TT>DlsMat</TT> type (requires Sundials)
	 <li> "DlsMat vector" for <TT>std::vector<DlsMat></TT> type (requires Sundials)
	 </ol>
       */

   //const std::map<std::string, std::string> types;
    }; // class WorkPiece
  } // namespace Modeling
} // namespace muq

#endif
