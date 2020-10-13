#ifndef ROOTFINDINGIVP_H_
#define ROOTFINDINGIVP_H_

#include "MUQ/Modeling/ODEBase.h"

namespace muq {
  namespace Modeling {
    /// A rootfinding initial value problem --- find the root of a function along an orbit of an ODE
    class RootfindingIVP : public ODEBase {
    public:

      /**
	 The first input is the initial state (at \f$t=0\f$).  It is also the first input to the both right hand side and the root muq::Modeling::WorkPiece.  The type must be the same in both sub-models (if it is known).  This must a N_Vector type (the vectors that Sundials uses).

	 The next set of inputs are the inputs to the right hand side muq::Modeling::WorkPiece.  If the right hand side input takes 2 inputs besides the state, these correspond to inputs 2 and 3 of the root finder.   Their types are known if the types are known by the rhs muq::Modeling::WorkPiece.

	 The next set of inputs are the inputs to the root muq::Modeling::WorkPiece.  If the root input takes 2 inputs besides the state, these correspond to inputs 4 and 5 of the root finder.  Their types are known if the types are known by the rhs muq::Modeling::WorkPiece.

	 An optional final input is a vector of times to save the state.  If the user gives muq::Modeling::RootfinderIVP this (optional) input it integrates until it finds a root, saving the state at each of these times.

	 The first output is the state at the root.  This output has the same type as the first input to the right hand side and the root (if either is known).  Note that the "the root" is the first state such that any one of the outputs from the root function is zero.  It is boost::none if no root was found.

	 The second output is the time at which the root (first output) was obtained.

	 The third output exists if the user gives muq::Modeling::RootfinderIVP times to save the state (optional final input).  This output is a vector --- std::vector<StateType> --- of states at the specified times.

	 Parameters options (for boost::property_tree::ptree):
	 <ol>
	 <li> All the parameter listed in muq::Modeling::ODEBase
	 <li> The maximum number of steps (<EM>Rootfinder.MaxSteps</EM>)
	 <ul>
	 <li> Defaults to \f$10^{10}\f$
	 </ul>
	 <li> The maximum amount of time to integrate (<EM>Rootfinder.MaxTime</EM>)
	 <ul>
	 <li> Defaults to \f$10^{3}\f$
	 </ul>
	 <li> The maximum number of error tests (<EM>Rootfinder.MaxErrorTests</EM>)
	 <ul>
	 <li> Defaults to \f$100\f$
	 </ul>
	 </ol>
	 @param[in] rhs A muq::Modeling::WorkPiece that evaluates the right hand side of the ODE
	 @param[in] root A muq::Modeling::WorkPiece whose outputs are double's --- we integrate the ODE until we find the first root of one of these outputs
	 @param[in] pt A boost::property_tree::ptree with options/tolerances for the ODE integrator
	 @param[in] algebra A muq::Utilities::AnyAlgebra used to manipulate the state and input parameters (defaults to the MUQ default)
       */
      RootfindingIVP(std::shared_ptr<WorkPiece> rhs, std::shared_ptr<WorkPiece> root, boost::property_tree::ptree const& pt, std::shared_ptr<muq::Utilities::AnyAlgebra> algebra = std::make_shared<muq::Utilities::AnyAlgebra>());

      virtual ~RootfindingIVP();
      
    private:

      /// Integrate the ODE until we find a root
      /**
	 @param[in] inputs The inputs (first: state, next group: rhs inputs, next group: root inputs, final: eval times (optional))
      */
      virtual void EvaluateImpl(ref_vector<boost::any> const& inputs) override;

      /// Compute the Jacobian of the state at the root
      /**
	 @param[in] wrtIn We are computing the derivative with repsect to this parameter
	 @param[in] wrtOut We are computing the derivative of this output
	 @param[in] inputs The inputs (first: state, next group: rhs inputs, next group: root inputs, final: eval times (optional))
      */
      virtual void JacobianImpl(unsigned int const wrtIn, unsigned int const wrtOut, ref_vector<boost::any> const& inputs) override;

      /// Run the CVODES integrator
      /**
	 @param[in] inputs The inputs (first: state, next group: rhs inputs, next group: root inputs, final: eval times (optional))
	 @param[in] wrtIn We are computing the derivative with repsect to this parameter
	 @param[in] wrtOut We are computing the derivative of this output
	 @param[in] mode Are we computing the Jacobian, Jacobian action, or Jacobian transpose action?
	 \return A vector of length root->numOutputs, 0 indicates that output is not zero (no root found), 1 indicates that output is zero (a root was found)
       */
      Eigen::VectorXi FindRoot(ref_vector<boost::any> const& inputs, int const wrtIn = -1, int const wrtOut = -1, DerivativeMode const& mode = DerivativeMode::Jac);

      /// Update the input and output types based on the rhs and root muq::Modeling::WorkPiece's
      /**
	 Add the root input types to the inputs and set the second output type.  Note the (optional) third output (the vector of state's at specified times) is std::vector<StateType> but we can't set this.
       */
      void UpdateInputOutputTypes();

      /// Evaluate the root function
      /**
	 @param[in] t The current time
	 @param[in] state The state
	 @param[out] root The roots
	 @param[in] user_data A pointer to an muq::Modeling::ODEData
       */
      static int EvaluateRoot(realtype t, N_Vector state, realtype *root, void *user_data);

      /// The root function
      std::shared_ptr<WorkPiece> root;

      /// The maximum number of steps the timesteper can take
      const unsigned int maxSteps;

      /// The maximum amount of time to integrate
      const double maxTime;

      /// The maximum number of error test failures
      const unsigned int maxErrorTests;

    };
  } // namespace Modeling
} // namespace muq

#endif
