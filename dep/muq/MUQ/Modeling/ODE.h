#ifndef ODE_H_
#define ODE_H_

#include "MUQ/Modeling/ODEBase.h"

namespace muq {
  namespace Modeling {
    class ODE : public ODEBase {
    public:
      /**
	 The first input is the initial state (at \f$t=0\f$).  It is also the first input to the right hand side muq::Modeling::ModPiece.

	 The next set of inputs are the inputs to the right hand side muq::Modeling::ModPiece.  If the right hand side input takes 2 inputs besides the state, these correspond to inputs 2 and 3 of the muq::Modeling::ODEBase.

	 The final input is the list of foutput times.
	 @param[in] rhs The right hand side of the ODE
	 @param[in] pt A boost::property_tree::ptree with options/tolerances for the ODE integrator
       */
#if MUQ_HAS_PARCER==1
      ODE(std::shared_ptr<ModPiece> const& rhs,  boost::property_tree::ptree pt, std::shared_ptr<parcer::Communicator> const& comm = nullptr);
#else
      ODE(std::shared_ptr<ModPiece> const& rhs,  boost::property_tree::ptree pt);
#endif

      virtual ~ODE();

    private:

      /**
      @param[in] rhs The right hand side of the ODE
   	  @param[in] pt A boost::property_tree::ptree with options/tolerances for the ODE integrator
      */
      static Eigen::VectorXi InputSizes(std::shared_ptr<ModPiece> const& rhs,  boost::property_tree::ptree pt);

      /**
      @param[in] rhs The right hand side of the ODE
      @param[in] pt A boost::property_tree::ptree with options/tolerances for the ODE integrator
      */
      static Eigen::VectorXi OutputSizes(std::shared_ptr<ModPiece> const& rhs,  boost::property_tree::ptree pt);

      /// Integrate the ODE forward in time
      /**
	 \f$M\f$ inputs:
	 <ul>
	 <li> The first \f$N\f$ are the inputs to muq::Modeling::ODE::rhs
	 <li> The second \fM-N\f$ are the times where we want to return the state (either vectors or scalars)
	 </ul>
	 @param[in] inputs The inputs (see description)
       */
      virtual void EvaluateImpl(ref_vector<Eigen::VectorXd> const& inputs) override;

      /// Evaluate the Jacobian of the state wrt each parameter
      /**
	 Returns the Jacobian of the state with respect to an input parameter at each timestep.  Returns a vector of Jacobians at each timestep if the output has more than one time.

	 \f$M\f$ inputs:
	 <ul>
	 <li> The first \f$N\f$ are the inputs to muq::Modeling::ODE::rhs
	 <li> The second \fM-N\f$ are the times where we want to return the state (either vectors or scalars)
	 </ul>
   @param[in] wrtOut We are computing the derivative of this output
	 @param[in] wrtIn We are computing the derivative with repsect to this parameter
	 @param[in] inputs The inputs (see description)
       */
      virtual void JacobianImpl(unsigned int const wrtOut, unsigned int const wrtIn, ref_vector<Eigen::VectorXd> const& inputs) override;

      /// Evaluate the action of the Jacobian transpose on a given vector wrt each parameter
      /**
   @param[in] wrtOut We are computing the derivative of this output
	 @param[in] wrtIn We are computing the derivative with repsect to this parameter
   @param[in] inputs The inputs
	 @param[in] vec The vector the Jacobian transpose is acting on at eact timestep
       */
      virtual void GradientImpl(unsigned int const wrtOut, unsigned int const wrtIn, ref_vector<Eigen::VectorXd> const& inputs, Eigen::VectorXd const& vec) override;

      /// Integrate forward in time, possibly keeping track fo forward sensitivities
      /**
	 @param[in] inputs The inputs to the right hand side and the output times
	 @param[in] wrtIn We are computing the derivative with repsect to this parameter
	 @param[in] vec The vector the Jacobian transpose is acting on at eact timestep
	 @param[in] mode Are we computing the Jacobian, Jacobian action, or Jacobian transpose action?
       */
      void Integrate(ref_vector<Eigen::VectorXd> const& inputs, int const wrtIn = -1,  N_Vector const& vec = nullptr, DerivativeMode const& mode = DerivativeMode::Jac);

      /// Integrate forward in time without computing derivative information
      /**
	 @param[in] cvode_mem The Sundials solver
	 @param[in,out] state The state --- begins at initial state, ends at final state
	 @param[in] outputTimes The output times the user has requested
       */
      void ForwardIntegration(void *cvode_mem, N_Vector& state, Eigen::VectorXd const& outputTimes);

      void BackwardIntegration(ref_vector<Eigen::VectorXd> const& inputs, unsigned int const wrtIn, Eigen::VectorXd const& vec);

      /// Integrate forward in time, computing derivative information
      /**
	 @param[in] cvode_mem The Sundials solver
	 @param[in,out] state The state --- begins at initial state, ends at final state
	 @param[in] wrtIn We are computing the derivative with repsect to this parameter
	 @param[in] outputTimes The output times corresponding to the output we are differentiating
	 @param[in] rhsInputs The inputs to the right hand side
       */
      void ForwardSensitivity(void *cvode_mem, N_Vector& state, unsigned int const wrtIn, Eigen::VectorXd const& outputTimes, ref_vector<Eigen::VectorXd> const& rhsInputs);

      static int AdjointQuad(realtype time, N_Vector state, N_Vector lambda, N_Vector quadRhs, void *user_data);
    };
  } // namespace Modeling
} // namespace muq

#endif
