#ifndef ODEBASE_H_
#define ODEBASE_H_

#include <boost/property_tree/ptree.hpp>

#include "MUQ/Modeling/ModPiece.h"
#include "MUQ/Modeling/ODEData.h"

namespace muq {
  namespace Modeling {

    /// A bass class to integrate ODE's
    class ODEBase : public ModPiece {
    public:

      /**
	    @param[in] rhs The right hand side of the ODE
      @param[in] inputSizes The input input sizes
      @param[in] outputSizes The output sizes
	    @param[in] pt A boost::property_tree::ptree with options/tolerances for the ODE integrator
      */
#if MUQ_HAS_PARCER==1
  ODEBase(std::shared_ptr<ModPiece> const& rhs, Eigen::VectorXi const& inputSizes, Eigen::VectorXi const& outputSizes, boost::property_tree::ptree const& pt, std::shared_ptr<parcer::Communicator> const& comm = nullptr);
#else
      ODEBase(std::shared_ptr<ModPiece> const& rhs, Eigen::VectorXi const& inputSizes, Eigen::VectorXi const& outputSizes, boost::property_tree::ptree const& pt);
#endif

      virtual ~ODEBase();

    protected:

      /// Are we computing the Jacobian, the action of the Jacobian, or the action of the Jacobian transpose
      enum DerivativeMode {
	/// The Jacobian
	Jac,
	/// The action of the Jacobian
	JacAction,
	/// The action of the Jacobian transpose
	JacTransAction
      };

      /// Check the return flag of a Sundials function
      /**
	 @param[in] flagvalue The value of the Sundials flag
	 @param[in] funcname The name of the Sundials function
	 @param[in] opt An option to determine how to check the flag, 0: check if flag is nullptr, 1: flag is an int, check if flag<0 (indicates Sundials error)
	 \return false: failure, true: success
       */
      bool CheckFlag(void* flagvalue, std::string const& funcname, unsigned int const opt) const;

      /// Alloc memory and set up options for the Sundials solver
      /**
	 @param[in] cvode_mem The Sundials solver
	 @param[in] state The initial state
	 @param[in] data An object that holds the RHS inputs and can evaluate the RHS
       */
      void CreateSolverMemory(void* cvode_mem, N_Vector const& state, std::shared_ptr<ODEData> data) const;

      int CreateSolverMemoryB(void* cvode_mem, double const timeFinal, N_Vector const& lambda, N_Vector const& nvGrad, std::shared_ptr<ODEData> data) const;

      /// Deal with Sundials errors
      /**
	 Sundials will call this function if it runs into a problem
	 @param[in] error_code Sundials error code
	 @param[in] module The name of the CVODES module reporting the error
	 @param[in] function The name of the function in which the error occured
	 @param[in] msg The error message
	 @param[in] user_data A pointer to an muq::Modeling::ODEData
       */
      static void ErrorHandler(int error_code, const char *module, const char *function, char *msg, void *user_data);

      /// Evaluate the right hand side
      /**
	 @param[in] time The current time
	 @param[in] state The current state
	 @param[out] deriv The derivative of the state with respect to time
	 @param[in] user_data A pointer to an muq::Modeling::ODEData
       */
      static int EvaluateRHS(realtype time, N_Vector state, N_Vector deriv, void *user_data);

      static int AdjointRHS(realtype time, N_Vector state, N_Vector lambda, N_Vector deriv, void *user_data);

      /// Evaluate the Jacobian of the right hand side
      /**
	 @param[in] N
	 @param[in] time The current time
	 @param[in] state The current state
	 @param[in] rhs The derivative of the state with respect to time
	 @param[out] jac The Jacobian of the right hand side with respect to the current state
	 @param[in] user_data A pointer to an muq::Modeling::ODEData
	 @param[in] tmp1
	 @param[in] tmp2
	 @param[in] tmp3
       */
      static int RHSJacobian(long int N, realtype time, N_Vector state, N_Vector rhs, DlsMat jac, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

      static int AdjointJacobian(long int N, realtype time, N_Vector state, N_Vector lambda, N_Vector rhs, DlsMat jac, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

      /// Evaluate the action of the Jacobian of the right hand side
      /**
	 @param[in] v The vector the Jacobian is acting on
	 @param[out] Jv The action of the Jacobian on v
	 @param[in] time The current time
	 @param[in] state The current state
	 @param[in] rhs The derivative of the state with respect to time
	 @param[in] user_data A pointer to an muq::Modeling::ODEData
	 @param[in] tmp
       */
      static int RHSJacobianAction(N_Vector v, N_Vector Jv, realtype time, N_Vector state, N_Vector rhs, void *user_data, N_Vector tmp);

      static int AdjointJacobianAction(N_Vector target, N_Vector output, realtype time, N_Vector state, N_Vector lambda, N_Vector adjRhs, void *user_data, N_Vector tmp);

      /// Sundials uses this function to compute the derivative of the state at each timestep
      /**
	 @param[in] Ns The number of sensitivities
	 @param[in] time The current time
	 @param[in] y The current state
	 @param[in] ydot The derivative of the current state with respect to time
	 @param[in] ys
	 @param[in] ySdot The sensitivties
	 @param[in] user_data A pointer to an muq::Modeling::ODEData
       */
      static int ForwardSensitivityRHS(int Ns, realtype time, N_Vector y, N_Vector ydot, N_Vector *ys, N_Vector *ySdot, void *user_data, N_Vector tmp1, N_Vector tmp2);

      static int AdjointQuad(realtype time, N_Vector state, N_Vector lambda, N_Vector quadRhs, void *user_data);

      /// Set up the solver for sensitivity information
      /**
	 @param[in] cvode_mem The Sundials solver
	 @param[in] paramSize The size of the input parameter we are differenating wrt
	 @param[in,out] sensState This will become the 'current' Jacobian
       */
      void SetUpSensitivity(void *cvode_mem, unsigned int const paramSize, N_Vector *sensState) const;

      /// Which linear solver should we use?
      enum LinearSolver {
	/// Dense solver
	Dense,
	/// SPGMR
	SPGMR,
	/// SPBCG
	SPBCG,
	/// SPTFQMR
	SPTFQMR
      };

      /// The right hand side of the ODE
      std::shared_ptr<ModPiece> rhs;

      /// Linear solver method
      LinearSolver slvr;

      /// The relative tolerance
      const double reltol;

      /// The absolute tolerance
      const double abstol;

      /// The maximum time step size
      const double maxStepSize;

      /// The maximum number of time steps
      const unsigned int maxNumSteps;

      /// Multistep method
      int multiStep;

      /// Nonlinear solver method
      int solveMethod;

      /// Is the RHS autonomous?
      const bool autonomous;

      /// Check point gap
      const unsigned int checkPtGap;

#if MUQ_HAS_PARCER==1
      /// The global size of the state vector
      const unsigned int globalSize = std::numeric_limits<unsigned int>::quiet_NaN();

      std::shared_ptr<parcer::Communicator> comm = nullptr;
#endif

    private:

    };

  } // namespace Modeling
} // namespace muq

#endif
