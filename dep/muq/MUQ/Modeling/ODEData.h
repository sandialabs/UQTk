#ifndef ODEDATA_H_
#define ODEDATA_H_

#include "MUQ/config.h"

#if MUQ_HAS_PARCER==1
#include <parcer/Communicator.h>
#endif

#include "MUQ/Modeling/ModPiece.h"

namespace muq {
  namespace Modeling {
    /// A helper class for Sundial's time integration
    /**
       Stores necessary muq::Modeling::WorkPiece's and thier inputs so Sundials can hold that information.
     */
    class ODEData {
    public:

      /// Construct basic ode data
      /**
	 @param[in] rhs A muq::Modeling::WorkPiece that evalautes the right hand side of an ODE
	 @param[in] inputs The inputs to the rhs --- the first is the state, the rest are constant in time
	 @param[in] autonomous Is the RHS autonomous?
	 @param[in] wrtIn The input we are computing the derivative wrt --- negative indicates no derivative is being computed
       */
#if MUQ_HAS_PARCER==1
      ODEData(std::shared_ptr<ModPiece> const& rhs, ref_vector<Eigen::VectorXd> const& refinputs, bool const autonomous, int const wrtIn, std::shared_ptr<parcer::Communicator> const& comm);
#else
      ODEData(std::shared_ptr<ModPiece> const& rhs, ref_vector<Eigen::VectorXd> const& refinputs, bool const autonomous, int const wrtIn);
#endif

      /// Construct with root function
      /**
	 @param[in] rhs A muq::Modeling::WorkPiece that evalautes the right hand side of an ODE
	 @param[in] root A muq::Modeling::WorkPiece that evalautes a function we are trying to find the root of along an orbit of the ODE
	 @param[in] inputs The inputs to the rhs --- the first is the state, the rest are constant in time
	 @param[in] autonomous Is the RHS autonomous?
	 @param[in] wrtIn The input we are computing the derivative wrt --- negative indicates no derivative is being computed
       */
      ODEData(std::shared_ptr<ModPiece> const& rhs, std::shared_ptr<ModPiece> const& root, ref_vector<Eigen::VectorXd> const& refinputs, bool const autonomous, int const wrtIn);

      virtual ~ODEData() = default;

      /// Update the time and state inputs
      void UpdateInputs(Eigen::VectorXd const& state, double const time);

      /// The right hand side of the ODE
      std::shared_ptr<ModPiece> rhs;

      /// A function we are trying to find the root of along an orbit of the ODE (nullptr if we are not doing a root finding problem)
      std::shared_ptr<ModPiece> root;

      /// The inputs to the rhs --- the first is the state, the rest are constant in time
      std::vector<Eigen::VectorXd> inputs;

      /// Is the RHS autonomous?
      const bool autonomous;

      /// The input we are computing the derivative wrt --- negative indicates no derivative is being computed
      const int wrtIn = -1;

#if MUQ_HAS_PARCER==1
      std::shared_ptr<parcer::Communicator> comm = nullptr;
#endif

    private:
    };
  } // namespace Modeling
} // namespace muq

#endif
