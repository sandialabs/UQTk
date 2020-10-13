#ifndef COSTFUNCTION_H_
#define COSTFUNCTION_H_

#include "MUQ/Utilities/VariadicMacros.h"

#include "MUQ/Modeling/ModPiece.h"

namespace muq {
  namespace Optimization {
    /// The cost function for an optimization routine
    /**
       The cost function has the form:
       \f{equation}{
       c = J(x; \theta_1, ..., \theta_1),
       \f}
       where, \f$c \in \mathbb{R}\f$.  The first input is always \f$x\f$ and the remaining inputs are \f$\theta_{1:n}\f$.
     */
    class CostFunction : public muq::Modeling::ModPiece {
    public:

      /**
	 @param[in] inputSizes The dimension of each input
       */
      CostFunction(Eigen::VectorXi const& inputSizes);

      virtual ~CostFunction();

      /// The value of the cost function
      /**
	 @param[in] input The inputs \f$x\f$, \f$\theta_{1:n}\f$
	 \return The value of the cost function
       */
      double Cost(muq::Modeling::ref_vector<Eigen::VectorXd> const& input);

      /// The value of the cost function
      /**
	 @param[in] args The inputs \f$x\f$, \f$\theta_{1:n}\f$
	 \return The value of the cost function
       */
      template<typename... Args>
	inline double Cost(Args const&... args) {
	return Evaluate(args...).at(0) (0);
      }

      /// The gradient of the cost function
      /**
	 @param[in] inputDimWrt Which input are we taking the derivative with respect to?
	 @param[in] input The inputs \f$x\f$, \f$\theta_{1:n}\f$
	 @param[in] sensitivity The sensitivity vector
	 \return The gradient of the cost function
       */
      Eigen::VectorXd const& Gradient(unsigned int const inputDimWrt,
                                      std::vector<Eigen::VectorXd> const& input,
                                      Eigen::VectorXd const& sensitivity);
      
      /// The gradient of the cost function
      /**
	 @param[in] inputDimWrt Which input are we taking the derivative with respect to?
	 @param[in] args The inputs \f$x\f$, \f$\theta_{1:n}\f$ and the sensitivity vector
	 \return The gradient of the cost function
       */
      template<typename... Args>
      inline Eigen::VectorXd const& Gradient(unsigned int const inputDimWrt,
                                             Args const&... args) {
	return ModPiece::Gradient(0, inputDimWrt, args...);
      }

      /// The Hessian of the cost function
      /**
         @param[in] inputDimWrt Which input are we taking the 2nd derivative with respect to?
         @param[in] input The inputs \f$x\f$, \f$\theta_{1:n}\f$
         \return The Hessian of the cost function
      */
      Eigen::MatrixXd Hessian(unsigned int const inputDimWrt,
                              std::vector<Eigen::VectorXd> const& input);


      /// The Hessian of the cost function using finite difference
      /**
         @param[in] inputDimWrt Which input are we taking the 2nd derivative with respect to?
         @param[in] input The inputs \f$x\f$, \f$\theta_{1:n}\f$
         \return The Hessian of the cost function
      */
      Eigen::MatrixXd HessianByFD(unsigned int const inputDimWrt,
                                  std::vector<Eigen::VectorXd> const& input);

      
      /// The Hessian of the cost function
      /**
         @param[in] inputDimWrt Which input are we taking the 2nd derivative with respect to?
         @param[in] input The inputs \f$x\f$, \f$\theta_{1:n}\f$
         @param[in] vec Vector to which the Hessian is applied
         \return The Hessian action on vec
      */
      Eigen::MatrixXd ApplyHessian(unsigned int const inputDimWrt,
                                   std::vector<Eigen::VectorXd> const& input,
                                   Eigen::VectorXd const& vec);

    private:

      /// The value of the cost function
      /**
	 Must be implemented by the user
	 @param[in] args The inputs \f$x\f$, \f$\theta_{1:n}\f$
	 \return The value of the cost function
       */
      virtual double CostImpl(muq::Modeling::ref_vector<Eigen::VectorXd> const& input) = 0;

      /// The value of the cost function
      /**
	 @param[in] args The inputs \f$x\f$, \f$\theta_{1:n}\f$
	 \return The value of the cost function
       */
      virtual void EvaluateImpl(muq::Modeling::ref_vector<Eigen::VectorXd> const& input) override;

      /// Compute the gradient of the cost function
      /**
	 @param[in] inputDimWrt Which input are we taking the derivative with respect to?
	 @param[in] args The inputs \f$x\f$, \f$\theta_{1:n}\f$ and the sensitivity vector
       */
      virtual void GradientImpl(unsigned int const outputDimWrt, unsigned int const inputDimWrt, muq::Modeling::ref_vector<Eigen::VectorXd> const& input, Eigen::VectorXd const& sensitivity) override;

      /// Compute the gradient of the cost function
      /**
	 Should be implemented by the user.
	 @param[in] inputDimWrt Which input are we taking the derivative with respect to?
	 @param[in] args The inputs \f$x\f$, \f$\theta_{1:n}\f$ and the sensitivity vector
       */
      virtual void GradientImpl(unsigned int const inputDimWrt, muq::Modeling::ref_vector<Eigen::VectorXd> const& input, Eigen::VectorXd const& sensitivity);
    };
  } // namespace Optimization
} // namespace muq

#endif
