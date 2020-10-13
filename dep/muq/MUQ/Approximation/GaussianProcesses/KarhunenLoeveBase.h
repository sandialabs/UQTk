#ifndef KARHUNENLOEVEBASE_H
#define KARHUNENLOEVEBASE_H

#include <Eigen/Core>

namespace muq
{
namespace Approximation
{
    class KarhunenLoeveBase
    {

    public:

      virtual Eigen::MatrixXd GetModes(Eigen::Ref<const Eigen::MatrixXd> const& pts) const = 0;

      /** @brief Evaluate the KL expansion at a new pt given known coefficients.
          @details Recall that a truncated KL expansion of a random process \f$u(x,\omega)\f$ takes the form
      \f[
      u(x,\omega) = \sum_{k=1}^N \phi_k(x) z_k(\omega),
      \f]
      where \f$\phi_k(x)\f$ is the \f$k^{th}\f$ KL mode and each \f$z_k(\omega)\f$ is an independent standard normal random variable. This function evaluates the expansion given a point \f$x\f$ and a vector containing each \f$z_k\f$.   Thus, by fixing the coefficient vector \f$z_k\f$, it is possible to evaluate the sample of the Gaussian process at many different points.
      @param[in] pt The point \f$x\f$ that where you want to evaluate the KL expansion.
      @param[in] coeffs The coefficients in the expansion.  The coefficients should be drawn as iid standard normal random variables to generate a sample of the Gaussian process.
      */
      virtual Eigen::VectorXd Evaluate(Eigen::Ref<const Eigen::MatrixXd> const& pts,
                                       Eigen::Ref<const Eigen::VectorXd> const& coeffs) const
      {
        return GetModes(pts) * coeffs;
      };

      virtual unsigned int NumModes() const = 0;

    }; // class KarhunenLoeveBase
  }
}


#endif
