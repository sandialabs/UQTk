#ifndef INVERSEGAMMA_H_
#define INVERSEGAMMA_H_

#include "MUQ/Modeling/Distributions/Distribution.h"

namespace muq {
  namespace Modeling {

    /** Defines the inverse gamma distribution, which has probability density function
    \f[
    \pi(x) = \frac{\beta^\alpha}{\Gamma(\alpha)}x^{-\alpha-1}\exp\left(-\frac{\beta}{x}\right),
    \f]
    where \f$\alpha\f$ and \f$\beta\f$ are parameters in the distribution.
    */
    class InverseGamma : public Distribution {
    public:

      InverseGamma(double       alphaIn,
                   double       betaIn);


      virtual ~InverseGamma() = default;


      const double alpha;
      const double beta;

    private:

      const double logConst; // log( beta^alpha / Gamma(alpha) )

      virtual double LogDensityImpl(ref_vector<Eigen::VectorXd> const& inputs) override;

      /// Sample the distribution
      virtual Eigen::VectorXd SampleImpl(ref_vector<Eigen::VectorXd> const& inputs) override;



    };
  } // namespace Modeling
} // namespace muq

#endif
