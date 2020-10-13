#ifndef MATERNKERNEL_H
#define MATERNKERNEL_H

#include "MUQ/Approximation/GaussianProcesses/KernelImpl.h"

#include <cmath>
#include <stdexcept>

#include <boost/property_tree/ptree_fwd.hpp>

#include <boost/math/constants/constants.hpp>



namespace muq
{
namespace Approximation
{


/**

@class MaternKernel
@ingroup CovarianceKernels
This class implements a kernel of the form
\f[
k(x,y) = \sigma^2 \frac{2^{1-\nu}}{\Gamma(\nu)}\left(\frac{\sqrt{2\nu}}{l}\tau\right)^\nu K_\nu\left(\frac{\sqrt{2\nu}}{l}\tau\right),
\f]
where \f$\Gamma(\cdot)\f$ is the gamma function, \f$K_\nu(\cdot)\f$ is a modified Bessel function of the second kind, \f$\nu\f$ is a smoothness parameter, \f$l\f$ is a lengthscale, and \f$\sigma^2\f$ is the variance.  Note that we only allow values of \f$\nu\f$ that take the form \f$i-0.5\f$ for some positive integer \f$i\f$.  Typical choices are \f$\nu=1/2\f$, \f$\nu=3/2\f$, and \f$\nu=5/2\f$.   In the limit, \f$\nu\rightarrow\infty\f$, this Matern kernel converges to the squared exponential kernel (muq::Approximation::SquaredExpKernel) and for \f$\nu=1/2\f$, the Matern kernel is equivalent to the exponential kernel associated with the Ornstein-Uhlenbeck process.

Note: the smoothness parameter \f$\nu\f$ is not optimized as a hyperparameter.
 */
class MaternKernel : public KernelImpl<MaternKernel>
{

public:

  MaternKernel(unsigned              dimIn,
               std::vector<unsigned> dimInds,
               double                sigma2In,
               double                lengthIn,
               double                nuIn) : MaternKernel(dimIn,
                                                          dimInds,
                                                          sigma2In,
                                                          lengthIn,
                                                          nuIn,
                                                          {0.0, std::numeric_limits<double>::infinity()},
                                                          {1e-10, std::numeric_limits<double>::infinity()})
  {};

    MaternKernel(unsigned              dimIn,
                 std::vector<unsigned> dimInds,
                 double                sigma2In,
                 double                lengthIn,
                 double                nuIn,
                 Eigen::Vector2d       sigmaBounds,
                 Eigen::Vector2d       lengthBounds);


    MaternKernel(unsigned        dimIn,
                 double          sigma2In,
                 double          lengthIn,
                 double          nuIn) : MaternKernel(dimIn,
                                                      sigma2In,
                                                      lengthIn,
                                                      nuIn,
                                                      {0.0, std::numeric_limits<double>::infinity()},
                                                      {1e-10, std::numeric_limits<double>::infinity()})
    {};

    MaternKernel(unsigned        dimIn,
                 double          sigma2In,
                 double          lengthIn,
                 double          nuIn,
                 Eigen::Vector2d sigmaBounds,
                 Eigen::Vector2d lengthBounds);


    virtual ~MaternKernel(){};


    template<typename ScalarType1, typename ScalarType2, typename ScalarType3>
    void FillBlockImpl(Eigen::Ref<const Eigen::Matrix<ScalarType1, Eigen::Dynamic, 1>> const& x1,
                       Eigen::Ref<const Eigen::Matrix<ScalarType1, Eigen::Dynamic, 1>> const& x2,
                       Eigen::Ref<const Eigen::Matrix<ScalarType2, Eigen::Dynamic, 1>> const& params,
                       Eigen::Ref<Eigen::Matrix<ScalarType3,Eigen::Dynamic, Eigen::Dynamic>>  block) const
    {
      int p = round(nu-0.5);

      ScalarType1 dist = (x1-x2).norm();

      block(0,0) = 0.0;
      for(int i=0; i<=p; ++i)
        block(0,0) +=  weights(i) * pow(sqrt(8.0*nu)*dist/params(1), p-i);

      block(0,0) *= params(0)*exp(-sqrt(2.0*nu)*dist / params(1)) * scale;
    }

    virtual std::tuple<std::shared_ptr<muq::Modeling::LinearSDE>, std::shared_ptr<muq::Modeling::LinearOperator>, Eigen::MatrixXd> GetStateSpace(boost::property_tree::ptree sdeOptions = boost::property_tree::ptree()) const override;

private:

    const double nu;
    const double scale; // std::pow(2.0, 1.0-nu)/boost::math::tgamma(nu)

    const Eigen::VectorXd weights;

    void CheckNu() const;

    static Eigen::VectorXd BuildWeights(int p);

    static inline int Factorial(int n){
      return (n == 1 || n == 0) ? 1 : Factorial(n - 1) * n;
    }
};

}
}


#endif
