#ifndef WHITENOISEKERNEL_H
#define WHITENOISEKERNEL_H

#include "MUQ/Approximation/GaussianProcesses/KernelImpl.h"


namespace muq
{
namespace Approximation
{

/**

@class WhiteNoiseKernel
@ingroup CovarianceKernels

This class implements a kernel of the form
\f[
k(x,y) = \sigma^2\delta(x,y)
\f]
where \f$\delta(x,y)=1\f$ if \f$x=y\f$ and \f$0\f$ otherwise.

 */
class WhiteNoiseKernel : public KernelImpl<WhiteNoiseKernel>
{

public:

    WhiteNoiseKernel(unsigned     dim,
                     const double sigma2In) : WhiteNoiseKernel(dim, sigma2In, {0.0, std::numeric_limits<double>::infinity()}){}

    WhiteNoiseKernel(unsigned     dim,
                     const double sigma2In,
                     const Eigen::Vector2d sigmaBounds) : KernelImpl<WhiteNoiseKernel>(dim, 1, 1)
    {
      paramBounds.resize(2,1);
      paramBounds(0) = sigmaBounds[0]; // lower bound on sigma2
      paramBounds(1) = sigmaBounds[1]; // upper bound on sigma2

      cachedParams.resize(1);
      cachedParams(0) = sigma2In;
    };

    virtual ~WhiteNoiseKernel(){};

    template<typename ScalarType1, typename ScalarType2, typename ScalarType3>
    void FillBlockImpl(Eigen::Ref<const Eigen::Matrix<ScalarType1, Eigen::Dynamic, 1>> const& x1,
                       Eigen::Ref<const Eigen::Matrix<ScalarType1, Eigen::Dynamic, 1>> const& x2,
                       Eigen::Ref<const Eigen::Matrix<ScalarType2, Eigen::Dynamic, 1>> const& params,
                       Eigen::Ref<Eigen::Matrix<ScalarType3,Eigen::Dynamic, Eigen::Dynamic>>  block) const
    {
      ScalarType1 squaredDist = (x1-x2).squaredNorm();
      block(0,0) = (squaredDist<1e-14) ? params(0) : 0.0;
    }
};


}
}

#endif
