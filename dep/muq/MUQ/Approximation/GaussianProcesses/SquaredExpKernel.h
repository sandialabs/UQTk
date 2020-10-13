#ifndef SQUAREDEXPKERNEL_H
#define SQUAREDEXPKERNEL_H

#include "MUQ/Approximation/GaussianProcesses/KernelImpl.h"


namespace muq
{
namespace Approximation
{


/**

@class SquaredExpKernel
@ingroup CovarianceKernels
This class implements a kernel of the form
\f[
k(x,y) = \sigma^2 \exp\left(-\frac{1}{2}\frac{|x-y|^2}{L^2}\right)
\f]
for some variance \f$\sigma^2\f$ and lengthscale \f$L\f$.

 */
class SquaredExpKernel : public KernelImpl<SquaredExpKernel>
{

public:

  SquaredExpKernel(unsigned              dimIn,
                   std::vector<unsigned> dimInds,
                   double                sigma2In,
                   double                lengthIn) : SquaredExpKernel(dimIn,
                                                                      dimInds,
                                                                      sigma2In,
                                                                      lengthIn,
                                                                      {0.0, std::numeric_limits<double>::infinity()},
                                                                      {1e-10, std::numeric_limits<double>::infinity()})
  {};

    SquaredExpKernel(unsigned              dimIn,
                     std::vector<unsigned> dimInds,
                     double                sigma2In,
                     double                lengthIn,
                     Eigen::Vector2d       sigmaBounds,
                     Eigen::Vector2d       lengthBounds) : KernelImpl<SquaredExpKernel>(dimIn, dimInds, 1, 2)
    {
      paramBounds.resize(2,2);
      paramBounds(0,0) = sigmaBounds(0);
      paramBounds(1,0) = sigmaBounds(1);
      paramBounds(0,1) = lengthBounds(0);
      paramBounds(1,1) = lengthBounds(1);

      cachedParams.resize(2);
      cachedParams(0) = sigma2In;
      cachedParams(1) = lengthIn;
    };

    SquaredExpKernel(unsigned        dimIn,
                     double          sigma2In,
                     double          lengthIn) : SquaredExpKernel(dimIn,
                                                                  sigma2In,
                                                                  lengthIn,
                                                                  {0.0, std::numeric_limits<double>::infinity()},
                                                                  {1e-10, std::numeric_limits<double>::infinity()} )
    {};

    SquaredExpKernel(unsigned        dimIn,
                     double          sigma2In,
                     double          lengthIn,
                     Eigen::Vector2d sigmaBounds,
                     Eigen::Vector2d lengthBounds) : KernelImpl<SquaredExpKernel>(dimIn, 1, 2)
    {
      paramBounds.resize(2,2);
      paramBounds(0,0) = sigmaBounds(0);
      paramBounds(1,0) = sigmaBounds(1);
      paramBounds(0,1) = lengthBounds(0);
      paramBounds(1,1) = lengthBounds(1);

      cachedParams.resize(2);
      cachedParams(0) = sigma2In;
      cachedParams(1) = lengthIn;
    };

    virtual ~SquaredExpKernel(){};

    template<typename ScalarType1, typename ScalarType2, typename ScalarType3>
    void FillBlockImpl(Eigen::Ref<const Eigen::Matrix<ScalarType1, Eigen::Dynamic, 1>> const& x1,
                       Eigen::Ref<const Eigen::Matrix<ScalarType1, Eigen::Dynamic, 1>> const& x2,
                       Eigen::Ref<const Eigen::Matrix<ScalarType2, Eigen::Dynamic, 1>> const& params,
                       Eigen::Ref<Eigen::Matrix<ScalarType3,Eigen::Dynamic, Eigen::Dynamic>>  block) const
    {
      ScalarType1 squaredDist = (x1-x2).squaredNorm();
      block(0,0) = params(0)*exp(-0.5*squaredDist/(params(1)*params(1)));
    }
};

}
}


#endif
