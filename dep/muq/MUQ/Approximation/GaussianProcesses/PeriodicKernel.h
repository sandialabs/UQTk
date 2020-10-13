#ifndef PERIODICKERNEL_H
#define PERIODICKERNEL_H

#include "MUQ/Approximation/GaussianProcesses/KernelImpl.h"


namespace muq
{
namespace Approximation
{

    /**

@class PeriodicKernel
@ingroup CovarianceKernels

\f[
k(x,y) = \exp\left[-\frac{2}{L^2} \sin^2\left(\frac{\pi}{p}(x-y)\right)\right]

 */
class PeriodicKernel : public KernelImpl<PeriodicKernel>
{

public:

  PeriodicKernel(unsigned dim,
                 std::vector<unsigned> dimInds,
                 const double sigma2In,
                 const double lengthIn,
                 const double periodIn) : PeriodicKernel(dim,
                                                         dimInds,
                                                         sigma2In,
                                                         lengthIn,
                                                         periodIn,
                                                         {0.0, std::numeric_limits<double>::infinity()},
                                                         {1e-10, std::numeric_limits<double>::infinity()},
                                                         {1e-10, std::numeric_limits<double>::infinity()}){};

     PeriodicKernel(unsigned dim,
            		    std::vector<unsigned> dimInds,
            		    const double sigma2In,
            		    const double lengthIn,
            		    const double periodIn,
                    const Eigen::Vector2d sigmaBounds,
	                  const Eigen::Vector2d lengthBounds,
		                const Eigen::Vector2d periodBounds);


  PeriodicKernel(unsigned dim,
          		   const double sigma2In,
          		   const double lengthIn,
          		   const double periodIn) : PeriodicKernel(dim,
                                                         sigma2In,
                                                         lengthIn,
                                                         periodIn,
                                                         {0.0, std::numeric_limits<double>::infinity()},
                                                         {1e-10, std::numeric_limits<double>::infinity()},
                                                         {1e-10, std::numeric_limits<double>::infinity()}){};

    PeriodicKernel(unsigned dim,
            		   const double sigma2In,
            		   const double lengthIn,
            		   const double periodIn,
                   const Eigen::Vector2d sigmaBounds,
            	     const Eigen::Vector2d lengthBounds,
            	     const Eigen::Vector2d periodBounds);


    virtual ~PeriodicKernel(){};

    template<typename ScalarType1, typename ScalarType2, typename ScalarType3>
    void FillBlockImpl(Eigen::Ref<const Eigen::Matrix<ScalarType1, Eigen::Dynamic, 1>> const& x1,
                       Eigen::Ref<const Eigen::Matrix<ScalarType1, Eigen::Dynamic, 1>> const& x2,
                       Eigen::Ref<const Eigen::Matrix<ScalarType2, Eigen::Dynamic, 1>> const& params,
                       Eigen::Ref<Eigen::Matrix<ScalarType3,Eigen::Dynamic, Eigen::Dynamic>>  block) const
    {
      ScalarType1 dist = (x1-x2).norm();
      block(0,0) = params(0) * exp(-2.0 * pow(sin(pi*dist/params(2)),2.0) / pow(params(1),2.0));
    }


    virtual std::tuple<std::shared_ptr<muq::Modeling::LinearSDE>, std::shared_ptr<muq::Modeling::LinearOperator>, Eigen::MatrixXd> GetStateSpace(boost::property_tree::ptree sdeOptions = boost::property_tree::ptree()) const override;

private:
    const double pi = 4.0 * atan(1.0); //boost::math::constants::pi<double>();

    void SetupBounds(Eigen::Vector2d const& sigmaBounds,
                     Eigen::Vector2d const& lengthBounds,
                     Eigen::Vector2d const& periodBounds);
};

}
}


#endif
