#ifndef CONSTANTKERNEL_H
#define CONSTANTKERNEL_H

#include "MUQ/Approximation/GaussianProcesses/KernelImpl.h"


namespace muq
{
namespace Approximation
{

/**

@class ConstantKernel
@ingroup CovarianceKernels

This class implements a constant kernel of the form
\f[
k(x,y) = \sigma^2
\f]
where, \f$\sigma^2\f$ is the variance.

 */
class ConstantKernel : public KernelImpl<ConstantKernel>
{

public:

    ConstantKernel(unsigned              dim,
	                 const double          sigma2In,
                   const Eigen::Vector2d sigmaBounds = {0.0, std::numeric_limits<double>::infinity()}) : ConstantKernel(dim, sigma2In*Eigen::MatrixXd::Ones(1,1), sigmaBounds){};

    ConstantKernel(unsigned              dim,
		               std::vector<unsigned> dimInds,
	                 const double          sigma2In,
                   const Eigen::Vector2d sigmaBounds = {0.0, std::numeric_limits<double>::infinity()}) : ConstantKernel(dim, dimInds, sigma2In*Eigen::MatrixXd::Ones(1,1), sigmaBounds){};


    ConstantKernel(unsigned               dim,
	                 Eigen::MatrixXd const& sigma2In,
                   const Eigen::Vector2d  sigmaBounds = {0.0, std::numeric_limits<double>::infinity()}) : KernelImpl<ConstantKernel>(dim, sigma2In.rows(), GetNumParams(sigma2In))
    {
    	paramBounds.resize(2,1);
    	paramBounds(0,0) = sigmaBounds(0);
    	paramBounds(1,0) = sigmaBounds(1);

      cachedParams.resize(numParams);
      int ind = 0;
      for(int i=0; i<sigma2In.rows(); ++i){
        for(int j=0; j<=i; ++j){
          cachedParams(ind) = sigma2In(i,j);
          ind++;
        }
      }
    };

    ConstantKernel(unsigned               dim,
		               std::vector<unsigned>  dimInds,
	                 Eigen::MatrixXd const& sigma2In,
                   const Eigen::Vector2d  sigmaBounds = {0.0, std::numeric_limits<double>::infinity()}) : KernelImpl<ConstantKernel>(dim, dimInds, sigma2In.rows(), GetNumParams(sigma2In))
    {
    	paramBounds.resize(2,1);
    	paramBounds(0,0) = sigmaBounds(0);
    	paramBounds(1,0) = sigmaBounds(1);

      cachedParams.resize(numParams);
      int ind = 0;
      for(int i=0; i<sigma2In.rows(); ++i){
        for(int j=0; j<=i; ++j){
          cachedParams(ind) = sigma2In(i,j);
          ind++;
        }
      }
    };

    virtual ~ConstantKernel(){};

    template<typename ScalarType1, typename ScalarType2, typename ScalarType3>
    void FillBlockImpl(Eigen::Ref<const Eigen::Matrix<ScalarType1, Eigen::Dynamic, 1>> const& x1,
                       Eigen::Ref<const Eigen::Matrix<ScalarType1, Eigen::Dynamic, 1>> const& x2,
                       Eigen::Ref<const Eigen::Matrix<ScalarType2, Eigen::Dynamic, 1>> const& params,
                       Eigen::Ref<Eigen::Matrix<ScalarType3,Eigen::Dynamic, Eigen::Dynamic>>  block) const
    {
      int ind = 0;
      for(int i=0; i<coDim; ++i){
        for(int j=0; j<=i; ++j){
          block(i,j) = params(ind);
          block(j,i) = block(i,j);
          ind++;
        }
      }
    }

private:

    static unsigned GetNumParams(Eigen::MatrixXd const& cov)
    {
	     return 0.5*cov.rows()*(cov.rows()+1);
    }
};

}
}


#endif
