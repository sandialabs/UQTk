#ifndef LINEARTRANSFORMKERNEL_H
#define LINEARTRANSFORMKERNEL_H

#include "MUQ/Approximation/GaussianProcesses/KernelBase.h"


namespace muq
{
namespace Approximation
{

/**

@class LinearTransformKernel
@ingroup CovarianceKernels

Given another kernel $k_2(x,y)$ and a linear transformation $A$, this class implements a kernel of the form
\f[
k(x,y) = A * k_2(x,y) * A^T.
\f]
 */
class LinearTransformKernel : public KernelBase
{

public:


    LinearTransformKernel(Eigen::MatrixXd       const& Ain,
		                      std::shared_ptr<KernelBase>  Kin) : KernelBase(Kin->inputDim, Ain.rows(), Kin->numParams),
                                                              A(Ain), K(Kin)
    {
	     assert(Ain.cols() == Kin->coDim);
       cachedParams = Kin->GetParams();
    };

    virtual ~LinearTransformKernel(){};


    virtual void FillBlock(Eigen::Ref<const Eigen::VectorXd> const& x1,
                           Eigen::Ref<const Eigen::VectorXd> const& x2,
                           Eigen::Ref<const Eigen::VectorXd> const& params,
                           Eigen::Ref<Eigen::MatrixXd>              block) const override
    {
      Eigen::MatrixXd temp(K->coDim,K->coDim);
      K->FillBlock(x1,x2, params, temp);
      block = A*temp*A.transpose();
    }

    virtual void FillPosDerivBlock(Eigen::Ref<const Eigen::VectorXd> const& x1,
                                   Eigen::Ref<const Eigen::VectorXd> const& x2,
                                   Eigen::Ref<const Eigen::VectorXd> const& params,
                                   std::vector<int>                  const& wrts,
                                   Eigen::Ref<Eigen::MatrixXd>              block) const override
    {
      Eigen::MatrixXd temp(K->coDim,K->coDim);
      K->FillPosDerivBlock(x1,x2, params, wrts, temp);
      block = A*temp*A.transpose();
    }

    virtual std::shared_ptr<KernelBase> Clone() const override{return std::make_shared<LinearTransformKernel>(A,K);};

private:
    Eigen::MatrixXd A;
    std::shared_ptr<KernelBase> K;

};


template<typename KernelType>
LinearTransformKernel TransformKernel(Eigen::MatrixXd const& A, KernelType const& K)
{
    return LinearTransformKernel(A,K.Clone());
}

template<typename KernelType,
         typename = typename std::enable_if<std::is_base_of<KernelBase, KernelType>::value>::type>
LinearTransformKernel operator*(Eigen::MatrixXd const& A, KernelType const& kernel)
{
    return LinearTransformKernel(A, kernel.Clone());
}

}
}

#endif
