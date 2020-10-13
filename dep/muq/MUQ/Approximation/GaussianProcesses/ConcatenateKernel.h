#ifndef CONCATENATEKERNEL_H
#define CONCATENATEKERNEL_H

#include "MUQ/Approximation/GaussianProcesses/KernelBase.h"


namespace muq
{
namespace Approximation
{

    /**
       @class ConcatenateKernel
       @ingroup CovarianceKernels
       @brief Concatenates two kernels together.
       @details Let \f$k_1(x,x^\prime)\f$ and \f$k_2(x,x^\prime)\f$ be two different covariance kernels with the same inputs.  This class describes a concatenated kernel of the form
\f[
k(x,x^\prime) = \left[\begin{array}{cc}k_1(x,x^\prime) & 0\\ 0 & k_2(x,x^\prime)\end{array}\right].
\f]
    */
    class ConcatenateKernel : public KernelBase
    {

    public:

        ConcatenateKernel(std::shared_ptr<KernelBase> const& kernel1In,
                          std::shared_ptr<KernelBase> const& kernel2In) : ConcatenateKernel(std::vector<std::shared_ptr<KernelBase>>({kernel1In, kernel2In})){};

        ConcatenateKernel(std::vector<std::shared_ptr<KernelBase>> const& kernelsIn);

        virtual ~ConcatenateKernel() = default;

        virtual std::shared_ptr<KernelBase> Clone() const override{return std::make_shared<ConcatenateKernel>(kernels);};

        virtual void FillBlock(Eigen::Ref<const Eigen::VectorXd> const& x1,
                               Eigen::Ref<const Eigen::VectorXd> const& x2,
                               Eigen::Ref<const Eigen::VectorXd> const& params,
                               Eigen::Ref<Eigen::MatrixXd>              block) const override;

        virtual void FillPosDerivBlock(Eigen::Ref<const Eigen::VectorXd> const& x1,
                                       Eigen::Ref<const Eigen::VectorXd> const& x2,
                                       Eigen::Ref<const Eigen::VectorXd> const& params,
                                       std::vector<int>                  const& wrts,
                                       Eigen::Ref<Eigen::MatrixXd>              block) const override;

        // template<typename VecType1, typename VecType2, typename MatType>
        //     inline void GetDerivative(VecType1 const& x1, VecType2 const& x2, int wrt, MatType & derivs) const
        // {
        //     assert(wrt < this->numParams);
        //
        //     // Initialize the derivative matrix to 0
        //     for(int j=0; j<derivs.cols(); ++j)
        //     {
        //         for(int i=0; i<derivs.rows(); ++i)
        //             derivs(i,j) = 0.0;
        //     }
        //
        //     if(wrt < kernel1.numParams )
        //     {
        //         auto block = GetBlock(derivs, 0, 0, kernel1.coDim, kernel1.coDim);
        //         return kernel1.GetDerivative(x1, x2, wrt, block);
        //     }
        //     else
        //     {
        //         auto block = GetBlock(derivs, kernel1.coDim, kernel1.coDim, kernel2.coDim, kernel2.coDim);
        //         return kernel2.GetDerivative(x1, x2, wrt-kernel1.numParams, block);
        //     }
        // }

    private:

        static unsigned int CountCoDims(std::vector<std::shared_ptr<KernelBase>> kernels);
        static unsigned int CountParams(std::vector<std::shared_ptr<KernelBase>> kernels);

        std::vector<std::shared_ptr<KernelBase>> kernels;
    };



template<typename KernelType1, typename KernelType2>
ConcatenateKernel Concatenate(KernelType1 const& kernel1, KernelType2 const& kernel2)
{
    return ConcatenateKernel(kernel1.Clone(), kernel2.Clone());
}




} // namespace Approximation
} // namespace muq


#endif
