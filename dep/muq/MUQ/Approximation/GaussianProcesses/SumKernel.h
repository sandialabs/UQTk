#ifndef SUMKERNEL_H
#define SUMKERNEL_H

#include "MUQ/Approximation/GaussianProcesses/KernelBase.h"


namespace muq
{
namespace Approximation
{

/**

@class SumKernel
@ingroup CovarianceKernels
\f[
k)x,y) = k_1(x,y) + k_2(x,y)
\f]

 */
class SumKernel : public KernelBase
{

public:
    SumKernel(std::shared_ptr<KernelBase> kernel1In,
              std::shared_ptr<KernelBase> kernel2In);

    virtual ~SumKernel(){};

    virtual std::shared_ptr<KernelBase> Clone() const override{return std::make_shared<SumKernel>(kernel1,kernel2);};

    virtual void FillBlock(Eigen::Ref<const Eigen::VectorXd> const& x1,
                           Eigen::Ref<const Eigen::VectorXd> const& x2,
                           Eigen::Ref<const Eigen::VectorXd> const& params,
                           Eigen::Ref<Eigen::MatrixXd>              block) const override;


    virtual void FillPosDerivBlock(Eigen::Ref<const Eigen::VectorXd> const& x1,
                                   Eigen::Ref<const Eigen::VectorXd> const& x2,
                                   Eigen::Ref<const Eigen::VectorXd> const& params,
                                   std::vector<int>                  const& wrts,
                                   Eigen::Ref<Eigen::MatrixXd>              block) const override;


    virtual std::tuple<std::shared_ptr<muq::Modeling::LinearSDE>, std::shared_ptr<muq::Modeling::LinearOperator>, Eigen::MatrixXd> GetStateSpace(boost::property_tree::ptree sdeOptions=boost::property_tree::ptree()) const override;

private:
    std::shared_ptr<KernelBase>  kernel1;
    std::shared_ptr<KernelBase> kernel2;
};


template<typename KernelType1, typename KernelType2, typename = typename std::enable_if<std::is_base_of<KernelBase, KernelType1>::value && std::is_base_of<KernelBase, KernelType2>::value, KernelType1>::type>
SumKernel operator+(KernelType1 const& k1, KernelType2 const& k2)
{
  return SumKernel(k1.Clone(), k2.Clone());
}

std::shared_ptr<SumKernel> operator+(std::shared_ptr<KernelBase> k1, std::shared_ptr<KernelBase> k2);

} // namespace Approximation
} // namespace muq


#endif
