#ifndef PRODUCTKERNEL_H
#define PRODUCTKERNEL_H

#include "MUQ/Approximation/GaussianProcesses/KernelBase.h"
#include "MUQ/Approximation/GaussianProcesses/KernelImpl.h"
#include "MUQ/Approximation/GaussianProcesses/PeriodicKernel.h"

#include "MUQ/Modeling/LinearSDE.h"

#include "MUQ/Modeling/LinearAlgebra/KroneckerProductOperator.h"
#include "MUQ/Modeling/LinearAlgebra/BlockDiagonalOperator.h"
#include "MUQ/Modeling/LinearAlgebra/BlockRowOperator.h"
#include "MUQ/Modeling/LinearAlgebra/SumOperator.h"
#include "MUQ/Modeling/LinearAlgebra/IdentityOperator.h"

#include "MUQ/Utilities/Exceptions.h"

namespace muq
{
namespace Approximation
{

/**

@class ProductKernel
#ingroup CovarianceKernels
\f[
k(x,y) = k_1(x,y)*k_2(x,y)
\f]

 */
class ProductKernel : public KernelBase
{

public:
    ProductKernel(std::shared_ptr<KernelBase> kernel1In,
                  std::shared_ptr<KernelBase> kernel2In);

    virtual ~ProductKernel(){};

    virtual void FillBlock(Eigen::Ref<const Eigen::VectorXd> const& x1,
                           Eigen::Ref<const Eigen::VectorXd> const& x2,
                           Eigen::Ref<const Eigen::VectorXd> const& params,
                           Eigen::Ref<Eigen::MatrixXd>              block) const override;

    virtual void FillPosDerivBlock(Eigen::Ref<const Eigen::VectorXd> const& x1,
                                   Eigen::Ref<const Eigen::VectorXd> const& x2,
                                   Eigen::Ref<const Eigen::VectorXd> const& params,
                                   std::vector<int>                  const& wrts,
                                   Eigen::Ref<Eigen::MatrixXd>              block) const override
    {


      Eigen::MatrixXd eval1(coDim,coDim);
      kernel1->FillBlock(x1,x2,params.head(kernel1->numParams),eval1);

      Eigen::MatrixXd eval2(coDim,coDim);
      kernel2->FillBlock(x1,x2,params.tail(kernel2->numParams),eval2);

      Eigen::MatrixXd deriv1_i(coDim,coDim);
      kernel1->FillPosDerivBlock(x1,x2,params.head(kernel1->numParams),{wrts.at(0)}, deriv1_i);

      Eigen::MatrixXd deriv2_i(coDim,coDim);
      kernel2->FillPosDerivBlock(x1,x2,params.tail(kernel2->numParams), {wrts.at(0)}, deriv2_i);

      if(wrts.size()==1){
        block = (deriv2_i.array() * eval1.array() + deriv1_i.array()*eval2.array()).matrix();

      }else if(wrts.size()==2){
        Eigen::MatrixXd deriv1_j(coDim,coDim);
        kernel1->FillPosDerivBlock(x1,x2,params.head(kernel1->numParams),{wrts.at(1)}, deriv1_j);

        Eigen::MatrixXd deriv2_j(coDim,coDim);
        kernel2->FillPosDerivBlock(x1,x2,params.tail(kernel2->numParams), {wrts.at(1)}, deriv2_j);

        Eigen::MatrixXd secDeriv1_ij(coDim,coDim);
        kernel1->FillPosDerivBlock(x1,x2,params.head(kernel1->numParams), wrts, secDeriv1_ij);

        Eigen::MatrixXd secDeriv2_ij(coDim,coDim);
        kernel2->FillPosDerivBlock(x1,x2,params.tail(kernel2->numParams), wrts, secDeriv2_ij);

        block = (secDeriv1_ij.array() * eval2.array() + deriv1_i.array()*deriv2_j.array() + deriv1_j.array()*deriv2_i.array() + eval1.array()*secDeriv2_ij.array()).matrix();

      }else{
        assert(false);
      }
    }

    virtual std::shared_ptr<KernelBase> Clone() const override{return std::make_shared<ProductKernel>(kernel1,kernel2);};

  //   template<typename VecType1, typename VecType2, typename MatType>
  //   inline void GetDerivative(VecType1 const& x1, VecType2 const& x2, int wrt, MatType & derivs) const
  //   {
	// assert(wrt < this->numParams);
  //
	// Eigen::MatrixXd temp1(kernel1.coDim, kernel1.coDim);
	// Eigen::MatrixXd temp2(kernel2.coDim, kernel2.coDim);
  //
	// if(wrt < kernel1.numParams )
	// {
	//     kernel1.GetDerivative(x1, x2, wrt, temp1);
	//     kernel2.EvaluateImpl(x1,x2, temp2);
	// }
	// else
	// {
	//     kernel1.EvaluateImpl(x1,x2, temp1);
	//     kernel2.GetDerivative(x1,x2, wrt-kernel1.numParams, temp2);
	// }
  //
	// if(kernel1.coDim==kernel2.coDim)
	// {
	//     derivs = Eigen::MatrixXd(temp1.array() * temp2.array());
	// }
	// else if(kernel1.coDim==1)
	// {
	//     derivs = temp1(0,0) * temp2;
	// }
	// else if(kernel2.coDim==1)
	// {
  //           derivs = temp2(0,0) * temp1;
	// }
	// else
	// {
	//     std::cerr << "\nERROR: Something unexpected happened with the dimensions of the kernels in this product.\n";
	//     assert(false);
	// }

    //}

    virtual std::vector<std::shared_ptr<KernelBase>> GetSeperableComponents() override;

    virtual std::tuple<std::shared_ptr<muq::Modeling::LinearSDE>, std::shared_ptr<muq::Modeling::LinearOperator>, Eigen::MatrixXd> GetStateSpace(boost::property_tree::ptree sdeOptions=boost::property_tree::ptree()) const override;

protected:
    std::shared_ptr<KernelBase> kernel1;
    std::shared_ptr<KernelBase> kernel2;

    std::tuple<std::shared_ptr<muq::Modeling::LinearSDE>, std::shared_ptr<muq::Modeling::LinearOperator>, Eigen::MatrixXd> GetProductStateSpace(std::shared_ptr<PeriodicKernel> const& kernel1,
                                                                                                                                                 	              std::shared_ptr<KernelBase>     const& kernel2,
                                                                                                                                                                boost::property_tree::ptree sdeOptions) const;

};


// Operator overload
template<typename KernelType1, typename KernelType2,
         typename = typename std::enable_if<std::is_base_of<KernelBase, KernelType1>::value && std::is_base_of<KernelBase, KernelType2>::value, KernelType1>::type>
ProductKernel operator*(KernelType1 const& k1, KernelType2 const& k2)
{
  return ProductKernel(k1.Clone(), k2.Clone());
}

std::shared_ptr<ProductKernel> operator*(std::shared_ptr<KernelBase> k1, std::shared_ptr<KernelBase> k2);

} // namespace Approximation
}// namespace muq


#endif
