#include "MUQ/Approximation/GaussianProcesses/SumKernel.h"

#include "MUQ/Modeling/LinearSDE.h"

#include "MUQ/Modeling/LinearAlgebra/LinearOperator.h"
#include "MUQ/Modeling/LinearAlgebra/BlockRowOperator.h"

using namespace muq::Approximation;

SumKernel::SumKernel(std::shared_ptr<KernelBase> kernel1In,
                     std::shared_ptr<KernelBase> kernel2In) : KernelBase(kernel1In->inputDim,
                                                                         kernel1In->coDim,
                                                                         kernel1In->numParams + kernel2In->numParams),
                                                   kernel1(kernel1In),
                                                   kernel2(kernel2In)
{
  assert(kernel2->inputDim == kernel2->inputDim);
  assert(kernel1->coDim == kernel2->coDim);
  cachedParams.resize(numParams);
  cachedParams.head(kernel1->numParams) = kernel1->GetParams();
  cachedParams.tail(kernel2->numParams) = kernel2->GetParams();
};


void SumKernel::FillBlock(Eigen::Ref<const Eigen::VectorXd> const& x1,
                       Eigen::Ref<const Eigen::VectorXd> const& x2,
                       Eigen::Ref<const Eigen::VectorXd> const& params,
                       Eigen::Ref<Eigen::MatrixXd>              block) const
{

    kernel1->FillBlock(x1,x2, params.head(kernel1->numParams), block);

    Eigen::MatrixXd temp(coDim, coDim);
    kernel2->FillBlock(x1,x2, params.tail(kernel2->numParams), temp);

    block += temp;
}


void SumKernel::FillPosDerivBlock(Eigen::Ref<const Eigen::VectorXd> const& x1,
                               Eigen::Ref<const Eigen::VectorXd> const& x2,
                               Eigen::Ref<const Eigen::VectorXd> const& params,
                               std::vector<int>                  const& wrts,
                               Eigen::Ref<Eigen::MatrixXd>              block) const
{
  kernel1->FillPosDerivBlock(x1,x2, params.head(kernel1->numParams), wrts, block);

  Eigen::MatrixXd temp(coDim, coDim);
  kernel2->FillPosDerivBlock(x1,x2, params.tail(kernel2->numParams), wrts, temp);

  block += temp;
}

std::tuple<std::shared_ptr<muq::Modeling::LinearSDE>, std::shared_ptr<muq::Modeling::LinearOperator>, Eigen::MatrixXd> SumKernel::GetStateSpace(boost::property_tree::ptree sdeOptions) const
{
    std::vector<std::shared_ptr<muq::Modeling::LinearSDE>> sdes(2);
    std::vector<std::shared_ptr<muq::Modeling::LinearOperator>> linOps(2);
    Eigen::MatrixXd pinf1, pinf2;

    // Get the statespace information from each component
    std::tie(sdes.at(0), linOps.at(0), pinf1) = kernel1->GetStateSpace(sdeOptions);
    std::tie(sdes.at(1), linOps.at(1), pinf2) = kernel2->GetStateSpace(sdeOptions);

    // Concantenate the sdes
    auto newSDE = muq::Modeling::LinearSDE::Concatenate(sdes,sdeOptions);

    // Concatenate the linear operators
    auto newObsOp = std::make_shared<muq::Modeling::BlockRowOperator>(linOps);

    // Set up the combined stationary covariance
    Eigen::MatrixXd newPinf = Eigen::MatrixXd::Zero(pinf1.rows() + pinf2.rows(), pinf1.cols() + pinf2.cols());
    newPinf.block(0,0,pinf1.rows(),pinf1.cols()) = pinf1;
    newPinf.block(pinf1.rows(),pinf1.cols(), pinf2.rows(), pinf2.cols()) = pinf2;

    return std::make_tuple(newSDE, newObsOp, newPinf);

}


std::shared_ptr<SumKernel> muq::Approximation::operator+(std::shared_ptr<KernelBase> k1, std::shared_ptr<KernelBase> k2)
{
  return std::make_shared<SumKernel>(k1, k2);
}
