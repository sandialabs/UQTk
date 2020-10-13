#include "MUQ/Approximation/GaussianProcesses/ConcatenateKernel.h"

using namespace muq::Approximation;

ConcatenateKernel::ConcatenateKernel(std::vector<std::shared_ptr<KernelBase>> const& kernelsIn) : KernelBase(kernelsIn.at(0)->inputDim,
                                                                                                             CountCoDims(kernelsIn),
                                                                                                             CountParams(kernelsIn)),
                                                                                                  kernels(kernelsIn)
{
  // Make sure all the input sizes are the same
  for(int i=1; i<kernels.size(); ++i)
    assert(kernels.at(i)->inputDim == kernels.at(0)->inputDim);

  // Set all the parameters from what was in the kernels before
  cachedParams.resize(numParams);
  int paramCnt = 0;
  for(int i=0; i<kernels.size(); ++i){
    cachedParams.segment(paramCnt, kernels.at(i)->numParams) = kernels.at(i)->GetParams();
    paramCnt += kernels.at(i)->numParams;
  }
}


void ConcatenateKernel::FillBlock(Eigen::Ref<const Eigen::VectorXd> const& x1,
                       Eigen::Ref<const Eigen::VectorXd> const& x2,
                       Eigen::Ref<const Eigen::VectorXd> const& params,
                       Eigen::Ref<Eigen::MatrixXd>              block) const
{
  block = Eigen::MatrixXd::Zero(coDim, coDim);
  int paramCnt = 0;
  int codimCnt = 0;

  for(int i=0; i<kernels.size(); ++i){

    kernels.at(i)->FillBlock(x1,
                             x2,
                             params.segment(paramCnt, kernels.at(i)->numParams),
                             block.block(codimCnt,codimCnt,kernels.at(i)->coDim, kernels.at(i)->coDim));

    codimCnt += kernels.at(i)->coDim;
    paramCnt += kernels.at(i)->numParams;
  }
}

void ConcatenateKernel::FillPosDerivBlock(Eigen::Ref<const Eigen::VectorXd> const& x1,
                               Eigen::Ref<const Eigen::VectorXd> const& x2,
                               Eigen::Ref<const Eigen::VectorXd> const& params,
                               std::vector<int>                  const& wrts,
                               Eigen::Ref<Eigen::MatrixXd>              block) const
{
  block = Eigen::MatrixXd::Zero(coDim, coDim);

  int paramCnt = 0;
  int codimCnt = 0;
  for(int i=0; i<kernels.size(); ++i){

      kernels.at(i)->FillPosDerivBlock(x1,
                                       x2,
                                       params.segment(paramCnt, kernels.at(i)->numParams),
                                       wrts,
                                       block.block(codimCnt,codimCnt,kernels.at(i)->coDim, kernels.at(i)->coDim));

      codimCnt += kernels.at(i)->coDim;
      paramCnt += kernels.at(i)->numParams;
  }
}

unsigned int ConcatenateKernel::CountCoDims(std::vector<std::shared_ptr<KernelBase>> kernelsIn)
{
  int cnt = 0;
  for(auto& kernel : kernelsIn)
    cnt += kernel->coDim;
  return cnt;
}
unsigned int ConcatenateKernel::CountParams(std::vector<std::shared_ptr<KernelBase>> kernelsIn)
{
  int cnt = 0;
  for(auto& kernel : kernelsIn)
    cnt += kernel->numParams;
  return cnt;
}
