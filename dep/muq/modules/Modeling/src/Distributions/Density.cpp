#include "MUQ/Modeling/Distributions/Density.h"

using namespace muq::Modeling;

DensityBase::DensityBase(Eigen::VectorXi const& inputSizes) : Distribution(inputSizes(0), inputSizes.tail(inputSizes.size()-1)),
                                                              ModPiece(inputSizes, Eigen::VectorXi::Ones(1))
{};

void DensityBase::EvaluateImpl(ref_vector<Eigen::VectorXd> const& inputs)
{
    outputs.resize(1);
    outputs.at(0) = LogDensity(inputs) * Eigen::VectorXd::Ones(1);
}

void DensityBase::GradientImpl(unsigned int                const  outputDimWrt,
                               unsigned int                const  inputDimWrt,
                               ref_vector<Eigen::VectorXd> const& input,
                               Eigen::VectorXd             const& sensitivity)
{
  gradient = GradLogDensity(inputDimWrt,input)*sensitivity(0);
}

void DensityBase::JacobianImpl(unsigned int                const  outputDimWrt,
                               unsigned int                const  inputDimWrt,
                               ref_vector<Eigen::VectorXd> const& input)
{
  jacobian = GradLogDensity(inputDimWrt,input).transpose();
}


void DensityBase::ApplyJacobianImpl(unsigned int                const  outputDimWrt,
                                    unsigned int                const  inputDimWrt,
                                    ref_vector<Eigen::VectorXd> const& input,
                                    Eigen::VectorXd             const& vec)
{
  jacobianAction = GradLogDensity(inputDimWrt, input).transpose()*vec;
}


Density::Density(std::shared_ptr<Distribution> distIn) : DensityBase(GetInputSizes(distIn)), dist(distIn) {}

double Density::LogDensityImpl(ref_vector<Eigen::VectorXd> const& inputs)
{
  return dist->LogDensity(inputs);
}

Eigen::VectorXd Density::GradLogDensityImpl(unsigned int wrt, ref_vector<Eigen::VectorXd> const& inputs)
{
  return dist->GradLogDensityImpl(wrt, inputs);
}

Eigen::VectorXd Density::SampleImpl(ref_vector<Eigen::VectorXd> const& inputs)
{
  return dist->SampleImpl(inputs);
}

Eigen::VectorXi Density::GetInputSizes(std::shared_ptr<Distribution> distIn)
{
  Eigen::VectorXi output(1+distIn->hyperSizes.size());
  output(0) = distIn->varSize;
  output.tail(distIn->hyperSizes.size()) = distIn->hyperSizes;
  return output;
}
