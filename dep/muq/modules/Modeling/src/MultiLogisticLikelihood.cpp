#include "MUQ/Modeling/MultiLogisticLikelihood.h"

using namespace muq::Modeling;

MultiLogisticLikelihood::MultiLogisticLikelihood(unsigned int numClassesIn,
                                                 Eigen::VectorXi const& dataIn) : ModPiece(numClassesIn*dataIn.size()*Eigen::VectorXi::Ones(1),
                                                                                         Eigen::VectorXi::Ones(1)),
                                                                                numClasses(numClassesIn),
                                                                                data(dataIn)
{
}


void MultiLogisticLikelihood::EvaluateImpl(muq::Modeling::ref_vector<Eigen::VectorXd> const& inputs)
{
  // Construct the scores matrix
  Eigen::VectorXd const& scoresVec = inputs.at(0).get();
  Eigen::Map<const Eigen::MatrixXd> scores(scoresVec.data(), numClasses, data.size());

  // Compute the partition function at each point
  Eigen::RowVectorXd logPartition = scores.array().exp().colwise().sum().log();

  Eigen::VectorXd logP(data.size());
  for(int i=0; i<data.size(); ++i)
    logP(i) = scores(data(i),i) - logPartition(i);

  outputs.resize(1);
  outputs.at(0).resize(1);
  outputs.at(0)(0) = logP.sum();
}

void MultiLogisticLikelihood::GradientImpl(unsigned int                const  outputDimWrt,
                                           unsigned int                const  inputDimWrt,
                                           ref_vector<Eigen::VectorXd> const& inputs,
                                           Eigen::VectorXd             const& sensitivity)
{
  Eigen::VectorXd const& scoresVec = inputs.at(0).get();
  Eigen::Map<const Eigen::MatrixXd> scores(scoresVec.data(), numClasses, data.size());

  // resize the gradient and create a matrix wrapper
  gradient.resize(inputSizes(0));
  Eigen::Map<Eigen::MatrixXd> gradMat(gradient.data(), numClasses, data.size());

  // Compute the partition function at each point
  Eigen::RowVectorXd logPartition = scores.array().exp().colwise().sum().log();

  // Compute the log probabilities
  Eigen::VectorXd logP(data.size());
  for(int i=0; i<data.size(); ++i)
    logP(i) = scores(data(i),i) - logPartition(i);

  double logLikely = logP.sum();

  // First, add the derivatives of the logP terms wrt the scores
  gradMat = Eigen::MatrixXd::Zero(numClasses, data.size());
  for(int i=0; i<data.size(); ++i){
    gradMat(data(i),i) += 1.0;
    gradMat.col(i) -= (1.0/exp(logPartition(i)))*scores.col(i).array().exp().matrix();
  }

  gradMat *= sensitivity(0);
}

void MultiLogisticLikelihood::JacobianImpl(unsigned int                const  outputDimWrt,
                                           unsigned int                const  inputDimWrt,
                                           ref_vector<Eigen::VectorXd> const& input)
{
  jacobian = Gradient(outputDimWrt, inputDimWrt, input, Eigen::VectorXd::Ones(1).eval()).transpose();
}

void MultiLogisticLikelihood::ApplyJacobianImpl(unsigned int                const  outputDimWrt,
                               unsigned int                const  inputDimWrt,
                               ref_vector<Eigen::VectorXd> const& input,
                               Eigen::VectorXd             const& vec)
{
  jacobianAction = Jacobian(outputDimWrt, inputDimWrt, input) * vec;
}
