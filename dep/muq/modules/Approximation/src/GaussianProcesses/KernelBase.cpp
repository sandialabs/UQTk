#include "MUQ/Approximation/GaussianProcesses/KernelBase.h"

using namespace muq::Approximation;

Eigen::MatrixXd KernelBase::Evaluate(Eigen::VectorXd const& x1,
																     Eigen::VectorXd const& x2) const
{
	Eigen::MatrixXd output(coDim, coDim);

  // Because Eigen doesn't support arbitrary indexing yet, we have to create
  // A new vector containing the individual components of x1 and x2
  Eigen::VectorXd x1Part(dimInds.size());
  Eigen::VectorXd x2Part(dimInds.size());

  for(int i=0; i<dimInds.size(); ++i){
    x1Part(i) = x1(dimInds.at(i));
    x2Part(i) = x2(dimInds.at(i));
  }

	FillBlock(x1Part, x2Part, cachedParams, output);
	return output;
}

Eigen::MatrixXd KernelBase::BuildCovariance(Eigen::MatrixXd const& x) const
{
	Eigen::MatrixXd output(coDim*x.cols(), coDim*x.cols());
	FillCovariance(x,output);
	return output;
}

Eigen::MatrixXd KernelBase::BuildCovariance(Eigen::MatrixXd const& x1,
	                                          Eigen::MatrixXd const& x2) const
{
	Eigen::MatrixXd output(coDim*x1.cols(), coDim*x2.cols());
	FillCovariance(x1,x2,output);
	return output;
}

void KernelBase::FillCovariance(Eigen::MatrixXd      const& x,
                                Eigen::Ref<Eigen::MatrixXd> output) const
{
	assert(output.rows()==x.cols()*coDim);
	assert(output.cols()==x.cols()*coDim);

	const int numPts = x.cols();

  Eigen::MatrixXd xParts(dimInds.size(), x.cols());
  for(int i=0; i<dimInds.size(); ++i)
    xParts.row(i) = x.row(dimInds.at(i));

	#pragma omp parallel for
	for(int i=0; i<numPts; ++i){

		for(int j=0; j<i; ++j){
			FillBlock(xParts.col(i), xParts.col(j), cachedParams, output.block(i*coDim, j*coDim, coDim, coDim));
			output.block(j*coDim, i*coDim, coDim, coDim) = output.block(i*coDim, j*coDim, coDim, coDim).transpose();
		}

		FillBlock(xParts.col(i), xParts.col(i), cachedParams, output.block(i*coDim,i*coDim,coDim,coDim));
	}

}

void KernelBase::FillCovariance(Eigen::MatrixXd      const& x1,
																Eigen::MatrixXd      const& x2,
															  Eigen::Ref<Eigen::MatrixXd> output) const
{
	int numRows = x1.cols();
	int numCols = x2.cols();

	assert(output.rows()==numRows*coDim);
	assert(output.cols()==numCols*coDim);

  Eigen::MatrixXd x1Parts(dimInds.size(), x1.cols());
  Eigen::MatrixXd x2Parts(dimInds.size(), x2.cols());
  for(int i=0; i<dimInds.size(); ++i){
    x1Parts.row(i) = x1.row(dimInds.at(i));
    x2Parts.row(i) = x2.row(dimInds.at(i));
  }
	#pragma omp parallel for
	for(int i=0; i<numRows; ++i){
		for(int j=0; j<numCols; ++j)
			FillBlock(x1Parts.col(i), x2Parts.col(j), cachedParams, output.block(i*coDim, j*coDim, coDim, coDim));
	}
}

void KernelBase::FillDerivCovariance(Eigen::MatrixXd             const& x1,
                                     Eigen::MatrixXd             const& x2,
                                     std::vector<int>            const& wrts,
                                     Eigen::Ref<Eigen::MatrixXd>        output) const
{
  int numRows = x1.cols();
	int numCols = x2.cols();

  output.resize(numRows*coDim,numCols*coDim);

  Eigen::MatrixXd x1Parts(dimInds.size(), x1.cols());
  Eigen::MatrixXd x2Parts(dimInds.size(), x2.cols());
  for(int i=0; i<dimInds.size(); ++i){
    x1Parts.row(i) = x1.row(dimInds.at(i));
    x2Parts.row(i) = x2.row(dimInds.at(i));
  }

	#pragma omp parallel for
	for(int i=0; i<numRows; ++i){
		for(int j=0; j<numCols; ++j){
			FillPosDerivBlock(x1Parts.col(i), x2Parts.col(j), cachedParams, wrts, output.block(i*coDim, j*coDim, coDim, coDim));
    }
	}
}

Eigen::MatrixXd KernelBase::GetPosDerivative(Eigen::VectorXd  const& x1,
                                             Eigen::VectorXd  const& x2,
                                             std::vector<int> const& wrts) const
{
  int numRows = x1.cols();
	int numCols = x2.cols();

  Eigen::MatrixXd output(numRows*coDim,numCols*coDim);

  FillDerivCovariance(x1,x2,wrts,output);

  return output;
}
