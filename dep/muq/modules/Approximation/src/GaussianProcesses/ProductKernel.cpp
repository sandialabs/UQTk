#include "MUQ/Approximation/GaussianProcesses/ProductKernel.h"

using namespace muq::Approximation;

ProductKernel::ProductKernel(std::shared_ptr<KernelBase> kernel1In,
							std::shared_ptr<KernelBase> kernel2In) : KernelBase(kernel1In->inputDim,
																																	std::max(kernel1In->coDim, kernel2In->coDim),
																																	kernel1In->numParams + kernel2In->numParams),
																											 kernel1(kernel1In),
																											 kernel2(kernel2In)
{
	 assert((kernel1->coDim==kernel2->coDim) | (kernel1->coDim==1) | (kernel2->coDim==1));

	 cachedParams.resize(numParams);
	 cachedParams.head(kernel1In->numParams) = kernel1In->GetParams();
	 cachedParams.tail(kernel2In->numParams) = kernel2In->GetParams();

};


void ProductKernel::FillBlock(Eigen::Ref<const Eigen::VectorXd> const& x1,
											        Eigen::Ref<const Eigen::VectorXd> const& x2,
											        Eigen::Ref<const Eigen::VectorXd> const& params,
											        Eigen::Ref<Eigen::MatrixXd>              block) const
{
	Eigen::MatrixXd temp1(kernel1->coDim, kernel1->coDim);
	Eigen::MatrixXd temp2(kernel2->coDim, kernel2->coDim);

	kernel1->FillBlock(x1, x2, params.head(kernel1->numParams), temp1);
	kernel2->FillBlock(x1, x2, params.tail(kernel2->numParams), temp2);

	if(kernel1->coDim==kernel2->coDim)
	{
			block = Eigen::MatrixXd(temp1.array() * temp2.array());
	}
	else if(kernel1->coDim==1)
	{
			block = temp1(0,0)*temp2;
	}
	else if(kernel2->coDim==1)
	{
			block = temp2(0,0)*temp1;
	}
	else
	{
			std::cerr << "\nERROR: Something unexpected happened with the dimensions of the kernels in this product.\n";
			assert(false);
	}
}


std::vector<std::shared_ptr<KernelBase>> ProductKernel::GetSeperableComponents()
{
		// Check if the dimensions of the components are distinct
		bool isSeperable = true;
		for(unsigned leftDim : kernel1->dimInds)
		{
				for(unsigned rightDim : kernel2->dimInds)
				{
						if(leftDim==rightDim)
						{
								isSeperable = false;
								break;
						}
				}

				if(!isSeperable)
						break;
		}

		if(isSeperable)
		{
				std::vector<std::shared_ptr<KernelBase>> output, output2;
				output = kernel1->GetSeperableComponents();
				output2 = kernel2->GetSeperableComponents();
				output.insert(output.end(), output2.begin(), output2.end());

				return output;
		}
		else
		{
				return std::vector<std::shared_ptr<KernelBase>>(1, this->Clone());
		}

};



std::tuple<std::shared_ptr<muq::Modeling::LinearSDE>, std::shared_ptr<muq::Modeling::LinearOperator>, Eigen::MatrixXd> ProductKernel::GetStateSpace(boost::property_tree::ptree sdeOptions) const
{
		auto periodicCast1 = std::dynamic_pointer_cast<PeriodicKernel>(kernel1);
		auto periodicCast2 = std::dynamic_pointer_cast<PeriodicKernel>(kernel2);

		if(periodicCast1 && (!periodicCast2)){
			return GetProductStateSpace(periodicCast1, kernel2, sdeOptions);
		}else if((!periodicCast1) && periodicCast2){
			return GetProductStateSpace(periodicCast2, kernel1, sdeOptions);
		}else{
      int status = 0;

      std::unique_ptr<char, void(*)(void*)> res1 {abi::__cxa_demangle(typeid(*kernel1).name(), NULL, NULL, &status), std::free};
      std::string type1 = res1.get();
      std::unique_ptr<char, void(*)(void*)> res2 {abi::__cxa_demangle(typeid(*kernel2).name(), NULL, NULL, &status), std::free};
      std::string type2 = res2.get();

      throw muq::NotImplementedError("ERROR in ProductKernel::GetStateSpace().  The GetStateSpace() function has not been implemented for these types: \"" + type1 + "\" and \"" + type2 + "\"");
		}
};



// See "Explicit Link Between Periodic
std::tuple<std::shared_ptr<muq::Modeling::LinearSDE>, std::shared_ptr<muq::Modeling::LinearOperator>, Eigen::MatrixXd> ProductKernel::GetProductStateSpace(std::shared_ptr<PeriodicKernel> const& kernel1,
                                                                                                                                             	              std::shared_ptr<KernelBase>     const& kernel2,
                                                                                                                                                            boost::property_tree::ptree sdeOptions) const
{

    auto periodicGP = kernel1->GetStateSpace(sdeOptions);
    auto periodicSDE = std::get<0>(periodicGP);

    auto periodicF = std::dynamic_pointer_cast<muq::Modeling::BlockDiagonalOperator>(periodicSDE->GetF());
    assert(periodicF);

    auto periodicL = std::dynamic_pointer_cast<muq::Modeling::BlockDiagonalOperator>(periodicSDE->GetL());
    assert(periodicL);

    auto otherGP = kernel2->GetStateSpace(sdeOptions);
    auto otherSDE = std::get<0>(otherGP);
    auto otherF = otherSDE->GetF();
    auto otherL = otherSDE->GetL();
    auto otherH = std::get<1>(otherGP);

    /// Construct the new F operator
    std::vector<std::shared_ptr<muq::Modeling::LinearOperator>> newBlocks( periodicF->GetBlocks().size() );
    for(int i=0; i<newBlocks.size(); ++i)
        newBlocks.at(i) = muq::Modeling::KroneckerSum(otherF, periodicF->GetBlock(i) );

    auto newF = std::make_shared<muq::Modeling::BlockDiagonalOperator>(newBlocks);

    /// Construct the new L operator
    for(int i=0; i<newBlocks.size(); ++i)
        newBlocks.at(i) = std::make_shared<muq::Modeling::KroneckerProductOperator>(otherL, periodicL->GetBlock(i) );

    auto newL = std::make_shared<muq::Modeling::BlockDiagonalOperator>(newBlocks);

    /// Construct the new H operator
    Eigen::MatrixXd Hblock(1,2);
    Hblock << 1.0, 0.0;

    for(int i=0; i<newBlocks.size(); ++i)
        newBlocks.at(i) = std::make_shared<muq::Modeling::KroneckerProductOperator>(otherH, muq::Modeling::LinearOperator::Create(Hblock) );

    auto newH = std::make_shared<muq::Modeling::BlockRowOperator>(newBlocks);

    // Construct Pinf
    Eigen::MatrixXd periodicP = std::get<2>(periodicGP);
    Eigen::MatrixXd otherP = std::get<2>(otherGP);

    Eigen::MatrixXd Pinf = Eigen::MatrixXd::Zero(periodicP.rows()*otherP.rows(), periodicP.cols()*otherP.cols());
    for(int i=0; i<newBlocks.size(); ++i)
        Pinf.block(2*i*otherP.rows(), 2*i*otherP.rows(), 2*otherP.rows(), 2*otherP.cols()) = muq::Modeling::KroneckerProduct(otherP, periodicP.block(2*i,2*i,2,2));

    // Construct Q
    Eigen::MatrixXd const& otherQ = otherSDE->GetQ();
    Eigen::MatrixXd Q = Eigen::MatrixXd::Zero(otherQ.rows()*periodicP.rows(), otherQ.cols()*periodicP.cols());
    for(int i=0; i<newBlocks.size(); ++i)
        Q.block(2*i*otherQ.rows(), 2*i*otherQ.rows(), 2*otherQ.rows(), 2*otherQ.cols()) = muq::Modeling::KroneckerProduct(otherQ, periodicP.block(2*i,2*i,2,2));

    // Construct the new statespace GP
    auto newSDE = std::make_shared<muq::Modeling::LinearSDE>(newF, newL, Q, sdeOptions);
    return std::make_tuple(newSDE, newH, Pinf);
}


std::shared_ptr<ProductKernel> muq::Approximation::operator*(std::shared_ptr<KernelBase> k1, std::shared_ptr<KernelBase> k2)
{
  return std::make_shared<ProductKernel>(k1,k2);
}
