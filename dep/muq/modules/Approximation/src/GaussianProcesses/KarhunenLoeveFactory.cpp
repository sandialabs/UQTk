#include "MUQ/Approximation/GaussianProcesses/KarhunenLoeveFactory.h"

#include "MUQ/Approximation/GaussianProcesses/ProductKernel.h"

using namespace boost::property_tree;
using namespace muq::Approximation;



KarhunenLoeveFactory::KarhunenLoeveFactory(std::shared_ptr<KernelBase> kernelIn) : kernel(kernelIn)
{
    // Pull out the separable parts of the kernel
    
};


void KarhunenLoeveFactory::Compute(Eigen::VectorXd const& p0,
                                   Eigen::VectorXd const& p1,
                                   Eigen::VectorXi const& ns,
                                   boost::property_tree::ptree options)
{

    
    


}


