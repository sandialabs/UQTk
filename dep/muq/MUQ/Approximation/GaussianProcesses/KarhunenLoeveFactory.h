#ifndef KARHUNENLOEVEFACTORY_H
#define KARHUNENLOEVEFACTORY_H


#include "MUQ/Approximation/GaussianProcesses/KernelBase.h"

#include <boost/property_tree/ptree.hpp>

#include <vector>
#include <memory>

namespace muq
{
namespace Approximation
{


    class KarhunenLoeveFactory
    {
        
    public:

        template<typename KernelType>
        KarhunenLoeveFactory(KernelType &kernelIn) : KarhunenLoeveFactory(kernelIn.GetPtr()){}
        
        KarhunenLoeveFactory(std::shared_ptr<KernelBase> kernelIn);
        
        void Compute(Eigen::VectorXd const& p0,
                     Eigen::VectorXd const& p1,
                     Eigen::VectorXi const& ns,
                     boost::property_tree::ptree options = boost::property_tree::ptree());

        Eigen::MatrixXd GetModes() const;
        Eigen::VectorXd GetWeights() const;

    private:
        std::shared_ptr<KernelBase> kernel;

        // Contains parts of the covariance kernel that can be separated based on dimension. i.e., the product of kernels acting in different dimensions
        std::vector<std::shared_ptr<KernelBase>> kernelParts;
        
        Eigen::MatrixXd modes;
        Eigen::VectorXd weights;


        // Extract the seperable components of the kernel (requires kernel to be set and will fill in the contents of kernelParts)
        static std::vector<std::shared_ptr<KernelBase>> SeparateKernel(std::shared_ptr<KernelBase> kernel);

    };



} // namespace Approximation
} // namespace muq


#endif
