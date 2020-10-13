#ifndef KARHUNENLOEVEEXPANSION_H
#define KARHUNENLOEVEEXPANSION_H

#include <Eigen/Core>

#include "MUQ/Approximation/GaussianProcesses/KernelBase.h"
#include "MUQ/Approximation/GaussianProcesses/KarhunenLoeveBase.h"

#include <boost/property_tree/ptree.hpp>

namespace muq
{
namespace Approximation
{

    /** @class KarhunenLoeveExpansion
        @ingroup GaussianProcesses
        @brief Used to compute and evaluate the Karhunen-Loeve decomposition of a zero mean Gaussian process.
        @sealso SeperableKarhunenLoeve
        @details

## Background ##

According to the Karhunen-Loeve theorem, under mild technical conditions, a zero mean (i.e., centered) stochastic process admits a decomposition of the form
\f[
u(x) = \sum_{k=1}^\infty \phi_k(x) z_k,
\f]
where \f$u(x)\f$ is the original stochastic process, \f$\phi_k(x)\f$ is deterministic basis function, and \f$z_k\f$ are uncorrelated random variables.  Notice that this Karhunen-Loeve decomposition is a continuous analog of standard principle component decompositions for finite dimensional random variables.

 For a Gaussian process, the basis functions \f$\phi_k(x)\f$ can be calculated from the covariance kernel defining the process.  In particular, each \f$z_k\f$ will be an indepedent standard normal random variable and \f$\phi_k(x)\f$ will be given by
\f[
\phi_k(x) = \sqrt{\lambda_k} e_k(x),
\f]
where \f$\lambda_k\f$ and \f$e_k(x)\f$ are the \f$k^{th}\f$ eigenvalue and eigenvector of the covariance kernel \f$k(x,x^\prime)\f$.  In particular, \f$\lambda_k\f$ and \f$e_k(x)\f$ satisfy
\f[
\int_\Omega k(x,x^\prime) e_k(x) dx = \lambda_k e_k(x^\prime).
\f]
The constructor of this class uses the Nystrom method to approximately solve for \f$\lambda_k\f$ and \f$e_k(x)\f$.  The method starts with a set of \f$N\f$ seed (quadrature) points \f$x_i\f$ and a set of weights \f$w_i\f$ that satisfy
\f[
\int_\Omega k(x,x^\prime) e_k(x) dx \approx \sum_{i=1}^N w_i k(x_i,x^\prime) e_k(x^\prime),
\f]
which allows us to approximately solve for \f$e_k(x)\f$ and \f$\lambda_k\f$ by replacing \f$x^\prime\f$ with $\f$x_i\f$ for each \f$i\f$ and solving the resulting \f$N\times N\f$ matrix eigenproblem.   Once the matrix problem is solved, we will have approximations of \f$e_k(x_i)\f$.  The value of the eigenfunction at other locations is approximated using the covariance kernel
\f[
e_k(x^\prime) = \frac{1}{\lambda_k} \sum_{i=1}^N w_i k(x_i, x^\prime) e_k(x_i).
\f]
Notice that the approximation quality largely depends on how well the seed points and weights approximate the integral.  For simple domains (e.g., rectangular), it is thus useful to choose the seed points and weights from efficient quadrature rules.

## Typical Usage ##
\code{.cpp}
#include "MUQ/Approximation/GaussianProcesses/CovarianceKernels.h"
#include "MUQ/Approximation/GaussianProcesses/KarhunenLoeveExpansion.h"

using namespace muq::Approximation;

// other setup...

int main(){

    // first, create the covariance kernel
    double sigma2 = 1.0;
    double length = 0.2;
    double nu = 3.0/2.0;

    auto kernel = std::make_shared<MaternKernel>(1, sigma2, length, nu);


    // Define some seed points and weights -- quadrature points and weights would be a great choice
    const int numSeeds = 10;

    Eigen::MatrixXd seedPts(1,numSeeds);

    for(int i=0; i<numSeeds; ++i)
        seedPts(0, i) = static_cast<double>(i)/static_cast<double>(numSeeds);

    Eigen::VectorXd seedWts = (ub-lb)/double(numSeeds)*Eigen::VectorXd::Ones(numSeeds);


    // Construct the Karhunen-Loeve decomposition based on a discrete eigenvalue problem at the seed points
    KarhunenLoeveExpansion kl(kernel, seedPts, seedWts);


    // Evaluate the KL modes at a bunch of other points
    const int numEval = 200;

    Eigen::MatrixXd evalPts(1,numEval);

    for(int i=0; i<numEval; ++i)
        evalPts(0, i) = static_cast<double>(i)/static_cast<double>(numEval);


    // Each column of the modes matrix contains a basis function from the KL expansion
    Eigen::MatrixXd modes = kl.GetModes(evalPts);

    return 0;
}
\endcode
    */
    class KarhunenLoeveExpansion : public KarhunenLoeveBase
    {

    public:

        KarhunenLoeveExpansion(std::shared_ptr<KernelBase> kernelIn,
                               Eigen::MatrixXd      const& seedPtsIn,
                               Eigen::VectorXd      const& seedWtsIn,
                               boost::property_tree::ptree options = boost::property_tree::ptree());

        /** Evaluates the KL modes at one or more locations.  Each column of the pts matrix contains a point where we want to evaluate the modes.  Each column of the output contains a mode.  Each row of the output corresponds to an input point.
        */
        virtual Eigen::MatrixXd GetModes(Eigen::Ref<const Eigen::MatrixXd> const& pts) const override;

        virtual unsigned int NumModes() const override;

    private:


        // Points used to discretize the KL modes.
        Eigen::MatrixXd seedPts;
        Eigen::VectorXd seedWts;

        // The covariance kernel used to construct this expansion
        std::shared_ptr<KernelBase> covKernel;

        // Values of the KL modes at the seed points.  Each column corresponds to a basis function and each row to a pt function
        Eigen::MatrixXd modeVecs;
        Eigen::VectorXd modeEigs;

    }; // class KarhuneLoeveExpansion



}// namespace Approximation
}// namespace muq


#endif
