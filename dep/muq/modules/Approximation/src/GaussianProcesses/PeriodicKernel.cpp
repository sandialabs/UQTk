#include "MUQ/Approximation/GaussianProcesses/PeriodicKernel.h"

#include "MUQ/Modeling/LinearAlgebra/LinearOperator.h"
#include "MUQ/Modeling/LinearAlgebra/EigenLinearOperator.h"
#include "MUQ/Modeling/LinearAlgebra/BlockDiagonalOperator.h"
#include "MUQ/Modeling/LinearAlgebra/IdentityOperator.h"

#include "MUQ/Approximation/GaussianProcesses/StateSpaceGP.h"


#include "MUQ/Modeling/LinearAlgebra/CompanionMatrix.h"
#include "MUQ/Modeling/LinearAlgebra/LyapunovSolver.h"

#include "MUQ/Modeling/LinearSDE.h"


#include <unsupported/Eigen/Polynomials>
#include <Eigen/SparseCore>

#include <boost/math/special_functions/bessel.hpp>
#include <boost/property_tree/ptree.hpp>

using namespace muq::Approximation;
using namespace muq::Modeling;


/** Implements the 2x2 block matrix in equation 28 of Solin and Sarkka */
class PeriodicKernel_F_block : public muq::Modeling::LinearOperator
{
public:
    PeriodicKernel_F_block(const double wjIn) : muq::Modeling::LinearOperator(2,2), wj(wjIn){};

    virtual Eigen::MatrixXd Apply(Eigen::Ref<const Eigen::MatrixXd> const& x) override
    {
        Eigen::MatrixXd output(2,x.cols());
        output.row(0) = -wj*x.row(1);
        output.row(1) = wj*x.row(0);
        return output;
    }

    virtual Eigen::MatrixXd ApplyTranspose(Eigen::Ref<const Eigen::MatrixXd> const& x) override
    {
        Eigen::MatrixXd output(2,x.cols());
        output.row(0) = wj*x.row(1);
        output.row(1) = -wj*x.row(0);
        return output;
    }

    virtual Eigen::MatrixXd GetMatrix() override
    {
        Eigen::MatrixXd output(2,2);
        output << 0.0, -wj,
                   wj, 0.0;

        return output;
    }

private:
    const double wj;

};

std::tuple<std::shared_ptr<muq::Modeling::LinearSDE>, std::shared_ptr<muq::Modeling::LinearOperator>, Eigen::MatrixXd> PeriodicKernel::GetStateSpace(boost::property_tree::ptree sdeOptions) const
{

  const double sigma2 = cachedParams(0);
  const double length = cachedParams(1);
  const double period = cachedParams(2);

    // This is the same as J in the paper
    const int numTerms = sdeOptions.get("PeriodicKernel.StateSpace.NumTerms",4);

    const double w0 = 2.0*pi/period;

    const double l2 = std::pow(length, 2.0);


    Eigen::VectorXd q2s(numTerms+1); // holds all of the q_j^2 values from eqn 27 in "Explicit Link Between Periodic Covariance Functions and State Space Models"

    q2s(0) = boost::math::cyl_bessel_i(0, 1.0/l2) / exp(1.0/l2);
    for(int i=1; i<numTerms+1; ++i)
        q2s(i) = 2.0 * boost::math::cyl_bessel_i(i, 1.0/l2) / exp(1.0/l2);


    // SET UP F
    std::vector<std::shared_ptr<LinearOperator>> fBlocks(numTerms+1);
    for(int i=0; i<numTerms+1; ++i)
        fBlocks.at(i) = std::make_shared<PeriodicKernel_F_block>(w0*i);

    auto F = std::make_shared<BlockDiagonalOperator>(fBlocks);

    // SET UP L
    std::vector<std::shared_ptr<LinearOperator>> lBlocks(numTerms+1);
    for(int i=0; i<numTerms+1; ++i)
        lBlocks.at(i) = std::make_shared<IdentityOperator>(2);

    auto L = std::make_shared<BlockDiagonalOperator>(lBlocks);

    // Set up Pinf
    Eigen::MatrixXd Pinf = Eigen::MatrixXd::Zero(2*(numTerms+1), 2*(numTerms+1));
    for(int i=0; i<numTerms+1; ++i)
        Pinf.block(2*i, 2*i, 2, 2) = q2s(i) * Eigen::MatrixXd::Identity(2,2);

    // Set up Q
    Eigen::MatrixXd Q = Eigen::MatrixXd::Zero(2*(numTerms+1), 2*(numTerms+1));

    // SET UP H
    std::vector<Eigen::Triplet<double>> Hcoeffs;
    for(int i=0; i<numTerms+1; ++i)
        Hcoeffs.push_back(Eigen::Triplet<double>(0, 2*i, 1.0));

    Eigen::SparseMatrix<double> Hmat(1, 2*(numTerms+1));
    Hmat.setFromTriplets(Hcoeffs.begin(), Hcoeffs.end());

    auto H = muq::Modeling::LinearOperator::Create(Hmat);

    // Create the SDE
    auto sde = std::make_shared<muq::Modeling::LinearSDE>(F,L,Q,sdeOptions);

    return std::make_tuple(sde, H, Pinf);

}




PeriodicKernel::PeriodicKernel(unsigned dim,
                               std::vector<unsigned> dimInds,
                               const double sigma2In,
                               const double lengthIn,
                               const double periodIn,
                               const Eigen::Vector2d sigmaBounds,
                               const Eigen::Vector2d lengthBounds,
                               const Eigen::Vector2d periodBounds) : KernelImpl<PeriodicKernel>(dim, dimInds, 1 , 3)
{
 SetupBounds(sigmaBounds, lengthBounds, periodBounds);

 cachedParams.resize(3);
 cachedParams(0) = sigma2In;
 cachedParams(1) = lengthIn;
 cachedParams(2) = periodIn;
};


PeriodicKernel::PeriodicKernel(unsigned dim,
                               const double sigma2In,
                               const double lengthIn,
                               const double periodIn,
                               const Eigen::Vector2d sigmaBounds,
                               const Eigen::Vector2d lengthBounds,
                               const Eigen::Vector2d periodBounds) : KernelImpl<PeriodicKernel>(dim, 1 , 3)
{
  SetupBounds(sigmaBounds, lengthBounds, periodBounds);

  cachedParams.resize(3);
  cachedParams(0) = sigma2In;
  cachedParams(1) = lengthIn;
  cachedParams(2) = periodIn;
};

void PeriodicKernel::SetupBounds(Eigen::Vector2d const& sigmaBounds,
                                 Eigen::Vector2d const& lengthBounds,
                                 Eigen::Vector2d const& periodBounds)
{
  paramBounds.resize(2,3);

  paramBounds(0,0) = sigmaBounds(0);
  paramBounds(1,0) = sigmaBounds(1);

  paramBounds(0,1) = lengthBounds(0);
  paramBounds(1,1) = lengthBounds(1);

  paramBounds(0,2) = periodBounds(0);
  paramBounds(1,2) = periodBounds(1);
};
