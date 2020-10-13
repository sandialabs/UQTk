#include "MUQ/Approximation/GaussianProcesses/MaternKernel.h"

#include "MUQ/Modeling/LinearAlgebra/LinearOperator.h"
#include "MUQ/Modeling/LinearAlgebra/EigenLinearOperator.h"

#include "MUQ/Approximation/GaussianProcesses/StateSpaceGP.h"


#include "MUQ/Modeling/LinearAlgebra/CompanionMatrix.h"
#include "MUQ/Modeling/LinearAlgebra/LyapunovSolver.h"

#include "MUQ/Modeling/LinearSDE.h"


#include <unsupported/Eigen/Polynomials>
#include <Eigen/SparseCore>

#include <boost/property_tree/ptree.hpp>

using namespace muq::Approximation;



MaternKernel::MaternKernel(unsigned              dimIn,
                           std::vector<unsigned> dimInds,
                           double                sigma2In,
                           double                lengthIn,
                           double                nuIn,
                           Eigen::Vector2d       sigmaBounds,
                           Eigen::Vector2d       lengthBounds) : KernelImpl<MaternKernel>(dimIn, dimInds, 1, 2),
                                                                 nu(nuIn),
                                                                 scale(boost::math::tgamma(nuIn+0.5)/boost::math::tgamma(2.0*nuIn)),
                                                                 weights(BuildWeights(nu-0.5))
{
    CheckNu();

    paramBounds.resize(2,2);
    paramBounds(0,0) = sigmaBounds(0);
    paramBounds(1,0) = sigmaBounds(1);

    paramBounds(0,1) = lengthBounds(0);
    paramBounds(1,1) = lengthBounds(1);

    cachedParams.resize(2);
    cachedParams(0) = sigma2In;
    cachedParams(1) = lengthIn;
};

MaternKernel::MaternKernel(unsigned        dimIn,
                           double          sigma2In,
                           double          lengthIn,
                           double          nuIn,
                           Eigen::Vector2d sigmaBounds,
                           Eigen::Vector2d lengthBounds) : KernelImpl<MaternKernel>(dimIn, 1, 2),
                                                           nu(nuIn),
                                                           scale(boost::math::tgamma(nuIn+0.5)/boost::math::tgamma(2.0*nuIn)),
                                                           weights(BuildWeights(nu-0.5))
{
    CheckNu();

    paramBounds.resize(2,2);
    paramBounds(0,0) = sigmaBounds(0);
    paramBounds(1,0) = sigmaBounds(1);

    paramBounds(0,1) = lengthBounds(0);
    paramBounds(1,1) = lengthBounds(1);

    cachedParams.resize(2);
    cachedParams(0) = sigma2In;
    cachedParams(1) = lengthIn;
};


void MaternKernel::CheckNu() const{

    if(nu<=0)
        throw std::invalid_argument("The value of nu must be greater than 0.");

    if((std::round(nu-0.5)-(nu-0.5)) > 4.0*std::numeric_limits<double>::epsilon())
        throw std::invalid_argument("The value of nu must take the form nu=i-0.5 for a positive integer i.");

};


std::tuple<std::shared_ptr<muq::Modeling::LinearSDE>, std::shared_ptr<muq::Modeling::LinearOperator>, Eigen::MatrixXd> MaternKernel::GetStateSpace(boost::property_tree::ptree sdeOptions) const
{
    double sigma2 = cachedParams(0);
    double length = cachedParams(1);
    int p = nu-0.5;

    double lambda = sqrt(2.0*nu)/length;

    double q = 2.0*sigma2*boost::math::constants::root_pi<double>() * std::pow(lambda, 2*p+1) * tgamma(p+1) / tgamma(p+0.5);

    Eigen::VectorXd roots = -lambda*Eigen::VectorXd::Ones(p+1);
    Eigen::VectorXd poly;

    Eigen::roots_to_monicPolynomial( roots, poly );

    poly = poly.head(poly.size()-1).eval();

    auto F = std::make_shared<muq::Modeling::CompanionMatrix>(-1.0*poly);

    std::vector<Eigen::Triplet<double>> Lcoeffs;
    Lcoeffs.push_back(Eigen::Triplet<double>(poly.size()-1,0,1.0));

    Eigen::SparseMatrix<double> Lmat(poly.size(), 1);
    Lmat.setFromTriplets(Lcoeffs.begin(), Lcoeffs.end());

    auto L = muq::Modeling::LinearOperator::Create(Lmat);

    Eigen::MatrixXd Q(1,1);
    Q(0,0) = q;


    // Set up the stochastic differential equation
    auto sde = std::make_shared<muq::Modeling::LinearSDE>(F, L, Q, sdeOptions);


    // Define the observation operator, which is just (1,0,...,0) in this case
    std::vector<Eigen::Triplet<double>> Hcoeffs;
    Hcoeffs.push_back(Eigen::Triplet<double>(0,0,1.0));

    Eigen::SparseMatrix<double> Hmat(1,poly.size());
    Hmat.setFromTriplets(Hcoeffs.begin(), Hcoeffs.end());

    auto H = muq::Modeling::LinearOperator::Create(Hmat);

    // Solve the continuous time Lyapunov equation to find the stationary covariance
    Q = L->Apply(L->Apply(q*Eigen::VectorXd::Ones(1)).transpose());

    Eigen::MatrixXd Pinf = muq::Modeling::LyapunovSolver<double>().compute(F->GetMatrix().transpose(), Q).matrixX().real();

    return std::make_tuple(sde, H, Pinf);
}


Eigen::VectorXd MaternKernel::BuildWeights(int p)
{
  Eigen::VectorXd output(p+1);
  for(int i=0; i<=p; ++i){

    if(i<p){
      output(i) = static_cast<double>(Factorial(p+i))/(static_cast<double>(Factorial(i)*Factorial(p-i)));
    }else{
      output(i) = static_cast<double>(Factorial(p+i)) / static_cast<double>(Factorial(i));
    }
  }
  return output;
}
