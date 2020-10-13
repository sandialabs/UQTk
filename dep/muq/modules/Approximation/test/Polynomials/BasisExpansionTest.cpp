#include "MUQ/Approximation/Polynomials/BasisExpansion.h"
#include "MUQ/Approximation/Polynomials/Legendre.h"
#include "MUQ/Approximation/Polynomials/Monomial.h"
#include "MUQ/Utilities/MultiIndices/MultiIndexFactory.h"

#include "gtest/gtest.h"

using namespace muq::Approximation;
using namespace muq::Utilities;



class Approximation_BasisExpansion : public::testing::Test {
public:
  Approximation_BasisExpansion() {

    // Define the bases used for each family.  Here, we use Legendre in all dims
    auto legendre = std::make_shared<Legendre>();
    bases = std::vector<std::shared_ptr<IndexedScalarBasis>>(numDims, legendre);
  }

  virtual ~Approximation_BasisExpansion() {};

  std::vector<std::shared_ptr<IndexedScalarBasis>> bases;
  const unsigned numDims = 2;

};

TEST_F(Approximation_BasisExpansion, Constructor1){

  BasisExpansion expansion(bases);

  // Evaluate with just the input
  Eigen::VectorXd evalPt = Eigen::VectorXd::Ones(numDims);
  Eigen::VectorXd output = expansion.Evaluate(evalPt).at(0);

  EXPECT_EQ(1, output.size());
  EXPECT_DOUBLE_EQ(0.0, output(0));
}

TEST_F(Approximation_BasisExpansion, Constructor2){

  std::shared_ptr<MultiIndexSet> multis = MultiIndexFactory::CreateTotalOrder(numDims, 2);
  EXPECT_EQ(6,multis->Size());

  BasisExpansion expansion(bases, multis);

  // Evaluate with just the input
  Eigen::VectorXd evalPt = Eigen::VectorXd::Ones(numDims);
  Eigen::VectorXd output = expansion.Evaluate(evalPt).at(0);

  EXPECT_EQ(1, output.size());
  EXPECT_DOUBLE_EQ(0.0, output(0));
}

TEST_F(Approximation_BasisExpansion, Constructor3){

  std::shared_ptr<MultiIndexSet> multis = MultiIndexFactory::CreateTotalOrder(numDims, 2);
  EXPECT_EQ(6,multis->Size());

  const unsigned numOut = 10;
  Eigen::MatrixXd coeffs = Eigen::MatrixXd::Zero(numOut, multis->Size());

  BasisExpansion expansion1(bases, multis, coeffs);

  // Evaluate with just the input
  Eigen::VectorXd evalPt = Eigen::VectorXd::Ones(numDims);
  Eigen::VectorXd output = expansion1.Evaluate(evalPt).at(0);

  EXPECT_EQ(numOut, output.size());
  for(int i=0; i<numOut; ++i)
    EXPECT_DOUBLE_EQ(0.0, output(i));

  BasisExpansion expansion2(bases, multis, coeffs, true);

  coeffs = Eigen::MatrixXd::Ones(numOut, multis->Size());
  Eigen::VectorXd coeffVec = Eigen::Map<Eigen::VectorXd>(coeffs.data(), coeffs.rows()*coeffs.cols());
  output = expansion2.Evaluate(evalPt, coeffVec).at(0);

  EXPECT_EQ(numOut, output.size());
  for(int i=0; i<numOut; ++i)
    EXPECT_DOUBLE_EQ(6.0, output(i));
}

TEST_F(Approximation_BasisExpansion, Linear){

  auto monomial = std::make_shared<Monomial>();
  bases = std::vector<std::shared_ptr<IndexedScalarBasis>>(1, monomial);

  std::shared_ptr<MultiIndexSet> multis = MultiIndexFactory::CreateTotalOrder(1, 1);
  EXPECT_EQ(2,multis->Size());

  Eigen::MatrixXd coeffs(1,2);
  coeffs << 1.0, 2.0; // intercept, slope

  BasisExpansion expansion(bases, multis, coeffs);

  // Evaluate with just the input
  Eigen::VectorXd evalPt = Eigen::VectorXd::Random(1);
  Eigen::VectorXd output = boost::any_cast<Eigen::VectorXd>(expansion.Evaluate(evalPt)[0]);

  EXPECT_EQ(1, output.size());
  EXPECT_DOUBLE_EQ(coeffs(0,0) + coeffs(0,1)*evalPt(0), output(0));

  Eigen::MatrixXd jac = boost::any_cast<Eigen::MatrixXd>(expansion.Jacobian(0,0,evalPt));
  EXPECT_EQ(1,jac.rows());
  EXPECT_EQ(1, jac.cols());
  EXPECT_EQ(coeffs(1), jac(0,0));

  Eigen::MatrixXd hess = expansion.SecondDerivative(0, 0, 0, evalPt, coeffs);

  EXPECT_EQ(1, hess.rows());
  EXPECT_EQ(1, hess.cols());
  EXPECT_DOUBLE_EQ(0.0, hess(0,0));

}

TEST_F(Approximation_BasisExpansion, Quadratic){

  auto monomial = std::make_shared<Monomial>();
  bases = std::vector<std::shared_ptr<IndexedScalarBasis>>(1, monomial);

  std::shared_ptr<MultiIndexSet> multis = MultiIndexFactory::CreateTotalOrder(1, 2);
  EXPECT_EQ(3,multis->Size());

  Eigen::MatrixXd coeffs(1,3);
  coeffs << 1.0, 2.0, 0.5; // intercept, slope

  BasisExpansion expansion(bases, multis, coeffs);

  // Evaluate with just the input
  Eigen::VectorXd evalPt = Eigen::VectorXd::Random(1);
  Eigen::VectorXd output = boost::any_cast<Eigen::VectorXd>(expansion.Evaluate(evalPt)[0]);

  EXPECT_EQ(1, output.size());
  EXPECT_DOUBLE_EQ(coeffs(0,0) + coeffs(0,1)*evalPt(0) + coeffs(0,2)*evalPt(0)*evalPt(0), output(0));

  Eigen::MatrixXd jac = boost::any_cast<Eigen::MatrixXd>(expansion.Jacobian(0,0,evalPt));
  EXPECT_EQ(1,jac.rows());
  EXPECT_EQ(1, jac.cols());
  EXPECT_EQ(coeffs(1) + 2.0*coeffs(2)*evalPt(0), jac(0,0));

  Eigen::MatrixXd hess = expansion.SecondDerivative(0, 0, 0, evalPt, coeffs);

  EXPECT_EQ(1, hess.rows());
  EXPECT_EQ(1, hess.cols());
  EXPECT_DOUBLE_EQ(2.0*coeffs(0,2), hess(0,0));
}
