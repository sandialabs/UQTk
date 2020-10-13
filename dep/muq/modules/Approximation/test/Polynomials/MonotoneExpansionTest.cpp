#include "MUQ/Approximation/Polynomials/BasisExpansion.h"
#include "MUQ/Approximation/Polynomials/MonotoneExpansion.h"

#include "MUQ/Approximation/Polynomials/Monomial.h"
#include "MUQ/Utilities/MultiIndices/MultiIndexFactory.h"

#include "gtest/gtest.h"

using namespace muq::Approximation;
using namespace muq::Utilities;


class Approximation_MonotoneExpansion1d : public::testing::Test {
public:
  Approximation_MonotoneExpansion1d() {

    auto monomial = std::make_shared<Monomial>();
    auto bases = std::vector<std::shared_ptr<IndexedScalarBasis>>(1, monomial);

    std::shared_ptr<MultiIndexSet> multis = MultiIndexFactory::CreateTotalOrder(1, 2);

    monoPart = std::make_shared<BasisExpansion>(bases, multis);
    expansion1 = std::make_shared<MonotoneExpansion>(monoPart);
    expansion2 = std::make_shared<MonotoneExpansion>(monoPart, true);

    coeffs.resize(4);
    coeffs << 0.0, -0.01, 0.1, -0.5; // intercept, slope

    expansion1->SetCoeffs(coeffs);
    expansion2->SetCoeffs(coeffs);
  }

  virtual ~Approximation_MonotoneExpansion1d() {};

  std::shared_ptr<BasisExpansion> monoPart;
  std::shared_ptr<MonotoneExpansion> expansion1, expansion2;
  Eigen::VectorXd coeffs;

};

TEST_F(Approximation_MonotoneExpansion1d, Evaluate){

  Eigen::VectorXd evalPt(1);

  int numSteps = 10;
  double ub = 2.0;
  double lb = -2.0;
  double dx = (ub-lb)/numSteps;

  evalPt << lb;

  double oldOutput = expansion1->Evaluate(evalPt).at(0)(0);

  for(int i=1; i<numSteps; ++i){
      evalPt(0) = lb + dx*double(i);
      double newOutput = expansion1->Evaluate(evalPt).at(0)(0);
      EXPECT_GT(newOutput, oldOutput);

      double monoEval = monoPart->Evaluate(evalPt).at(0)(0);
      double trueMono = coeffs(1) + coeffs(2)*evalPt(0) + coeffs(3)*evalPt(0)*evalPt(0);
      EXPECT_DOUBLE_EQ(trueMono, monoEval);

      double truth = coeffs(1)*coeffs(1)*evalPt(0)
                     + coeffs(1)*coeffs(2)*std::pow(evalPt(0),2.0)
                     + (1.0/3.0)*(2.0*coeffs(1)*coeffs(3)+coeffs(2)*coeffs(2))*std::pow(evalPt(0),3.0)
                     + 0.5*coeffs(2)*coeffs(3)*std::pow(evalPt(0),4.0)
                     + 0.2*coeffs(3)*coeffs(3)*std::pow(evalPt(0),5.0);
      EXPECT_NEAR(truth, newOutput, 1e-4);


      Eigen::MatrixXd jac = expansion1->Jacobian(0,0,evalPt);
      EXPECT_GT(jac(0,0), 0.0);

      oldOutput = newOutput;
  }
}


TEST_F(Approximation_MonotoneExpansion1d, JacobianFD){

  Eigen::VectorXd evalPt(1);

  int numSteps = 10;
  double ub = 2.0;
  double lb = -2.0;
  double dx = (ub-lb)/numSteps;

  evalPt << lb;

  const double eps = 1e-3;

  for(int i=1; i<numSteps; ++i){

      evalPt(0) = lb + dx*double(i);
      double f1 = expansion1->Evaluate(evalPt).at(0)(0);

      Eigen::MatrixXd jac = expansion1->Jacobian(0,0,evalPt);

      evalPt(0) += eps;
      double f2 = expansion1->Evaluate(evalPt).at(0)(0);

      EXPECT_NEAR((f2-f1)/eps, jac(0,0), 3e-3);
  }
}

TEST_F(Approximation_MonotoneExpansion1d, CoeffJacobian){

  Eigen::VectorXd evalPt(1);
  evalPt << 0.25;

  const double eps = 1e-3;

  double f1 = expansion1->Evaluate(evalPt).at(0)(0);
  Eigen::MatrixXd jac = expansion2->Jacobian(1, 0, evalPt, coeffs);

  for(int i=0; i<coeffs.size(); ++i){
    Eigen::VectorXd newCoeffs = coeffs;
    newCoeffs(i) += eps;

    double f2 = expansion2->Evaluate(evalPt, newCoeffs).at(0)(0);

    EXPECT_NEAR((f2-f1)/eps, jac(0,i), 1e-3);
  }
}

TEST_F(Approximation_MonotoneExpansion1d, Determinant){

  Eigen::VectorXd evalPt(1);
  evalPt << 0.25;

  double logDet1 = expansion1->LogDeterminant(evalPt);
  Eigen::VectorXd grad = expansion2->GradLogDeterminant(evalPt,coeffs);

  const double eps = 1e-7;

  for(int i=0; i<coeffs.size(); ++i){
    Eigen::VectorXd newCoeffs = coeffs;
    newCoeffs(i) += eps;
    double logDet2 = expansion2->LogDeterminant(evalPt,newCoeffs);

    EXPECT_NEAR((logDet2-logDet1)/eps, grad(i), 1e-3);
  }
}


TEST_F(Approximation_MonotoneExpansion1d, EvaluateWithCoeffs){

  Eigen::VectorXd evalPt(1);

  Eigen::VectorXd newCoeffs(4);
  newCoeffs << 0.1, 0.5, 2.0, -0.1;
  coeffs = newCoeffs;

  int numSteps = 10;
  double ub = 2.0;
  double lb = -2.0;
  double dx = (ub-lb)/numSteps;

  evalPt << lb;
  double oldOutput = expansion2->Evaluate(evalPt, newCoeffs).at(0)(0);

  for(int i=1; i<numSteps; ++i){
      evalPt(0) = lb + dx*double(i);
      double newOutput = expansion2->Evaluate(evalPt, newCoeffs).at(0)(0);
      EXPECT_GT(newOutput, oldOutput);

      double truth = coeffs(0) + coeffs(1)*coeffs(1)*evalPt(0)
                     + coeffs(1)*coeffs(2)*std::pow(evalPt(0),2.0)
                     + (1.0/3.0)*(2.0*coeffs(1)*coeffs(3)+coeffs(2)*coeffs(2))*std::pow(evalPt(0),3.0)
                     + 0.5*coeffs(2)*coeffs(3)*std::pow(evalPt(0),4.0)
                     + 0.2*coeffs(3)*coeffs(3)*std::pow(evalPt(0),5.0);
      EXPECT_NEAR(truth, newOutput, 5e-2);

      oldOutput = newOutput;
  }
}

class Approximation_MonotoneExpansion2d : public::testing::Test {
public:
  Approximation_MonotoneExpansion2d() {

    auto monomial = std::make_shared<Monomial>();

    // Build the general pieces
    auto bases = std::vector<std::shared_ptr<IndexedScalarBasis>>(2, monomial);
    std::shared_ptr<MultiIndexSet> constantMulti = MultiIndexFactory::CreateTotalOrder(2, 0);
    std::shared_ptr<MultiIndexSet> multis = MultiIndexFactory::CreateTotalOrder(2, 2, 0, std::make_shared<DimensionLimiter>(0,1));
    Eigen::MatrixXd coeffs(1,3);
    coeffs << 0.0, 0.5, -0.1;
    generalParts.push_back(std::make_shared<BasisExpansion>(bases, constantMulti));
    generalParts.push_back(std::make_shared<BasisExpansion>(bases, multis, coeffs));

    // Build the monotone pieces
    coeffs << 0.1, 1.0, 0.5;
    monoParts.push_back(std::make_shared<BasisExpansion>(bases, multis, coeffs));

    multis = MultiIndexFactory::CreateTotalOrder(2, 2);
    coeffs = Eigen::MatrixXd::Random(1,multis->Size());
    monoParts.push_back(std::make_shared<BasisExpansion>(bases, multis, coeffs));

    // Construct the expansion
    expansion1 = std::make_shared<MonotoneExpansion>(generalParts, monoParts);
    expansion2 = std::make_shared<MonotoneExpansion>(generalParts, monoParts, true);
  };

  virtual ~Approximation_MonotoneExpansion2d() {};

  std::vector<std::shared_ptr<BasisExpansion>> generalParts;
  std::vector<std::shared_ptr<BasisExpansion>> monoParts;

  std::shared_ptr<MonotoneExpansion> expansion1, expansion2;
};

TEST_F(Approximation_MonotoneExpansion2d, Evaluate){

  Eigen::VectorXd evalPt = Eigen::VectorXd::Zero(2);

  int numSteps = 10;
  double ub = 2.0;
  double lb = -2.0;
  double dx = (ub-lb)/numSteps;

  evalPt(1) = lb;
  Eigen::VectorXd oldOutput = expansion1->Evaluate(evalPt).at(0);

  for(int i=1; i<numSteps; ++i){
    evalPt(1) = lb + dx*double(i);
    Eigen::VectorXd newOutput = expansion1->Evaluate(evalPt).at(0);
    EXPECT_GT(newOutput(1), oldOutput(1));

    oldOutput = newOutput;
  }

}

TEST_F(Approximation_MonotoneExpansion2d, Inverse){

  Eigen::VectorXd tgtPt = Eigen::Vector2d{0.1, 0.5};
  Eigen::VectorXd refPt = expansion1->EvaluateForward(tgtPt);
  Eigen::VectorXd invPt = expansion1->EvaluateInverse(refPt);

  EXPECT_NEAR(tgtPt(0), invPt(0), 1e-8);
  EXPECT_NEAR(tgtPt(1), invPt(1), 1e-8);

  tgtPt = Eigen::Vector2d{-9.2, 5.5};
  refPt = expansion1->EvaluateForward(tgtPt);
  invPt = expansion1->EvaluateInverse(refPt);

  EXPECT_NEAR(tgtPt(0), invPt(0), 1e-8);
  EXPECT_NEAR(tgtPt(1), invPt(1), 1e-8);
}



TEST_F(Approximation_MonotoneExpansion2d, Jacobian){

  Eigen::VectorXd evalPt(2);
  evalPt << 0.1, 0.2;

  Eigen::VectorXd f1 = expansion1->Evaluate(evalPt).at(0);
  Eigen::MatrixXd jac = expansion1->Jacobian(0,0, evalPt);

  EXPECT_GT(jac(0,0),0.0);
  EXPECT_GT(jac(1,1),0.0);

  const double eps = 1e-6;
  Eigen::VectorXd newPt = evalPt;
  newPt(0) += eps;
  Eigen::VectorXd f2 = expansion1->Evaluate(newPt).at(0);

  for(int i=0; i<f2.size(); ++i)
    EXPECT_NEAR((f2(i)-f1(i))/eps, jac(i,0), 1e-3);

  newPt = evalPt;
  newPt(1) += eps;
  f2 = expansion1->Evaluate(newPt).at(0);

  for(int i=0; i<f2.size(); ++i)
    EXPECT_NEAR((f2(i)-f1(i))/eps, jac(i,1), 1e-3);
}

TEST_F(Approximation_MonotoneExpansion2d, CoeffsJacobian){

  Eigen::VectorXd evalPt(2);
  evalPt << 0.1, 0.2;

  Eigen::VectorXd coeffs = expansion2->GetCoeffs();

  const double eps = 1e-3;

  Eigen::VectorXd f1 = expansion2->Evaluate(evalPt, coeffs).at(0);
  Eigen::MatrixXd jac = expansion2->Jacobian(1,0,evalPt,coeffs);

  for(int i=0; i<coeffs.size(); ++i){
    Eigen::VectorXd newCoeffs = coeffs;
    newCoeffs(i) += eps;

    Eigen::VectorXd f2 = expansion2->Evaluate(evalPt, newCoeffs).at(0);

    for(int j=0; j<f2.size(); ++j)
      EXPECT_NEAR((f2(j)-f1(j))/eps, jac(j,i), 1e-3);
  }
}

TEST_F(Approximation_MonotoneExpansion2d, Determinant){

  Eigen::VectorXd evalPt(2);
  evalPt << 0.1, 0.2;

  Eigen::VectorXd coeffs = expansion1->GetCoeffs();

  double logDet1 = expansion1->LogDeterminant(evalPt);

  Eigen::VectorXd grad = expansion1->GradLogDeterminant(evalPt);

  const double eps = 1e-7;

  for(int i=0; i<coeffs.size(); ++i){
    Eigen::VectorXd newCoeffs = coeffs;
    newCoeffs(i) += eps;
    double logDet2 = expansion1->LogDeterminant(evalPt,newCoeffs);

    EXPECT_NEAR((logDet2-logDet1)/eps, grad(i), 1e-3);
  }

}
