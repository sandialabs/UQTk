#include <gtest/gtest.h>

#include "MUQ/Approximation/Quadrature/AdaptiveSmolyakQuadrature.h"
#include "MUQ/Approximation/Quadrature/SmolyakQuadrature.h"
#include "MUQ/Approximation/Quadrature/GaussQuadrature.h"
#include "MUQ/Approximation/Quadrature/ClenshawCurtisQuadrature.h"
#include "MUQ/Approximation/Quadrature/GaussPattersonQuadrature.h"
#include "MUQ/Approximation/Quadrature/ExponentialGrowthQuadrature.h"

#include "MUQ/Approximation/Polynomials/Legendre.h"

#include "MUQ/Modeling/CwiseOperators/CwiseUnaryOperator.h"
#include "MUQ/Utilities/MultiIndices/MultiIndexFactory.h"

#include <Eigen/Core>

using namespace muq::Approximation;
using namespace muq::Utilities;
using namespace muq::Modeling;


TEST(Quadrature, AdaptiveSmolyak_StaticGaussQuad) {

  unsigned int dim = 2;

  // Define the model
  auto model = std::make_shared<CosOperator>(dim);

  // Define the quadrature rules
  auto poly1d = std::make_shared<Legendre>();
  auto quad1d = std::make_shared<GaussQuadrature>(poly1d);

  AdaptiveSmolyakQuadrature smolyQuad(model, {quad1d, quad1d});

  unsigned int maxOrder = 7;
  auto multiSet = MultiIndexFactory::CreateTotalOrder(dim, maxOrder);

  Eigen::VectorXd res = smolyQuad.Compute(multiSet);

  EXPECT_NEAR(4.0*std::sin(1.0), res(0), 1e-10);
  EXPECT_NEAR(4.0*std::sin(1.0), res(1), 1e-10);

  EXPECT_LT(smolyQuad.Error(), 1e-12);
}


TEST(Quadrature, AdaptiveSmolyak_AdaptiveGaussQuad) {

  unsigned int dim = 2;

  // Define the model
  auto model = std::make_shared<CosOperator>(dim);

  // Define the quadrature rules
  auto poly1d = std::make_shared<Legendre>();
  auto quad1d = std::make_shared<GaussQuadrature>(poly1d);

  boost::property_tree::ptree options;
  options.put("ShouldAdapt",true);
  options.put("ErrorTol",1e-12);

  AdaptiveSmolyakQuadrature smolyQuad(model, {quad1d, quad1d});

  unsigned int maxOrder = 1;
  auto multiSet = MultiIndexFactory::CreateTotalOrder(dim, maxOrder);

  Eigen::VectorXd res = smolyQuad.Compute(multiSet, options);

  EXPECT_NEAR(4.0*std::sin(1.0), res(0), 1e-10);
  EXPECT_NEAR(4.0*std::sin(1.0), res(1), 1e-10);

  EXPECT_LT(smolyQuad.Error(), 1e-12);
}

TEST(Quadrature, AdaptiveSmolyak_AdaptiveClenshawCurtis) {

  unsigned int dim = 2;

  // Define the model
  auto model = std::make_shared<CosOperator>(dim);

  // Define the quadrature rules
  auto quad1d = std::make_shared<ClenshawCurtisQuadrature>();

  boost::property_tree::ptree options;
  options.put("ShouldAdapt",true);
  options.put("ErrorTol",1e-12);

  AdaptiveSmolyakQuadrature smolyQuad(model, {quad1d, quad1d});

  unsigned int maxOrder = 1;
  auto multiSet = MultiIndexFactory::CreateTotalOrder(dim, maxOrder);

  Eigen::VectorXd res = smolyQuad.Compute(multiSet, options);

  EXPECT_NEAR(4.0*std::sin(1.0), res(0), 1e-10);
  EXPECT_NEAR(4.0*std::sin(1.0), res(1), 1e-10);

  EXPECT_LT(smolyQuad.Error(), 1e-12);
}


TEST(Quadrature, AdaptiveSmolyak_AdaptiveGaussPatterson) {

  unsigned int dim = 2;

  // Define the model
  auto model = std::make_shared<CosOperator>(dim);

  // Define the quadrature rules
  auto quad1d = std::make_shared<GaussPattersonQuadrature>();

  boost::property_tree::ptree options;
  options.put("ShouldAdapt",true);
  options.put("ErrorTol",1e-12);

  AdaptiveSmolyakQuadrature smolyQuad(model, {quad1d, quad1d});

  unsigned int maxOrder = 1;
  auto multiSet = MultiIndexFactory::CreateTotalOrder(dim, maxOrder);

  Eigen::VectorXd res = smolyQuad.Compute(multiSet, options);

  EXPECT_NEAR(4.0*std::sin(1.0), res(0), 1e-10);
  EXPECT_NEAR(4.0*std::sin(1.0), res(1), 1e-10);

  EXPECT_LT(smolyQuad.Error(), 1e-12);
}

TEST(Quadrature, AdaptiveSmolyak_AdaptiveExponential) {

  unsigned int dim = 2;

  // Define the model
  auto model = std::make_shared<CosOperator>(dim);

  // Define the quadrature rules
  auto poly1d = std::make_shared<Legendre>();
  auto baseQuad = std::make_shared<GaussQuadrature>(poly1d);
  auto quad1d = std::make_shared<ExponentialGrowthQuadrature>(baseQuad);

  boost::property_tree::ptree options;
  options.put("ShouldAdapt",true);
  options.put("ErrorTol",1e-12);

  AdaptiveSmolyakQuadrature smolyQuad(model, {quad1d, quad1d});

  unsigned int maxOrder = 1;
  auto multiSet = MultiIndexFactory::CreateTotalOrder(dim, maxOrder);

  Eigen::VectorXd res = smolyQuad.Compute(multiSet, options);

  EXPECT_NEAR(4.0*std::sin(1.0), res(0), 1e-10);
  EXPECT_NEAR(4.0*std::sin(1.0), res(1), 1e-10);

  EXPECT_LT(smolyQuad.Error(), 1e-12);
}
