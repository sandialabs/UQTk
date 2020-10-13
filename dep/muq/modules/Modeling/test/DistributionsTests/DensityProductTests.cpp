#include "MUQ/Modeling/Distributions/DensityProduct.h"
#include "MUQ/Modeling/Distributions/Gaussian.h"
#include "MUQ/Modeling/Distributions/UniformBox.h"
#include "MUQ/Modeling/WorkGraph.h"
#include "MUQ/Modeling/WorkGraphPiece.h"

#include "MUQ/Utilities/AnyHelpers.h"

#include <gtest/gtest.h>

using namespace muq::Modeling;
using namespace muq::Utilities;

TEST(Distribution_DensityProduct, EvaluateDensity) {
  std::shared_ptr<DensityBase> prod = std::make_shared<DensityProduct>(2);

  Eigen::VectorXd dens1(1);
  dens1 << 1.0;

  Eigen::VectorXd dens2(1);
  dens2 << -2.0;

  double res = prod->LogDensity(dens1,dens2);
  EXPECT_DOUBLE_EQ(-1.0, res);
}


TEST(Distribution_DensityProduct, EvaluateOnGraph) {
  auto graph = std::make_shared<WorkGraph>();

  Eigen::VectorXd mu(1);
  mu << 0.0;

  Eigen::MatrixXd bnds(1,2);
  bnds << 1.0, 2.0;

  std::shared_ptr<DensityBase> gaussDens = std::make_shared<Gaussian>(mu)->AsDensity();
  std::shared_ptr<DensityBase> uniformDens = std::make_shared<UniformBox>(bnds)->AsDensity();

  graph->AddNode(gaussDens, "Gaussian");
  graph->AddNode(uniformDens, "Uniform");

  auto prodPiece = std::make_shared<DensityProduct>(2);

  graph->AddNode(prodPiece, "Product");
  graph->AddEdge("Gaussian", 0, "Product", 0);
  graph->AddEdge("Uniform", 0, "Product", 1);

  auto graphPiece = graph->CreateWorkPiece("Product");

  Eigen::VectorXd x1(1);
  x1 << 0.0;
  Eigen::VectorXd x2(1);
  x2 << 1.5;

  Eigen::VectorXd res = AnyConstCast(graphPiece->Evaluate(x1,x2).at(0));

  Eigen::VectorXd v1(1);
  v1(0) = gaussDens->LogDensity(x1);

  Eigen::VectorXd v2(1);
  v2(0) = uniformDens->LogDensity(x2);

  EXPECT_DOUBLE_EQ(prodPiece->LogDensity(v1,v2), res(0));
}
