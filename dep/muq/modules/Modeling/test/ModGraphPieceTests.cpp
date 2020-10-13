#include "gtest/gtest.h"

#include "MUQ/Modeling/WorkGraph.h"
#include "MUQ/Modeling/ConstantVector.h"
#include "MUQ/Modeling/IdentityPiece.h"

#include "MUQ/Modeling/ModPiece.h"
#include "MUQ/Modeling/ModGraphPiece.h"

#include "WorkPieceTestClasses.h"

using namespace muq::Modeling;
using namespace std;

// define a simple single input forward model
class SquareMod : public ModPiece {
public:

  /** Constructor taking vector dimension and resizing the State.*/
  SquareMod(int dim) : ModPiece(dim*Eigen::VectorXi::Ones(1), dim*Eigen::VectorXi::Ones(1)) {}

  virtual ~SquareMod() = default;

private:

  virtual void EvaluateImpl(ref_vector<Eigen::VectorXd> const& input) override
  {
    outputs.resize(1);
    outputs.at(0) = input.at(0).get().array().square();
  }

  virtual void GradientImpl(unsigned int                const  outputDimWrt,
                            unsigned int                const  inputDimWrt,
                            ref_vector<Eigen::VectorXd> const& input,
                            Eigen::VectorXd             const& sensitivity) override
  {
    gradient = 2.0 * sensitivity.array() * input.at(0).get().array();
  }

  virtual void JacobianImpl(unsigned int                const  outputDimWrt,
                            unsigned int                const  inputDimWrt,
                            ref_vector<Eigen::VectorXd> const& input) override
  {
    jacobian = 2.0*input.at(0).get().asDiagonal();
  }

  virtual void ApplyJacobianImpl(unsigned int                const  outputDimWrt,
                                 unsigned int                const  inputDimWrt,
                                 ref_vector<Eigen::VectorXd> const& input,
                                 Eigen::VectorXd             const& vec) override
  {
    jacobianAction = 2.0*input.at(0).get().array()*vec.array();
  };

};


// define a simple single input forward model
class SinSumMod : public ModPiece {
public:

  /** Constructor taking vector dimension and resizing the State.*/
  SinSumMod(int dim) : ModPiece(dim*Eigen::VectorXi::Ones(2), dim*Eigen::VectorXi::Ones(1)) {}

  virtual ~SinSumMod() = default;

private:

  virtual void EvaluateImpl(ref_vector<Eigen::VectorXd> const& input) override
  {
    Eigen::VectorXd out = Eigen::VectorXd::Zero(outputSizes(0));

    for (unsigned int j = 0; j < input.size(); ++j) {
      for (int i = 0; i < outputSizes(0); ++i) {
        out[i] += sin(input.at(j)(i));
      }
    }

    outputs.resize(1);
    outputs.at(0) = out;
  }

  virtual void GradientImpl(unsigned int                const  outputDimWrt,
                            unsigned int                const  inputDimWrt,
                            ref_vector<Eigen::VectorXd> const& input,
                            Eigen::VectorXd             const& sensitivity) override
  {
    gradient = Eigen::VectorXd::Zero(outputSizes(0));

    for (int i = 0; i < outputSizes(0); ++i)
      gradient(i) = sensitivity(i) * cos(input.at(inputDimWrt)(i));
  }

  virtual void JacobianImpl(unsigned int                const  outputDimWrt,
                            unsigned int                const  inputDimWrt,
                            ref_vector<Eigen::VectorXd> const& input) override
  {
    jacobian = Eigen::MatrixXd::Zero(outputSizes(0), inputSizes(inputDimWrt));
    for (int i = 0; i < outputSizes(0); ++i)
      jacobian(i, i) = cos(input.at(inputDimWrt)(i));
  }

  virtual void ApplyJacobianImpl(unsigned int                const  outputDimWrt,
                                 unsigned int                const  inputDimWrt,
                                 ref_vector<Eigen::VectorXd> const& input,
                                 Eigen::VectorXd             const& vec) override
  {
    jacobianAction = Eigen::VectorXd::Zero(outputSizes(0));
    for(int i=0; i<outputSizes(0); ++i)
      jacobianAction(i) = cos(input.at(inputDimWrt)(i))*vec(i);
  };

};


TEST(Modelling_ModGraphPiece, MatchInputs)
{
  auto myGraph = make_shared<WorkGraph>();

  // add nodes
  myGraph->AddNode(make_shared<SinSumMod>(2), "f2");
  myGraph->AddNode(make_shared<SinSumMod>(2), "f1");
  myGraph->AddNode(make_shared<SquareMod>(2), "x");
  myGraph->AddNode(make_shared<SquareMod>(2), "y1");
  myGraph->AddNode(make_shared<SquareMod>(2), "y2");

  // add connectivity
  myGraph->AddEdge("x", 0, "f2", 0);
  myGraph->AddEdge("y1", 0, "f1", 0);
  myGraph->AddEdge("y2", 0, "f1", 1);
  myGraph->AddEdge("f1", 0, "f2", 1);

  auto piece1 = myGraph->CreateModPiece("f1");
  auto piece2 = myGraph->CreateModPiece("f2");

  EXPECT_EQ(2, piece1->inputSizes.size());
  EXPECT_EQ(3, piece2->inputSizes.size());

  std::vector<int> sharedIns = piece1->MatchInputs(piece2);
  EXPECT_EQ(3, sharedIns.size());
  EXPECT_EQ(0, sharedIns.at(0));
  EXPECT_EQ(1, sharedIns.at(1));
  EXPECT_EQ(-1, sharedIns.at(2));

  sharedIns = piece2->MatchInputs(piece1);
  EXPECT_EQ(2, sharedIns.size());
  EXPECT_EQ(0, sharedIns.at(0));
  EXPECT_EQ(1, sharedIns.at(1));
}


TEST(Modelling_ModGraphPiece, BasicTest)
{
  auto myGraph = make_shared<WorkGraph>();

  // add nodes
  myGraph->AddNode(make_shared<SinSumMod>(2), "f2");
  myGraph->AddNode(make_shared<SinSumMod>(2), "f1");
  myGraph->AddNode(make_shared<SquareMod>(2), "x");
  myGraph->AddNode(make_shared<SquareMod>(2), "y1");
  myGraph->AddNode(make_shared<SquareMod>(2), "y2");

  // add connectivity
  myGraph->AddEdge("x", 0, "f2", 0);
  myGraph->AddEdge("y1", 0, "f1", 0);
  myGraph->AddEdge("y2", 0, "f1", 1);
  myGraph->AddEdge("f1", 0, "f2", 1);

  boost::any anyVal1 = Eigen::VectorXd::Constant(2,0.1).eval();
  boost::any anyVal2 = Eigen::VectorXd::Constant(2,0.2).eval();

  myGraph->BindNode("y1", {anyVal1});
  myGraph->BindNode("y2", {anyVal2});

  auto graphMod = myGraph->CreateModPiece("f2");

  // make sure this modpiece is the size we expect
  EXPECT_EQ(1, graphMod->inputSizes.size());
  EXPECT_EQ(2, graphMod->inputSizes(0));
  EXPECT_EQ(2, graphMod->outputSizes(0));

  // evaluation testing
  Eigen::VectorXd input  = Eigen::VectorXd::Ones(2);
  Eigen::VectorXd output = graphMod->Evaluate(input).at(0);

  EXPECT_DOUBLE_EQ(sin(input[0] * input[0]) + sin(sin(0.1) + sin(0.2)), output[0]);
  EXPECT_DOUBLE_EQ(sin(input[1] * input[1]) + sin(sin(0.1) + sin(0.2)), output[1]);

  // // gradient testing (same as J^T*x)
  // Eigen::VectorXd grad = GraphMod->Gradient(input, Eigen::VectorXd::Ones(2), 0);
  //
  // EXPECT_DOUBLE_EQ(2.0 * input[0] * cos(input[0] * input[0]), grad[0]);
  // EXPECT_DOUBLE_EQ(2.0 * input[1] * cos(input[1] * input[1]), grad[1]);
  //
  // // Jacobian action testing (same as J*x)
  // Eigen::VectorXd input2 = input + 1e-2 * Eigen::VectorXd::Random(2);
  // Eigen::VectorXd linOut = GraphMod->JacobianAction(input, input2, 0);
  // EXPECT_DOUBLE_EQ(2.0 * input[0] * cos(input[0] * input[0]) * input2[0], linOut[0]);
  // EXPECT_DOUBLE_EQ(2.0 * input[1] * cos(input[1] * input[1]) * input2[1], linOut[1]);
  //
  // // Jacobian testing
  // Eigen::MatrixXd jacOut = GraphMod->Jacobian(input, 0);
  // EXPECT_DOUBLE_EQ(2.0 * input[0] * cos(input[0] * input[0]), jacOut(0, 0));
  // EXPECT_DOUBLE_EQ(0,                                         jacOut(1, 0));
  // EXPECT_DOUBLE_EQ(0,                                         jacOut(0, 1));
  // EXPECT_DOUBLE_EQ(2.0 * input[1] * cos(input[1] * input[1]), jacOut(1, 1));
}


// TEST(Modelling_ModGraphPiece, NodeOrdering)
// {
//   auto myGraph = make_shared<ModGraph>();
//
//   // add nodes
//   myGraph->AddNode(make_shared<sinSumMod>(2), "f1");
//   myGraph->AddNode(make_shared<squareMod>(2), "z");
//   myGraph->AddNode(make_shared<squareMod>(2), "x1");
//
//   // add connectivity
//   myGraph->AddEdge("z", "f1", 0);
//   myGraph->AddEdge("x1", "f1", 1);
//
//   myGraph->writeGraphViz("results/tests/GraphViz/NodeOrderTest.pdf");
//
//   auto GraphMod = ModGraphPiece::Create(myGraph, "f1");
//   GraphMod->writeGraphViz("results/tests/GraphViz/NodeOrderPieceTest.pdf");
//
//   std::vector<std::string> newOrder(2);
//   newOrder[0] = "z";
//   newOrder[1] = "x1";
//   auto GraphMod2 = ModGraphPiece::Create(myGraph, "f1", newOrder);
//   GraphMod->writeGraphViz("results/tests/GraphViz/NodeOrderPieceTest2.pdf");
// }
//
//
// TEST(Modelling_ModGraphPiece, DiamondTest)
// {
//   auto myGraph = make_shared<ModGraph>();
//
//   // add nodes
//   myGraph->AddNode(make_shared<sinSumMod>(2), "f1");
//   myGraph->AddNode(make_shared<squareMod>(2), "x");
//   myGraph->AddNode(make_shared<squareMod>(2), "y1");
//   myGraph->AddNode(make_shared<squareMod>(2), "y2");
//
//   // add connectivity
//   myGraph->AddEdge("x", "y1", 0);
//   myGraph->AddEdge("x", "y2", 0);
//   myGraph->AddEdge("y1", "f1", 0);
//   myGraph->AddEdge("y2", "f1", 1);
//   myGraph->writeGraphViz("results/tests/GraphViz/DiamondTest.pdf");
//   auto GraphMod = ModGraphPiece::Create(myGraph, "f1");
//   GraphMod->writeGraphViz("results/tests/GraphViz/DiamondPieceTest.pdf");
//
//   // make sure this modpiece is the size we expect
//   EXPECT_EQ(1, GraphMod->inputSizes.size());
//   EXPECT_EQ(2, GraphMod->inputSizes[0]);
//   EXPECT_EQ(2, GraphMod->outputSize);
//
//   // evaluation testing
//   Eigen::VectorXd input  = 0.5 * Eigen::VectorXd::Ones(2);
//   Eigen::VectorXd output = GraphMod->Evaluate(input);
//
//   EXPECT_DOUBLE_EQ(2.0 * sin(pow(input[0], 4.0)), output[0]);
//   EXPECT_DOUBLE_EQ(2.0 * sin(pow(input[1], 4.0)), output[1]);
//
//
//   // gradient testing
//   Eigen::VectorXd grad = GraphMod->Gradient(input, Eigen::VectorXd::Ones(2), 0);
//
//   EXPECT_DOUBLE_EQ(8.0 * pow(input[0], 3) * cos(pow(input[0], 4.0)), grad[0]);
//   EXPECT_DOUBLE_EQ(8.0 * pow(input[1], 3) * cos(pow(input[1], 4.0)), grad[1]);
//
//   // jacobian action testing
//   Eigen::VectorXd input2 = input + 1e-2 * Eigen::VectorXd::Random(2);
//   Eigen::VectorXd linOut = GraphMod->JacobianAction(input, input2, 0);
//   EXPECT_DOUBLE_EQ(8.0 * pow(input[0], 3) * cos(pow(input[0], 4.0)) * input2[0], linOut[0]);
//   EXPECT_DOUBLE_EQ(8.0 * pow(input[1], 3) * cos(pow(input[1], 4.0)) * input2[1], linOut[1]);
//
//   // jacobian testing
//   Eigen::MatrixXd jacOut = GraphMod->Jacobian(input, 0);
//   EXPECT_DOUBLE_EQ(8.0 * pow(input[0], 3) * cos(pow(input[0], 4.0)), jacOut(0, 0));
//   EXPECT_DOUBLE_EQ(0,                                                jacOut(1, 0));
//   EXPECT_DOUBLE_EQ(0,                                                jacOut(0, 1));
//   EXPECT_DOUBLE_EQ(8.0 * pow(input[1], 3) * cos(pow(input[1], 4.0)), jacOut(1, 1));
//
//   // Hessian testing
//   Eigen::MatrixXd hessOut = GraphMod->Hessian(input, Eigen::VectorXd::Ones(2), 0);
//   const double    diag    = -32.0 *
//                             pow(input[0],
//                                 6.0) * sin(pow(input[0], 4.0)) + 24 * pow(input[0], 2.0) * cos(pow(input[0], 4.0));
//
//   EXPECT_DOUBLE_EQ(hessOut(0, 0), diag);
//   EXPECT_DOUBLE_EQ(hessOut(1, 1), diag);
//   EXPECT_DOUBLE_EQ(hessOut(1, 0), 0);
//   EXPECT_DOUBLE_EQ(hessOut(0, 1), 0);
// }
//
// TEST(Modelling_ModGraphPiece, UnionTest)
// {
//   auto a = make_shared<ModGraph>();
//   auto b = make_shared<ModGraph>();
//   auto sourceNode = make_shared<squareMod>(2);
//
//   // add nodes
//   a->AddNode(make_shared<sinSumMod>(2), "af2");
//   a->AddNode(make_shared<sinSumMod>(2), "af1");
//   a->AddNode(sourceNode, "x");
//   a->AddNode(make_shared<squareMod>(2), "ay1");
//   a->AddNode(make_shared<squareMod>(2), "ay2");
//
//   // add connectivity
//   a->AddEdge("x", "af2", 0);
//   a->AddEdge("ay1", "af1", 0);
//   a->AddEdge("ay2", "af1", 1);
//   a->AddEdge("af1", "af2", 1);
//
//     // add nodes
//   b->AddNode(make_shared<sinSumMod>(2), "bf2");
//   b->AddNode(make_shared<sinSumMod>(2), "bf1");
//   b->AddNode(sourceNode, "x");
//   b->AddNode(make_shared<squareMod>(2), "by1");
//   b->AddNode(make_shared<squareMod>(2), "by2");
//
//   // add connectivity
//   b->AddEdge("x", "bf2", 0);
//   b->AddEdge("by1", "bf1", 0);
//   b->AddEdge("by2", "bf1", 1);
//   b->AddEdge("bf1", "bf2", 1);
//
//   auto unionGraph = ModGraph::FormUnion(a,b);
//   unionGraph->writeGraphViz("results/tests/UnionTest.pdf");
//
//   EXPECT_EQ(5, unionGraph->NumInputs()); //this tests that the "x" node is not duplicated
//   EXPECT_EQ(2, unionGraph->NumOutputs()); //and we have both model outputs
// }
//
// TEST(Modelling_ModGraphPiece, UnionNameClashDeath)
// {
//   auto a = make_shared<ModGraph>();
//   auto b = make_shared<ModGraph>();
//
//   // add nodes - x is distinct but has the same name
//   a->AddNode(make_shared<sinSumMod>(2), "af2");
//   a->AddNode(make_shared<sinSumMod>(2), "af1");
//   a->AddNode(make_shared<squareMod>(2), "x");
//
//   b->AddNode(make_shared<squareMod>(2), "x");
//   b->AddNode(make_shared<squareMod>(2), "by1");
//   b->AddNode(make_shared<squareMod>(2), "by2");
//
//   ASSERT_DEATH(ModGraph::FormUnion(a,b), "(result->GetNodeModel(currentName) == b->ModelGraph[v]->piece)*");
//
// }
//
// TEST(Modelling_ModGraphPiece, UnionEdgeClashDeath)
// {
// 	//Both graphs share a node correctly, but both try to provide an input, so we don't know how to
// 	//uniquely resolve it and hence assert out
//   auto a = make_shared<ModGraph>();
//   auto b = make_shared<ModGraph>();
//   auto sourceNode = make_shared<squareMod>(2);
//
//   // add nodes
//   a->AddNode(make_shared<sinSumMod>(2), "af2");
//   a->AddNode(make_shared<sinSumMod>(2), "af1");
//   a->AddNode(sourceNode, "x");
//   a->AddNode(make_shared<squareMod>(2), "ay1");
//   a->AddNode(make_shared<squareMod>(2), "ay2");
//
//   // add connectivity
//   a->AddEdge("x", "af2", 0);
//   a->AddEdge("ay1", "af1", 0);
//   a->AddEdge("ay1", "x", 0); //clashes for the first input of x
//   a->AddEdge("ay2", "af1", 1);
//   a->AddEdge("af1", "af2", 1);
//
//     // add nodes
//   b->AddNode(make_shared<sinSumMod>(2), "bf2");
//   b->AddNode(make_shared<sinSumMod>(2), "bf1");
//   b->AddNode(sourceNode, "x");
//   b->AddNode(make_shared<squareMod>(2), "by1");
//   b->AddNode(make_shared<squareMod>(2), "by2");
//
//   // add connectivity
//   b->AddEdge("x", "bf2", 0);
//   b->AddEdge("by1", "bf1", 0);
//   b->AddEdge("by2", "x", 0); //clashes for the first input of x
//   b->AddEdge("by2", "bf1", 1);
//   b->AddEdge("bf1", "bf2", 1);
//
//   ASSERT_DEATH(ModGraph::FormUnion(a,b), "(result->ModelGraph[*e_result]->GetDim() != b->ModelGraph[e]->GetDim())*");
// }
