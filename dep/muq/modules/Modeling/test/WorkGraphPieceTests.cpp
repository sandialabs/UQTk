#include "gtest/gtest.h"

#include "MUQ/Modeling/WorkGraph.h"
#include "MUQ/Modeling/ConstantPiece.h"
#include "MUQ/Modeling/IdentityPiece.h"

#include "MUQ/Modeling/WorkPiece.h"
#include "MUQ/Modeling/WorkGraphPiece.h"
#include "MUQ/Modeling/LinearAlgebra/AnyAlgebra.h"

#include "WorkPieceTestClasses.h"

using namespace muq::Modeling;

TEST(WorkGraphPiece, FixedInOutNum) {
  // the input types
  std::vector<std::string> inTypes({typeid(std::string).name(), typeid(double).name(), typeid(std::shared_ptr<AnObject>).name()});
  // the output types
  std::vector<std::string> outTypes({typeid(std::string).name(), typeid(double).name()});

  // create the test WorkPiece
  auto test0 = std::make_shared<FixedInOutMod>(inTypes, outTypes);
  auto test1 = std::make_shared<FixedInOutMod>(inTypes, outTypes);

  // create a version of "AnObject"
  auto obj = std::make_shared<AnObject>(2.0);
  obj->flag = true;

  // make a constant parameter
  auto test2 = std::make_shared<ConstantPiece>(obj, 1);
  auto test3 = std::make_shared<IdentityPiece>(outTypes);

  // create and empty graph
  auto graph = std::make_shared<WorkGraph>();

  // add WorkPieces to the graph
  graph->AddNode(test0, "test 0");
  graph->AddNode(test1, "test 1");
  graph->AddNode(test2, "test 2");
  graph->AddNode(test3, "test 3");
  EXPECT_EQ(graph->NumNodes(), 4);

  // connect test0 to test1
  graph->AddEdge("test 1", 0, "test 0", 0);
  graph->AddEdge("test 1", 1, "test 0", 1);
  graph->AddEdge("test 2", 0, "test 1", 2);
  graph->AddEdge("test 3", 0, "test 1", 0);
  graph->AddEdge("test 3", 1, "test 1", 1);
  EXPECT_EQ(graph->NumEdges(), 5);

  graph->Visualize("modules/Modeling/test/WorkGraphVisualizations/FixedInOutNum_WorkPiece.pdf");

  // create a WorkPiece whose outs are from "test 0"
  auto work = graph->CreateWorkPiece("test 0");

  // output of the WorkGraphPiece
  auto outs = work->Evaluate(obj, (std::string)"string", 2.0);

  // check the outputs
  EXPECT_EQ(outs.size(), 2);
  EXPECT_TRUE(boost::any_cast<std::string>(outs[0]).compare("string")==0);
  EXPECT_DOUBLE_EQ(boost::any_cast<double>(outs[1]), 2.0);
}

/*TEST(WorkGraphPiece, Constant) {
  auto graph = std::make_shared<WorkGraph>();

    auto ip = std::make_shared<IdentityPiece>(1);

    double a = 1.0;
    double b = 2.0;
    auto input = std::vector<boost::any>({boost::any(a), boost::any(b)});

    auto cp = std::make_shared<ConstantPiece>(input);

    graph->AddNode(cp,"x");
    graph->AddNode(ip,"y");
    graph->AddEdge("x",0,"y",0);
    graph->CreateWorkPiece("y");
    }*/
