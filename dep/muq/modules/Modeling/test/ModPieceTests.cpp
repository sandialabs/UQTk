// google testing library
#include "gtest/gtest.h"

// #include "MUQ/Modelling/Model.h"
#include "MUQ/Modeling/ModPiece.h"

// #include "MUQ/Modelling/ModParameter.h"
// #include "MUQ/Modelling/AnalyticFunctions/SymbolicHeaders.h"

using namespace muq::Modeling;

// define a simple single input forward model
class SinMod : public ModPiece {
public:

  /** Constructor taking vector dimension and resizing the State.*/
  SinMod(int dim) : ModPiece(dim*Eigen::VectorXi::Ones(1), dim*Eigen::VectorXi::Ones(1)) {}

  virtual ~SinMod() = default;

private:

  virtual void EvaluateImpl(ref_vector<Eigen::VectorXd> const& input) override
  {
    outputs.resize(1);
    outputs.at(0) = sin(input.at(0).get().array());
  }

  virtual void GradientImpl(unsigned int                const  outputDimWrt,
                            unsigned int                const  inputDimWrt,
                            ref_vector<Eigen::VectorXd> const& input,
                            Eigen::VectorXd             const& sensitivity) override
  {
    gradient = sensitivity.array() * cos(input.at(0).get().array());
  }

  virtual void JacobianImpl(unsigned int                const  outputDimWrt,
                            unsigned int                const  inputDimWrt,
                            ref_vector<Eigen::VectorXd> const& input) override
  {
    jacobian = Eigen::MatrixXd(input.at(0).get().array().cos().matrix().asDiagonal());
  }
};

// define a simple two input forward model
class SumMod : public ModPiece {
public:

  /** Constructor taking vector dimension and resizing the State.*/
  SumMod(int dim) : ModPiece(dim*Eigen::VectorXi::Ones(2), dim*Eigen::VectorXi::Ones(1)){};

  virtual ~SumMod() = default;

private:

  virtual void EvaluateImpl(ref_vector<Eigen::VectorXd> const& input) override
  {
    outputs.resize(1);
    outputs.at(0) = input.at(0).get() + input.at(1).get();
  }

  virtual void GradientImpl(unsigned                    const  outputDimWrt,
                            unsigned                    const  inputDimWrt,
                            ref_vector<Eigen::VectorXd> const& input,
                            Eigen::VectorXd             const& sensitivity) override
  {
    gradient = sensitivity;
  }
};

// define a simple single input model with different input/output dimensions
class FooMod : public ModPiece {
public:

  /** Constructor taking vector dimensions */
  FooMod(int dimIn, int dimOut) : ModPiece(dimIn * Eigen::VectorXi::Ones(1), dimOut* Eigen::VectorXi::Ones(1)){};

  virtual ~FooMod() = default;

private:

  virtual void EvaluateImpl(ref_vector<Eigen::VectorXd> const& input) override
  {
    Eigen::ArrayXd shifted(inputSizes[0]);

    for (int i = 0; i < inputSizes[0]; ++i) {
      shifted(i) = input[0]((i + 1) % inputSizes[0]);
    }

    Eigen::VectorXd output(outputSizes(0));
    for (int i = 0; i < outputSizes(0); ++i) {
      output[i] = (input[0].get().array() * shifted).sum() + (input[0].get().array() * input[0].get().array() * input[0].get().array()).sum();
    }

    outputs.resize(1);
    outputs.at(0) = output;
  }

  virtual void GradientImpl(unsigned int                const  outputDimWrt,
                            unsigned int                const  inputDimWrt,
                            ref_vector<Eigen::VectorXd> const& input,
                            Eigen::VectorXd             const& sensitivity) override
  {
    Eigen::ArrayXd shifted(inputSizes[0]);

    for (int i = 0; i < inputSizes[0]; ++i) {
      shifted(i) = input[0]((i + 1) % inputSizes[0]);
    }

    gradient = Eigen::VectorXd::Zero(inputSizes[0]);
    for (int i = 0; i < inputSizes[0]; ++i) {
      int j = ((i - 1) < 0) ? inputSizes[0] - 1 : i - 1;
      gradient[i] = sensitivity.dot((input[0](j) + shifted(i) + 3.0 * input[0](i) * input[0](
                                       i)) * Eigen::VectorXd::Ones(outputSizes(0)));
    }

  }

  // virtual Eigen::MatrixXd HessianImpl(std::vector<Eigen::VectorXd> const& input,
  //                                     Eigen::VectorXd const             & sensitivity,
  //                                     int const                           inputDimWrt) override
  // {
  //   Eigen::MatrixXd Hessian = Eigen::MatrixXd::Zero(inputSizes[0], inputSizes[0]);
  //
  //   for (int out = 0; out < outputSize; ++out) {
  //     for (int i = 0; i < inputSizes[0]; ++i) {
  //       Hessian(i, i) += 6.0 * input[0](i) * sensitivity(out);
  //
  //       int j1 = ((i - 1) < 0) ? inputSizes[0] - 1 : i - 1;
  //       Hessian(i, j1) += 1.0 * sensitivity(out);
  //
  //       int j2 = ((i + 1) == inputSizes[0]) ? 0 : i + 1;
  //       Hessian(i, j2) += 1.0 * sensitivity(out);
  //     }
  //   }
  //
  //   return Hessian;
  // }
};

// // sum up the input vector
// class ReduceMod : public ModPiece<>{
// public:
//     ReduceMod(int dim){
//         DimIn[0] = dim;
//         DimOut = 1;
//         State.resize(DimOut);
//     };
//
//
//     virtual void UpdateBase(const Eigen::VectorXd &VectorIn0){
//         State[0] = VectorIn0.sum();
//     };
//
//     virtual void GradBase(const Eigen::VectorXd &VectorIn, const Eigen::VectorXd &SensIn, Eigen::VectorXd &GradVec,
// int dim){
//         GradVec = SensIn[0]*Eigen::VectorXd::Ones(DimIn[dim]);
//     };
//
//     // this is just a random function that we may want to call -- it will just add 1 to the state
//     virtual void AddOne(){
//         State += Eigen::VectorXd::Ones(DimOut);
//     };
//
// };

TEST(ModPieceTest, NamingTest)
{
  // set the vector size
  int dim = 10;

  // create the models
  auto testSinMod = std::make_shared<SinMod>(dim);
  auto testFooMod = std::make_shared<FooMod>(dim, dim);

  EXPECT_STREQ(testSinMod->Name().substr(0,6).c_str(), "SinMod");
  EXPECT_STREQ(testFooMod->Name().substr(0,6).c_str(), "FooMod");

  testSinMod->SetName("RenamedSinMod");
  EXPECT_STREQ(testSinMod->Name().c_str(), "RenamedSinMod");
}

TEST(ModPieceTest, UniqueNames)
{
  int dim = 10;

  // create the forward model object
  auto mod1 = std::make_shared<SinMod>(dim);
  auto mod2 = std::make_shared<SinMod>(dim);

  //We just want to make sure that even though we didn't provide names, and their types are identical,
  //the default names they are given are different.
  EXPECT_NE(0, mod1->Name().compare(mod2->Name()));
}

TEST(ModPieceTest, FiniteDifference)
{
  // set the vector size
  int dim = 10;

  // test the model evaluation
  Eigen::VectorXd vecIn  = Eigen::VectorXd::Ones(dim);
  Eigen::VectorXd sensIn = Eigen::VectorXd::Ones(dim);

  // create the model
  auto testSinMod = std::make_shared<SinMod>(dim);

  std::vector<Eigen::VectorXd> vecVec(1,vecIn);

  Eigen::VectorXd gradVec  = testSinMod->GradientByFD(0,0, vecVec, sensIn);
  Eigen::VectorXd trueGrad = vecIn.array().cos();

  EXPECT_EQ(dim+1, testSinMod->GetNumCalls("Evaluate"));

  for(int i=0; i<dim; ++i)
    EXPECT_NEAR(trueGrad(i), gradVec(i),1e-4);

  int numPrevEvals = testSinMod->GetNumCalls("Evaluate");

  Eigen::MatrixXd fdJac  = testSinMod->JacobianByFD(0,0,vecVec);

  Eigen::MatrixXd trueJac = testSinMod->Jacobian(0,0,vecVec);
  int numJacEvals = testSinMod->GetNumCalls("Evaluate");
  EXPECT_EQ(dim+1,numJacEvals - numPrevEvals);

  for(int j=0; j<dim; ++j){
    for(int i=0; i<dim; ++i){
      EXPECT_NEAR(fdJac(i,j), trueJac(i,j), 1e-4);
    }
  }

  {
    Eigen::VectorXd targetVec = 0.8*Eigen::VectorXd::Ones(dim);
    Eigen::VectorXd jacX = testSinMod->ApplyJacobianByFD(0,0, vecVec, targetVec);
    EXPECT_EQ(2, testSinMod->GetNumCalls("Evaluate")-numJacEvals);

    Eigen::VectorXd trueJacX = testSinMod->ApplyJacobian(0,0, vecVec,targetVec);

    for(int i=0; i<jacX.size(); ++i)
      EXPECT_NEAR(jacX(i), trueJacX(i),1e-4);
  }

}

TEST(ModellingModelTest, SingleEval)
{
  // set the vector size
  int dim = 10;

  // create the model
  auto testMod = std::make_shared<SinMod>(dim);

  // test the model evaluation
  Eigen::VectorXd vecIn = Eigen::VectorXd::Ones(dim);


  Eigen::VectorXd vecOut     = testMod->Evaluate(vecIn).at(0);
  Eigen::VectorXd trueOutput = vecIn.array().sin();
  for(int i=0; i<vecOut.size(); ++i)
    EXPECT_NEAR(trueOutput(i), vecOut(i), 1e-14);

  //try it again
  vecOut = testMod->Evaluate(vecIn).at(0);
  for(int i=0; i<vecOut.size(); ++i)
    EXPECT_NEAR(trueOutput(i), vecOut(i), 1e-14);

  //and now at a different point
  vecIn      = Eigen::VectorXd::Ones(dim) * 2;
  vecOut     = testMod->Evaluate(vecIn).at(0);
  trueOutput = vecIn.array().sin();
  for(int i=0; i<vecOut.size(); ++i)
    EXPECT_NEAR(trueOutput(i), vecOut(i), 1e-14);

}

TEST(ModPieceTest, SingleGrad)
{
  // set the vector size
  int dim = 10;

  // create the model
  auto testMod = std::make_shared<SinMod>(dim);

  // test the model evaluation
  Eigen::VectorXd vecIn = Eigen::VectorXd::Ones(dim);

  Eigen::VectorXd sensIn   = Eigen::VectorXd::Ones(dim);
  Eigen::VectorXd gradVec  = testMod->Gradient(0, 0, vecIn, sensIn);
  Eigen::VectorXd trueGrad = vecIn.array().cos();

  for(int i=0; i<gradVec.size(); ++i)
    EXPECT_NEAR(trueGrad(i), gradVec(i), 1e-14);

  //try it again
  gradVec = testMod->Gradient(0,0,vecIn, sensIn);
  for(int i=0; i<gradVec.size(); ++i)
    EXPECT_NEAR(trueGrad(i), gradVec(i), 1e-14);

  //and now at a different point
  vecIn    = Eigen::VectorXd::Ones(dim) * 2;
  gradVec  = testMod->Gradient(0,0,vecIn, sensIn);
  trueGrad = vecIn.array().cos();
  for(int i=0; i<gradVec.size(); ++i)
    EXPECT_NEAR(trueGrad(i), gradVec(i), 1e-14);
}

TEST(ModPieceTest, SingleJacobian)
{
  // set the vector size
  int dim = 3;

  // create the model
  auto testMod = std::make_shared<SinMod>(dim);

  // test the model evaluation
  Eigen::VectorXd vecIn(3);

  vecIn << 1, 2, .3;

  Eigen::MatrixXd computedJacobian = testMod->Jacobian(0,0,vecIn);
  Eigen::MatrixXd correctJacobian(3, 3);
  correctJacobian << 0.540302305868140, 0, 0,
                     0,  -0.416146836547142, 0,
                     0, 0, 0.955336489125606;

  for(int j=0; j<dim; ++j){
    for(int i=0; i<dim; ++i)
      EXPECT_NEAR(correctJacobian(i,j), computedJacobian(i,j), 1e-14);
  }

  //do it again
  computedJacobian = testMod->Jacobian(0,0,vecIn);
  for(int j=0; j<dim; ++j){
    for(int i=0; i<dim; ++i)
      EXPECT_NEAR(correctJacobian(i,j), computedJacobian(i,j), 1e-14);
  }
}
