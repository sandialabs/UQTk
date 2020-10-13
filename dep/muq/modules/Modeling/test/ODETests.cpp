#include "gtest/gtest.h"

#include <Eigen/Core>

#include "MUQ/Modeling/ODE.h"

namespace pt = boost::property_tree;
using namespace muq::Modeling;

class RHS : public ModPiece {
public:

  /// Constructor
  inline RHS() : ModPiece(Eigen::Vector2i(2, 1), Eigen::VectorXi::Constant(1, 2)) {}

  inline virtual ~RHS() = default;

private:

  inline virtual void EvaluateImpl(ref_vector<Eigen::VectorXd> const& inputs) override {
    // get the state vector
    const Eigen::VectorXd& state = inputs[0];

    // get the parameter (spring constant)
    const double k = inputs[1].get() [0];

    // set the output
    outputs.resize(1);
    outputs[0] = Eigen::VectorXd(2);

    outputs[0](0) = state(1);
    outputs[0](1) = -k*state(0);
  }

  inline virtual void ApplyJacobianImpl(unsigned int const wrtOut, unsigned int const wrtIn, ref_vector<Eigen::VectorXd> const& inputs, Eigen::VectorXd const& vec) override {
    // there is only one output
    assert(wrtOut==0);

    if( wrtIn==0 ) { // wrt the state
      assert(vec.size()==2);

      // get the parameter (spring constant)
      const double k = inputs[1].get() [0];

      jacobianAction = Eigen::VectorXd(2);

      jacobianAction(0) = vec(1);
      jacobianAction(1) = -k*vec(0);
    } else if( wrtIn==1 ) {
      assert(vec.size()==1);

      // get the state vector
      const Eigen::VectorXd& state = inputs[0];

      jacobianAction = Eigen::VectorXd(2);

      jacobianAction(0) = 0.0;
      jacobianAction(1) = -state(0)*vec(0);
    }
  }

    inline virtual void JacobianImpl(unsigned int const wrtOut, unsigned int const wrtIn, ref_vector<Eigen::VectorXd> const& inputs) override {
    // there is only one output
    assert(wrtOut==0);

    if( wrtIn==0 ) { // wrt the state
      // get the parameter (spring constant)
      const double k = inputs[1].get() [0];

      jacobian = Eigen::MatrixXd::Zero(2, 2);

      jacobian(0, 1) = 1.0;
      jacobian(1, 0) = -k;
    } else if( wrtIn==1 ) {
      // get the state vector
      const Eigen::VectorXd& state = inputs[0];

      jacobian = Eigen::MatrixXd::Zero(2, 1);
      jacobian(1, 0) = -state(0);
    }
  }
};

// test the right hand side WorkPiece
TEST(ODEExample, RHS) {
  // create the right hand side
  auto rhs = std::make_shared<RHS>();

  // inputs
  const Eigen::VectorXd state = Eigen::Vector2d(1.2, 0.8);
  const Eigen::VectorXd k = Eigen::VectorXd::Constant(1, 2.25);

  { // test evaluate
    const std::vector<Eigen::VectorXd>& result = rhs->Evaluate(state, k);
    EXPECT_EQ(result.size(), 1);

    EXPECT_EQ(result[0].size(), 2);
    EXPECT_EQ(result[0](0), state(1));
    EXPECT_EQ(result[0](1), -k(0)*state(0));
  }

  // expected jacobians
  Eigen::MatrixXd expectedJac0 = Eigen::MatrixXd::Zero(2, 2);
  expectedJac0(0, 1) = 1.0;
  expectedJac0(1, 0) = -k(0);

  Eigen::MatrixXd expectedJac1 = Eigen::MatrixXd::Zero(2, 1);
  expectedJac1(1, 0) = -state(0);

  { // test jacobian
    const Eigen::MatrixXd& jac0 = rhs->Jacobian(0, 0, state, k);

    EXPECT_EQ(jac0.rows(), 2); // check the number of rows
    EXPECT_EQ(jac0.cols(), 2); // check the number of cols
    for( unsigned int i=0; i<2; ++i ) {
      for( unsigned int j=0; j<2; ++j ) {
	       EXPECT_DOUBLE_EQ(jac0(i, j), expectedJac0(i, j));
      }
    }

    const Eigen::MatrixXd& jac1 = rhs->Jacobian(0, 1, state, k);

    EXPECT_EQ(jac1.rows(), 2); // check the number of rows
    EXPECT_EQ(jac1.cols(), 1); // check the number of cols
    for( unsigned int i=0; i<2; ++i ) {
      EXPECT_DOUBLE_EQ(jac1(i, 0), expectedJac1(i,0));
    }
  }

  { // test jacobian action
    const Eigen::VectorXd vec0 = Eigen::VectorXd::Ones(2);

    const Eigen::VectorXd& jacact0 = rhs->ApplyJacobian(0, 0, state, k, vec0);

    // expected value
    const Eigen::Vector2d expectedJacact0 = expectedJac0*vec0;

    EXPECT_EQ(jacact0.size(), 2);
    for( unsigned int i=0; i<2; ++i ) {
      EXPECT_DOUBLE_EQ(jacact0(i), expectedJacact0(i));
    }

    const Eigen::VectorXd vec1 = Eigen::VectorXd::Ones(1);

    const Eigen::VectorXd& jacact1 = rhs->ApplyJacobian(0, 1, state, k, vec1);

    // expected value
    const Eigen::Vector2d expectedJacact1 = expectedJac1*vec1;

    EXPECT_EQ(jacact1.size(), 2);
    for( unsigned int i=0; i<2; ++i ) {
      EXPECT_DOUBLE_EQ(jacact1(i), expectedJacact1(i));
    }
  }
}

/// A class to test the behavior of WorkPiece with various input/output types/numbers
class ODETests : public::testing::Test {
public:

  /// Default constructor
  ODETests() {
    // create the right hand side
    rhs = std::make_shared<RHS>();

    // solver options
    pt.put<double>("ODE.RelativeTolerance", 1.0e-10);
    pt.put<double>("ODE.AbsoluteTolerance", 1.0e-10);
    pt.put<double>("ODE.MaxStepSize", 1.0);
    pt.put<unsigned int>("ODE.NumObservations", outTimes.size());
    pt.put<unsigned int>("ODE.MaxNumSteps", 1000);
  }

  /// Default destructor
  virtual ~ODETests() = default;

  virtual void TearDown() override {
    // create the ODE integrator
    auto ode = std::make_shared<ODE>(rhs, pt.get_child("ODE"));

    // integrate the ODE
    const std::vector<Eigen::VectorXd>& result = ode->Evaluate(ic, (Eigen::VectorXd)Eigen::VectorXd::Constant(1, k), outTimes);

    // check the result for the first vector of times
    for( unsigned int t=0; t<outTimes.size(); ++t ) {
      const double time = outTimes(t);

      // check evaluate values
      EXPECT_NEAR(result[0](2*t), ic(0)*std::cos(std::sqrt(k)*time)+ic(1)/std::sqrt(k)*std::sin(std::sqrt(k)*time), 1.0e-6);
      EXPECT_NEAR(result[0](2*t+1), -std::sqrt(k)*ic(0)*std::sin(std::sqrt(k)*time)+ic(1)*std::cos(std::sqrt(k)*time), 1.0e-6);
    }

    // compute jacobians of the initial condition
    const Eigen::MatrixXd& jac00 = ode->Jacobian(0, 0, ic, (Eigen::VectorXd)Eigen::VectorXd::Constant(1, k), outTimes);
    Eigen::MatrixXd jac00true = Eigen::MatrixXd::Zero(2*outTimes.size(), 2);

    EXPECT_EQ(jac00.rows(), 2*outTimes.size()); // rows
    EXPECT_EQ(jac00.cols(), 2); // cols
    for( unsigned int i=0; i<outTimes.size(); ++i ) {
      const double time = outTimes(i);

      jac00true(2*i, 0) = std::cos(std::sqrt(k)*time);
      jac00true(2*i, 1) = std::sin(std::sqrt(k)*time)/std::sqrt(k);
      jac00true(2*i+1, 0) = -std::sqrt(k)*std::sin(std::sqrt(k)*time);
      jac00true(2*i+1, 1) = std::cos(std::sqrt(k)*time);

      // check jacobian wrt initial conditions
      EXPECT_NEAR(jac00(2*i, 0), jac00true(2*i, 0), 1.0e-6);
      EXPECT_NEAR(jac00(2*i, 1), jac00true(2*i, 1), 1.0e-6);
      EXPECT_NEAR(jac00(2*i+1, 0), jac00true(2*i+1, 0), 1.0e-6);
      EXPECT_NEAR(jac00(2*i+1, 1), jac00true(2*i+1, 1), 1.0e-6);
    }

    // compute the gradient wrt initial conditions
    const Eigen::VectorXd sens00 = Eigen::VectorXd::Ones(2*outTimes.size());
    const Eigen::VectorXd& grad00 = ode->Gradient(0, 0, ic, (Eigen::VectorXd)Eigen::VectorXd::Constant(1, k), outTimes, sens00);
    const Eigen::VectorXd expectedGrad00 = jac00true.transpose()*sens00;

    EXPECT_EQ(grad00.size(), expectedGrad00.size());
    EXPECT_EQ(grad00.size(), 2);
    EXPECT_NEAR(grad00(0), expectedGrad00(0), 1.0e-6);
    EXPECT_NEAR(grad00(0), expectedGrad00(0), 1.0e-6);

    // compute jacobians of the parameter k
    const Eigen::MatrixXd& jac01 = ode->Jacobian(0, 1, ic, (Eigen::VectorXd)Eigen::VectorXd::Constant(1, k), outTimes);
    Eigen::MatrixXd jac01true = Eigen::MatrixXd::Zero(2*outTimes.size(), 1);

    EXPECT_EQ(jac01.rows(), 2*outTimes.size()); // rows
    EXPECT_EQ(jac01.cols(), 1); // cols
    for( unsigned int i=0; i<outTimes.size(); ++i ) {
      const double time = outTimes(i);

      jac01true(2*i, 0) = 0.5*(-ic(0)/std::sqrt(k)*time*std::sin(std::sqrt(k)*time)+ic(1)/k*time*std::cos(std::sqrt(k)*time)-ic(1)/std::pow(k, 1.5)*std::sin(std::sqrt(k)*time));
      jac01true(2*i+1, 0) = 0.5*(-ic(0)*time*std::cos(std::sqrt(k)*time)-ic(0)/std::sqrt(k)*std::sin(std::sqrt(k)*time)-ic(1)/std::sqrt(k)*time*std::sin(std::sqrt(k)*time));

      // check jacobian wrt parameter k
      EXPECT_NEAR(jac01(2*i, 0), jac01true(2*i, 0), 1.0e-6);
      EXPECT_NEAR(jac01(2*i+1, 0), jac01true(2*i+1, 0), 1.0e-6);
    }

    // compute the gradient wrt the parameter k
    const Eigen::VectorXd sens01 = Eigen::VectorXd::Ones(2*outTimes.size());
    const Eigen::VectorXd& grad01 = ode->Gradient(0, 1, ic, (Eigen::VectorXd)Eigen::VectorXd::Constant(1, k), outTimes, sens01);
    const Eigen::VectorXd expectedGrad01 = jac01true.transpose()*sens01;

    EXPECT_EQ(grad01.size(), expectedGrad01.size());
    EXPECT_EQ(grad01.size(), 1);
    EXPECT_NEAR(grad01(0), expectedGrad01(0), 1.0e-5);

    // compute jacobians of the times
    const Eigen::MatrixXd& jac02 = ode->Jacobian(0, 2, ic, (Eigen::VectorXd)Eigen::VectorXd::Constant(1, k), outTimes);
    Eigen::MatrixXd jac02true = Eigen::MatrixXd::Zero(2*outTimes.size(), outTimes.size());

    EXPECT_EQ(jac02.rows(), 2*outTimes.size()); // rows
    EXPECT_EQ(jac02.cols(), outTimes.size()); // cols
    for( unsigned int i=0; i<outTimes.size(); ++i ) {
      const double time = outTimes(i);

      jac02true(2*i, i) = -std::sqrt(k)*ic(0)*std::sin(std::sqrt(k)*time)+ic(1)*std::cos(std::sqrt(k)*time);
      jac02true(2*i+1, i) = -k*(ic(0)*std::cos(std::sqrt(k)*time)+ic(1)/std::sqrt(k)*std::sin(std::sqrt(k)*time));

      // check jacobian wrt output times
      EXPECT_NEAR(jac02(2*i, i), jac02true(2*i, i), 1.0e-6);
      EXPECT_NEAR(jac02(2*i+1, i), jac02true(2*i+1, i), 1.0e-6);
    }

    // compute the gradient wrt the output times
    const Eigen::VectorXd sens02 = Eigen::VectorXd::Ones(2*outTimes.size());
    const Eigen::VectorXd& grad02 = ode->Gradient(0, 2, ic, (Eigen::VectorXd)Eigen::VectorXd::Constant(1, k), outTimes, sens02);
    const Eigen::VectorXd expectedGrad02 = jac02true.transpose()*sens02;

    EXPECT_EQ(grad02.size(), expectedGrad02.size());
    EXPECT_EQ(grad02.size(), outTimes.size());
    for( unsigned int i=0; i<outTimes.size(); ++i ) {
      EXPECT_NEAR(grad02(i), expectedGrad02(i), 1.0e-5);
    }
  }

  /// The right hand side
  std::shared_ptr<RHS> rhs;

  /// Options for the ODE integrator
  pt::ptree pt;

  /// The initial condition
  Eigen::VectorXd ic = Eigen::Vector2d(2.0, 0.5);

  /// The spring constant
  const double k = 0.12;

  /// The output times
  const Eigen::VectorXd outTimes = Eigen::VectorXd::LinSpaced(10, 0.0, 1.0);

  private:
};

TEST_F(ODETests, BDFNewtonMethod) {
  pt.put<std::string>("ODE.MultistepMethod", "BDF");
  pt.put<std::string>("ODE.NonlinearSolver", "Newton");
  pt.put<std::string>("ODE.LinearSolver", "Dense");
}

TEST_F(ODETests, BDFIterMethod) {
  pt.put<std::string>("ODE.MultistepMethod", "BDF");
  pt.put<std::string>("ODE.NonlinearSolver", "Iter");
  pt.put<std::string>("ODE.LinearSolver", "Dense");
}

TEST_F(ODETests, AdamsNewtonMethod) {
  pt.put<std::string>("ODE.MultistepMethod", "Adams");
  pt.put<std::string>("ODE.NonlinearSolver", "Newton");
  pt.put<std::string>("ODE.LinearSolver", "Dense");
}

TEST_F(ODETests, AdamsIterMethod) {
  pt.put<std::string>("ODE.MultistepMethod", "Adams");
  pt.put<std::string>("ODE.NonlinearSolver", "Iter");
  pt.put<std::string>("ODE.LinearSolver", "Dense");
}

TEST_F(ODETests, SPGMR) {
  pt.put<std::string>("ODE.MultistepMethod", "BDF");
  pt.put<std::string>("ODE.NonlinearSolver", "Newton");
  pt.put<std::string>("ODE.LinearSolver", "SPGMR");
}

TEST_F(ODETests, SPBCG) {
  pt.put<std::string>("ODE.MultistepMethod", "Adams");
  pt.put<std::string>("ODE.NonlinearSolver", "Newton");
  pt.put<std::string>("ODE.LinearSolver", "SPBCG");
}

TEST_F(ODETests, SPTFQMR) {
  pt.put<std::string>("ODE.MultistepMethod", "Adams");
  pt.put<std::string>("ODE.NonlinearSolver", "Iter");
  pt.put<std::string>("ODE.LinearSolver", "SPTFQMR");
}
