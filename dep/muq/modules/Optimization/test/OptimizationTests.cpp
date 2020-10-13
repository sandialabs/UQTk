#include <gtest/gtest.h>

#include "MUQ/Optimization/NLoptOptimizer.h"

#include "RosenbrockFunction.h"

namespace pt = boost::property_tree;
using namespace muq::Modeling;
using namespace muq::Optimization;

class OptimizationTests : public::testing::Test {
public:

  inline OptimizationTests() {
    pt.put("Optimization.Ftol.AbsoluteTolerance", 1.0e-14);
    pt.put("Optimization.Ftol.RelativeTolerance", 1.0e-14);
    pt.put("Optimization.Xtol.AbsoluteTolerance", 1.0e-14);
    pt.put("Optimization.Xtol.RelativeTolerance", 1.0e-14);
    pt.put("Optimization.MaxEvaluations", 1000); // max number of cost function evaluations
  }

  inline virtual ~OptimizationTests() {}
  
  inline void TearDown() {
    auto cost = std::make_shared<RosenbrockFunction>();
    
    const Eigen::VectorXd x = Eigen::Vector2d(0.85, 1.2);

    auto opt =
      std::make_shared<NLoptOptimizer>(cost, pt.get_child("Optimization"));

    std::vector<Eigen::VectorXd> inputs;
    inputs.push_back(x);

    std::pair<Eigen::VectorXd, double> soln = opt->Solve(inputs);
    
    EXPECT_EQ(soln.first.size(), 2);
    EXPECT_NEAR(soln.first(0), 1.0, 1.0e-4);
    EXPECT_NEAR(soln.first(1), 1.0, 1.0e-4);
    EXPECT_NEAR(soln.second, 0.0, 1.0e-10);

  }
  
  pt::ptree pt;
  
private:
};

TEST_F(OptimizationTests, COBYLA) {
  pt.put("Optimization.Algorithm", "COBYLA");
}

TEST_F(OptimizationTests, BOBYQA) {
  pt.put("Optimization.Algorithm", "BOBYQA");
}

TEST_F(OptimizationTests, NEWUOA) {
  pt.put("Optimization.Algorithm", "NEWUOA");
}

TEST_F(OptimizationTests, PRAXIS) {
  pt.put("Optimization.Algorithm", "PRAXIS");
}

TEST_F(OptimizationTests, NM) {
  pt.put("Optimization.Algorithm", "NM");
}

TEST_F(OptimizationTests, SBPLX) {
  pt.put("Optimization.Algorithm", "SBPLX");
}

TEST_F(OptimizationTests, MMA) {
  pt.put("Optimization.Algorithm", "MMA");
}

TEST_F(OptimizationTests, SLSQP) {
  pt.put("Optimization.Algorithm", "SLSQP");
}

TEST_F(OptimizationTests, LBFGS) {
  pt.put("Optimization.Algorithm", "LBFGS");
}

TEST_F(OptimizationTests, PreTN) {
  pt.put("Optimization.Algorithm", "PreTN");
}

TEST_F(OptimizationTests, LMVM) {
  pt.put("Optimization.Algorithm", "LMVM");
}

class InequalityConstraint : public muq::Modeling::ModPiece {
  public: 
  inline InequalityConstraint() :
    ModPiece(2*Eigen::VectorXi::Ones(1), Eigen::VectorXi::Ones(1)) {}

  virtual inline ~InequalityConstraint() {}

 private:

  inline virtual
  void
  EvaluateImpl(muq::Modeling::ref_vector<Eigen::VectorXd> const& input) override {

    const Eigen::VectorXd& xc = input[0];
    const Eigen::VectorXd b = Eigen::VectorXd::Constant(1, 2.0);

    std::vector<Eigen::VectorXd> output;
    outputs.resize(outputSizes[0]);
    outputs[0] = ((b(0)-xc(0))*Eigen::VectorXd::Ones(1));

  }


  inline virtual 
  void
  JacobianImpl(unsigned int const outputDimWrt,
               unsigned int const inputDimWrt,
               muq::Modeling::ref_vector<Eigen::VectorXd> const& input) override {


    assert(inputDimWrt==0);

    jacobian = Eigen::MatrixXd::Zero(inputSizes[0], outputSizes[0]);
    jacobian(0,0) = -1.0;
    
  }
  
};


class EqualityConstraint : public ModPiece {
  public: 
  inline EqualityConstraint() :
    ModPiece(2*Eigen::VectorXi::Ones(1), Eigen::VectorXi::Ones(1)) {}

  virtual inline ~EqualityConstraint() {}

 private:

  inline virtual
  void EvaluateImpl(muq::Modeling::ref_vector<Eigen::VectorXd> const& input) override {

    const Eigen::VectorXd& xc = input[0];

    outputs.resize(outputSizes[0]);
    outputs[0] = (xc(1)-xc(0)*xc(0)-1.0)*Eigen::VectorXd::Ones(1);

  }


  inline virtual
  void JacobianImpl(unsigned int const outputDimWrt,
                    unsigned int const inputDimWrt,
                    muq::Modeling::ref_vector<Eigen::VectorXd> const& input) override {

    assert(inputDimWrt==0);
    
    const Eigen::VectorXd& xc = input[0];

    jacobian.resize(inputSizes[0], outputSizes[0]);
    jacobian(0,0) = -2.0*xc(0);
    jacobian(1,0) = 1.0;

  }

};



class ConstrainedOptimizationTests : public::testing::Test {
public:

  inline ConstrainedOptimizationTests() {
    pt.put("Optimization.Ftol.AbsoluteTolerance", 1.0e-14);
    pt.put("Optimization.Ftol.RelativeTolerance", 1.0e-14);
    pt.put("Optimization.Xtol.AbsoluteTolerance", 1.0e-14);
    pt.put("Optimization.Xtol.RelativeTolerance", 1.0e-14);
    pt.put("Optimization.ConstraintTolerance", 1.0e-14);
    pt.put("Optimization.MaxEvaluations", 100000); // max number of cost function evaluations
  }

  inline virtual ~ConstrainedOptimizationTests() {}
  
  inline void TearDown() {
    auto cost = std::make_shared<RosenbrockFunction>();
    auto ineqconstraint = std::make_shared<InequalityConstraint>();
    auto eqconstraint = std::make_shared<EqualityConstraint>();
    
    const Eigen::VectorXd x = Eigen::Vector2d(2.0, 5.0);
    const Eigen::VectorXd a = Eigen::VectorXd::Constant(1, 5.0);

	    
    auto opt = std::make_shared<NLoptOptimizer>(cost, pt.get_child("Optimization"));
     
    opt->AddInequalityConstraint(ineqconstraint);
    opt->AddEqualityConstraint(eqconstraint);
    
    std::vector<Eigen::VectorXd> inputs;
    inputs.push_back(x);
    
    std::pair<Eigen::VectorXd, double> soln1 = opt->Solve(inputs);

    EXPECT_EQ(soln1.first.size(), 2);
    EXPECT_NEAR(soln1.first(0), 2.0, 1.0e-4);
    EXPECT_NEAR(soln1.first(1), 5.0, 1.0e-4);
    EXPECT_NEAR(soln1.second, 6.0, 1.0e-10);

  }
  
  pt::ptree pt;
  
private:
};

TEST_F(ConstrainedOptimizationTests, COBYLA) {
  pt.put("Optimization.Algorithm", "COBYLA");
}

TEST_F(ConstrainedOptimizationTests, SLSQP) {
  pt.put("Optimization.Algorithm", "SLSQP");
}
