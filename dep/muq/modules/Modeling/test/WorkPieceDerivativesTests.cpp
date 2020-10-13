#include "gtest/gtest.h"

#include <Eigen/Core>

#include "MUQ/Modeling/WorkGraph.h"
#include "MUQ/Modeling/WorkGraphPiece.h"

using namespace muq::Modeling;

/// A linear function to test the finite difference derivatives
class Linear : public WorkPiece {
public:

  Linear(Eigen::MatrixXd const& Q, Eigen::VectorXd const& b) :
    WorkPiece(std::vector<std::string>({typeid(double).name(), typeid(Eigen::VectorXd).name()}), std::vector<std::string>({typeid(std::string).name(), typeid(Eigen::VectorXd).name()})),
    Q(Q), b(b)
  {}

  virtual ~Linear() {}

private:

  virtual void EvaluateImpl(ref_vector<boost::any> const& inputs) override {
    outputs.resize(2);

    // a double to scale the linear part
    const double a = boost::any_cast<double>(inputs[0]);

    // constant reference to the input vector
    const Eigen::VectorXd& in = boost::any_cast<const Eigen::VectorXd>(inputs[1]);

    // the first output is a string
    outputs[0] = (std::string)"string";

    // compute the linear function
    outputs[1] = (Eigen::VectorXd)(a*Q*in + b);
  }

  /// Matrix for the linear part
  const Eigen::MatrixXd Q;

  /// Vector for the constant part
  const Eigen::VectorXd b;
};

/// A quadratic function to test the derivatives
class Quadratic : public WorkPiece {
public:

  Quadratic(Eigen::MatrixXd const& Q, Eigen::VectorXd const& a, Eigen::VectorXd const& b) :
    WorkPiece(std::vector<std::string>(1, typeid(Eigen::VectorXd).name()), std::vector<std::string>(1, typeid(double).name())),
    Q(Q), a(a), b(b)
  {}

  virtual ~Quadratic() {}

private:

  virtual void EvaluateImpl(ref_vector<boost::any> const& inputs) override {
    outputs.resize(1);

    // constant reference to the input vector
    const Eigen::VectorXd& in = boost::any_cast<Eigen::VectorXd>(inputs[0]);
    // compute the quadratic function
    outputs[0] = (in.transpose()*Q*in + a.transpose()*in + b) (0);
  }

  virtual void JacobianImpl(unsigned int const wrtIn, unsigned int const wrtOut, ref_vector<boost::any> const& inputs) override {
    // constant reference to the input vector
    const Eigen::VectorXd& in = boost::any_cast<const Eigen::VectorXd>(inputs[0]);

    // compute the Jacobian
    jacobian = (Eigen::MatrixXd)(2.0*in.transpose()*Q + a.transpose());
  }

  virtual void JacobianActionImpl(unsigned int const wrtIn, unsigned int const wrtOut, boost::any const& vec, ref_vector<boost::any> const& inputs) override {
    // constant reference to the input vector
    const Eigen::VectorXd& in = boost::any_cast<const Eigen::VectorXd>(inputs[0]);

    // constant reference to the vector we are applying the Jacobian to
    const Eigen::VectorXd& appvec = boost::any_cast<const Eigen::VectorXd&>(vec);

    // compute the Jacobian
    jacobianAction = (2.0*in.transpose()*Q*appvec + a.transpose()*appvec) (0);
  }

  virtual void JacobianTransposeActionImpl(unsigned int const wrtIn, unsigned int const wrtOut, boost::any const& vec, ref_vector<boost::any> const& inputs) override {
    // constant reference to the input vector
    const Eigen::VectorXd& in = boost::any_cast<const Eigen::VectorXd>(inputs[0]);

    // constant reference to the vector we are applying the Jacobian to
    const double appvec = boost::any_cast<const double>(vec);

    // compute the Jacobian
    jacobianTransposeAction = (Eigen::VectorXd)(2.0*Q.transpose()*in*appvec + a*appvec);
  }

  /// Matrix for the quadratic part
  const Eigen::MatrixXd Q;

  /// Vector for the linear part
  const Eigen::VectorXd a;

  /// Vector for the constant part
  const Eigen::VectorXd b;
};

/// A class to test the behavior of WorkPiece with various input/output types/numbers
class WorkPieceDerivativesTests : public::testing::Test {
public:

  /// Default constructor
  WorkPieceDerivativesTests() {
    // a random matrix (for the quadratic term)
    Q = Eigen::MatrixXd::Random(N, N);

    // a random vector (for the linear term)
    a = Eigen::VectorXd::Random(N);

    // a random vector (for the constant term)
    b = Eigen::VectorXd::Random(1);

    // create a linear polynomial
    lin = std::make_shared<Linear>(Q, a);

    // create a quadratic polynomial
    quad = std::make_shared<Quadratic>(Q, a, b);
  }

  /// Default destructor
  virtual ~WorkPieceDerivativesTests() {}

  /// The size of the system
  const unsigned int N = 5;

  /// Matrix for the quadratic part
  Eigen::MatrixXd Q;

  /// Vector for the linear part
  Eigen::VectorXd a;

  /// Vector for the constant part
  Eigen::VectorXd b;

  /// A quadratic function
  std::shared_ptr<Quadratic> quad;

  /// A linear function
  std::shared_ptr<Linear> lin;

private:
};

TEST_F(WorkPieceDerivativesTests, LinearFunction) {
  // check input/output sizes
  EXPECT_EQ(lin->numInputs, 2);
  EXPECT_EQ(lin->numOutputs, 2);

  // choose random inputs
  const Eigen::VectorXd in = Eigen::VectorXd::Random(N);
  const double scalar = 3.5;

  { // test evaluate
    // evaluate
    auto result = lin->Evaluate(scalar, in);

    // make sure we get the expected result
    EXPECT_EQ(result.size(), 2);
    EXPECT_EQ(boost::any_cast<std::string>(result[0]).compare("string"), 0);
    const Eigen::VectorXd& vec = boost::any_cast<Eigen::VectorXd>(result[1]);

    // the expected vector
    const Eigen::VectorXd expectedVec = scalar*Q*in+a;

    EXPECT_EQ(vec.size(), N);
    EXPECT_EQ(expectedVec.size(), N);
    for( unsigned int i=0; i<N; ++i ) {
      EXPECT_DOUBLE_EQ(vec(i), expectedVec(i));
    }
  }

  { // test the jacobian
    // compute the Jacobian and get a reference to it
    auto jac = lin->Jacobian(1, 1, scalar, in);
    const Eigen::MatrixXd& jacref = boost::any_cast<const Eigen::MatrixXd&>(jac);

    EXPECT_EQ(jacref.rows(), N);
    EXPECT_EQ(jacref.cols(), N);
    for( unsigned int i=0; i<N; ++i ) {
      for( unsigned int j=0; j<N; ++j ) {
	// its linear so FD should be exact, but the error is very small ...
	EXPECT_NEAR(jacref(i,j), scalar*Q(i,j), 1.0e-8);
      }
    }
  }

  // a vector to apply the jacobian to
  const Eigen::VectorXd apply = Eigen::VectorXd::Random(N);

  { // test the action of the jacobian
    // compute the action of the jacobian and get a references to it
    auto jacAction = lin->JacobianAction(1, 1, apply, scalar, in);
    const Eigen::VectorXd& jacActionref = boost::any_cast<const Eigen::VectorXd&>(jacAction);

    // compute the exected action of the jacobian
    const Eigen::VectorXd expectedJacAction = scalar*Q*apply;

    EXPECT_EQ(jacActionref.size(), N);
    for( unsigned int i=0; i<N; ++i ) {
      // its linear so FD should be exact, but the error is very small ...
      EXPECT_NEAR(jacActionref(i), expectedJacAction(i), 1.0e-8);
    }
  }

  { // test the action of the jacobian transpose
    // compute the action of the jacobian transpose and get a references to it
    auto jacTransAction = lin->JacobianTransposeAction(1, 1, apply, scalar, in);
    const Eigen::VectorXd& jacTransActionref = boost::any_cast<const Eigen::VectorXd&>(jacTransAction);

    // compute the exected action of the jacobian transpose
    const Eigen::VectorXd expectedJacTransAction = scalar*Q.transpose()*apply;

    EXPECT_EQ(jacTransActionref.size(), N);
    for( unsigned int i=0; i<N; ++i ) {
      // its linear so FD should be exact, but the error is very small ...
      EXPECT_NEAR(jacTransActionref(i), expectedJacTransAction(i), 1.0e-8);
    }
  }
}

TEST_F(WorkPieceDerivativesTests, QuadraticFunction) {
  // check input/output sizes
  EXPECT_EQ(quad->numInputs, 1);
  EXPECT_EQ(quad->numOutputs, 1);

  // choose random inputs
  const Eigen::VectorXd in = Eigen::VectorXd::Random(N);

  { // test evaluate
    // evaluate and make sure we get the expected result
    auto result = quad->Evaluate(in);
    EXPECT_DOUBLE_EQ(boost::any_cast<double>(result[0]), (in.transpose()*Q*in + a.transpose()*in + b) (0));
  }

  { // test jacobian
    // evaluate the jacobian
    auto jacBoost = quad->Jacobian(0, 0, in);
    const Eigen::MatrixXd& jac = boost::any_cast<Eigen::MatrixXd>(jacBoost);

    // compute the expected jacobian
    const Eigen::MatrixXd jacExpected = 2.0*in.transpose()*Q + a.transpose();

    // make sure the jacobians match
    EXPECT_EQ(jac.rows(), 1);
    EXPECT_EQ(jacExpected.rows(), 1);
    EXPECT_EQ(jac.cols(), N);
    EXPECT_EQ(jacExpected.cols(), N);
    for( unsigned int i=0; i<N; ++i ) {
      EXPECT_DOUBLE_EQ(jac(0,i), jacExpected(0,i));
    }
  }

  { // test jacobian action
    // a random vector to apply the Jacobian to
    const Eigen::VectorXd vec = Eigen::VectorXd::Random(N);

    // evaluate the jacobian action
    auto jacAction = quad->JacobianAction(0, 0, vec, in);

    // compute the expected jacobian action
    const double jacActionExpected = (2.0*in.transpose()*Q*vec + a.transpose()*vec) (0);

    // make sure the jacobian action matches
    EXPECT_DOUBLE_EQ(boost::any_cast<double>(jacActionExpected), jacActionExpected);
  }

  { // test jacobian action transpose
    // a random vector to apply the Jacobian transpose to
    const Eigen::VectorXd vec = Eigen::VectorXd::Random(1);

    // evaluate the jacobian action
    auto jacTransposeActionBoost = quad->JacobianTransposeAction(0, 0, vec(0), in);
    const Eigen::VectorXd& jacTransposeAction = boost::any_cast<Eigen::VectorXd>(jacTransposeActionBoost);

    // compute the expected jacobian action
    const Eigen::VectorXd jacTransposeActionExpected = 2.0*Q.transpose()*in*vec + a*vec;

    // make sure the jacobian action matches
    EXPECT_EQ(jacTransposeAction.rows(), N);
    EXPECT_EQ(jacTransposeActionExpected.rows(), N);
    EXPECT_EQ(jacTransposeAction.cols(), 1);
    EXPECT_EQ(jacTransposeActionExpected.cols(), 1);
    for( unsigned int i=0; i<N; ++i ) {
      EXPECT_DOUBLE_EQ(jacTransposeAction(i,0), jacTransposeActionExpected(i,0));
    }
  }
}



/// A 1D polynomial
class Polynomial : public WorkPiece {
public:

  /**
     @param[in] coefficients Coefficients for the polynomial, the first is the constant coefficient, second is linear, ect...
   */
  Polynomial(std::vector<double> const& coefficients) : WorkPiece(std::vector<std::string>({typeid(double).name()}), std::vector<std::string>({typeid(double).name()})), coefficients(coefficients) {
    // needs to be at least constant
    assert(coefficients.size()>0);
  }

  virtual ~Polynomial() {}

private:

  virtual void EvaluateImpl(ref_vector<boost::any> const& inputs) override {
    // get the variable
    const double x = boost::any_cast<double>(inputs[0]);

    // resize the outputs
    outputs.resize(1);

    // compute the polynomial
    outputs[0] = coefficients[0];
    double& result = boost::any_cast<double&>(outputs[0]);
    for( unsigned int i=1; i<coefficients.size(); ++i ) {
      result += coefficients[i]*std::pow(x, (double)i);
    }
  }

  virtual void JacobianImpl(unsigned int const wrtIn, unsigned int const wrtOut, ref_vector<boost::any> const& inputs) override {
    // get the variable
    const double x = boost::any_cast<double>(inputs[0]);

    // compute the polynomial's derivative
    jacobian = 0.0;
    double& jac = boost::any_cast<double&>(*jacobian); // * operator because it is boost::optional
    for( unsigned int i=1; i<coefficients.size(); ++i ) {
      jac += (double)i*coefficients[i]*std::pow(x, (double)i-1.0);
    }
  }

  virtual void JacobianActionImpl(unsigned int const wrtIn, unsigned int const wrtOut, boost::any const& vec, ref_vector<boost::any> const& inputs) override {
    // get the variable
    const double x = boost::any_cast<double>(inputs[0]);

    // get the vector we are applying the jacobian to (scalar since this is a 1D function)
    const double v = boost::any_cast<double>(vec);

    // compute the polynomial's derivative
    jacobianAction = 0.0;
    double& jac = boost::any_cast<double&>(*jacobianAction); // * operator because it is boost::optional
    for( unsigned int i=1; i<coefficients.size(); ++i ) {
      jac += (double)i*coefficients[i]*std::pow(x, (double)i-1.0);
    }

    // multiply by the vector
    jac *= v;
  }

  virtual void JacobianTransposeActionImpl(unsigned int const wrtIn, unsigned int const wrtOut, boost::any const& vec, ref_vector<boost::any> const& inputs) override {
    // get the variable
    const double x = boost::any_cast<double>(inputs[0]);

    // get the vector we are applying the jacobian transpose to (scalar since this is a 1D function)
    const double v = boost::any_cast<double>(vec);

    // compute the polynomial's derivative
    jacobianTransposeAction = 0.0;
    double& jac = boost::any_cast<double&>(*jacobianTransposeAction); // * operator because it is boost::optional
    for( unsigned int i=1; i<coefficients.size(); ++i ) {
      jac += (double)i*coefficients[i]*std::pow(x, (double)i-1.0);
    }

    // multiply by the vector
    jac *= v;
  }

  /// Coefficients for the polynomial, the first is the constant coefficient, second is linear, ect...
  std::vector<double> coefficients;
};

/// A model
class Model : public WorkPiece {
public:

  Model() : WorkPiece(std::vector<std::string>(2, typeid(double).name()), std::vector<std::string>({typeid(Eigen::VectorXd).name(), typeid(double).name()})) {}

  virtual ~Model() {}

private:

  virtual void EvaluateImpl(ref_vector<boost::any> const& inputs) override {
    // get the double inputs
    const double a0 = boost::any_cast<double>(inputs[0]);
    const double a1 = boost::any_cast<double>(inputs[1]);

    // resize the outputs
    outputs.resize(2);

    outputs[0] = (Eigen::VectorXd)(a0*Eigen::VectorXd::LinSpaced(4, 0.0, 1.0));
    outputs[1] = a0*a1;
  }

  virtual void JacobianImpl(unsigned int const wrtIn, unsigned int const wrtOut, ref_vector<boost::any> const& inputs) override {
    // get the double inputs
    const double a0 = boost::any_cast<double>(inputs[0]);
    const double a1 = boost::any_cast<double>(inputs[1]);

    if( wrtIn==0 ) {
      if( wrtOut==0 ) {
	jacobian = (Eigen::MatrixXd)Eigen::VectorXd::LinSpaced(4, 0.0, 1.0);
      } else if( wrtOut==1 ) {
	jacobian = a1;
      }
    } else if( wrtIn==1 ) {
      if( wrtOut==0 ) {
	jacobian = (Eigen::MatrixXd)Eigen::MatrixXd::Zero(4, 1);
      } else if( wrtOut==1 ) {
	jacobian = a0;
      }
    }
  }

  virtual void JacobianActionImpl(unsigned int const wrtIn, unsigned int const wrtOut, boost::any const& vec, ref_vector<boost::any> const& inputs) override {
    // get the double inputs
    const double a0 = boost::any_cast<double>(inputs[0]);
    const double a1 = boost::any_cast<double>(inputs[1]);

    // the vector the jacobian is acting on (the input dimension is always 1 since both inputs are scalar)
    const double v = boost::any_cast<double>(vec);

    if( wrtIn==0 ) {
      if( wrtOut==0 ) {
	jacobianAction = (Eigen::VectorXd)(Eigen::VectorXd::LinSpaced(4, 0.0, 1.0)*v);
      } else if( wrtOut==1 ) {
	jacobianAction = a1*v;
      }
    } else if( wrtIn==1 ) {
      if( wrtOut==0 ) {
	jacobianAction = (Eigen::VectorXd)Eigen::VectorXd::Zero(4);
      } else if( wrtOut==1 ) {
	jacobianAction = a0*v;
      }
    }
  }

  virtual void JacobianTransposeActionImpl(unsigned int const wrtIn, unsigned int const wrtOut, boost::any const& vec, ref_vector<boost::any> const& inputs) override {
    // get the double inputs
    const double a0 = boost::any_cast<double>(inputs[0]);
    const double a1 = boost::any_cast<double>(inputs[1]);

    if( wrtIn==0 ) {
      if( wrtOut==0 ) {
	// the vector the jacobian is acting on (the output dimension is 4)
	const Eigen::VectorXd& v = boost::any_cast<const Eigen::VectorXd&>(vec);
	assert(v.size()==4);

	jacobianTransposeAction = (Eigen::VectorXd::LinSpaced(4, 0.0, 1.0).transpose()*v) (0);
      } else if( wrtOut==1 ) {
	// the vector the jacobian is acting on (the output dimension is 1)
	const double v = boost::any_cast<double>(vec);

	jacobianTransposeAction = a1*v;
      }
    } else if( wrtIn==1 ) {
      if( wrtOut==0 ) {
	jacobianTransposeAction = 0.0;
      } else if( wrtOut==1 ) {
	// the vector the jacobian is acting on (the output dimension is 1)
	const double v = boost::any_cast<double>(vec);

	jacobianTransposeAction = a0*v;
      }
    }
  }
};

TEST(WorkGraphPieceDerivativesTests, GraphDerivatives) {
  // the size of the system
  const unsigned int N = 5;

  // a random matrix (for the quadratic term)
  const Eigen::MatrixXd Q = Eigen::MatrixXd::Random(N, N);

  // a random vector (for the linear term)
  const Eigen::VectorXd a = Eigen::VectorXd::Random(N);

  // a random vector (for the constant term)
  const Eigen::VectorXd b = Eigen::VectorXd::Random(1);

  // create a quadratic polynomial
  const auto quad = std::make_shared<Quadratic>(Q, a, b);

  // create polynomial models
  std::vector<double> coeff0({1.0e-1, 2.0e-2, 3.0e-3, 4.0e-4, 4.0e-5});
  const auto poly0 = std::make_shared<Polynomial>(coeff0);
  std::vector<double> coeff1({4.0e-1, 3.0e-2, 2.0e-3});
  const auto poly1 = std::make_shared<Polynomial>(coeff1);

  // create a model
  const auto mod = std::make_shared<Model>();

  // create a work graph
  auto graph = std::make_shared<WorkGraph>();

  // add the models
  graph->AddNode(quad, "model 0");
  graph->AddNode(poly0, "model 2");
  graph->AddNode(poly1, "model 3");
  graph->AddNode(mod, "model 4");

  // connect the models
  graph->AddEdge("model 0", 0, "model 2", 0);
  graph->AddEdge("model 0", 0, "model 3", 0);
  graph->AddEdge("model 2", 0, "model 4", 0);
  graph->AddEdge("model 3", 0, "model 4", 1);

  // create the model we want to take the derivative of
  const auto graphmod = graph->CreateWorkPiece("model 4");

  // inputs to the model
  const double scalar = 2.5;
  const Eigen::VectorXd invec = Eigen::VectorXd::Random(N);

  // linear
  const Eigen::VectorXd l = scalar*Q*invec+a;

  // quadradic
  const double q = (l.transpose()*Q*l + a.transpose()*l + b) (0);
  const Eigen::MatrixXd q_jac = 2.0*l.transpose()*Q + a.transpose(); // jacobian

    // polynomials
  double d0 = coeff0[0];
  double d0_jac = 0.0; // jacobian
  for( unsigned int i=1; i<coeff0.size(); ++i ) {
    d0 += coeff0[i]*std::pow(q, (double)i);
    d0_jac += (double)i*coeff0[i]*std::pow(q, (double)i-1);
  }

  double d1 = coeff1[0];
  double d1_jac = 0.0; // jacobian
  for( unsigned int i=1; i<coeff1.size(); ++i ) {
    d1 += coeff1[i]*std::pow(q, (double)i);
    d1_jac += (double)i*coeff1[i]*std::pow(q, (double)i-1);
  }

  { // test evaluate
    //const auto result = graphmod->Evaluate(scalar, invec);
    const auto result = graphmod->Evaluate((Eigen::VectorXd)(scalar*Q*invec+a));
    const Eigen::VectorXd& resultVec = boost::any_cast<const Eigen::VectorXd&>(result[0]);

    // expected graphmod result
    const Eigen::VectorXd expectedVec = d0*Eigen::VectorXd::LinSpaced(4, 0.0, 1.0);
    const double expectedDouble = d0*d1;

    // make sure the results match
    EXPECT_EQ(resultVec.size(), 4);
    for( unsigned int i=0; i<4; ++i ) {
      EXPECT_DOUBLE_EQ(resultVec(i), expectedVec(i));
    }
    EXPECT_DOUBLE_EQ(boost::any_cast<double>(result[1]), expectedDouble);
  }

  // expected model Jacobians
  const Eigen::MatrixXd model_jac0 = Eigen::VectorXd::LinSpaced(4, 0.0, 1.0)*d0_jac*q_jac;
  const Eigen::MatrixXd model_jac1 = (d0*d1_jac + d1*d0_jac)*q_jac;

  { // test Jacobian
    // check first Jacobian
    const auto jacBoost0 = graphmod->Jacobian(0, 0, l);
    const Eigen::MatrixXd& jac0 = boost::any_cast<const Eigen::MatrixXd&>(jacBoost0);

    EXPECT_EQ(jac0.rows(), 4);
    EXPECT_EQ(jac0.cols(), N);
    EXPECT_EQ(model_jac0.rows(), 4);
    EXPECT_EQ(model_jac0.cols(), N);
    for( unsigned int i=0; i<4; ++i ) {
      for( unsigned int j=0; j<N; ++j ) {
	EXPECT_DOUBLE_EQ(jac0(i,j), model_jac0(i,j));
      }
    }

    // check second Jacobian
    const auto jacBoost1 = graphmod->Jacobian(0, 1, l);
    const Eigen::MatrixXd& jac1 = boost::any_cast<const Eigen::MatrixXd&>(jacBoost1);

    EXPECT_EQ(jac1.rows(), 1);
    EXPECT_EQ(jac1.cols(), N);
    EXPECT_EQ(model_jac1.rows(), 1);
    EXPECT_EQ(model_jac1.cols(), N);
    for( unsigned int i=0; i<N; ++i ) {
      EXPECT_DOUBLE_EQ(jac1(0,i), model_jac1(0,i));
    }
  }

  { // test JacobianAction
    // apply the Jacobians to this vector
    const Eigen::VectorXd vec = Eigen::VectorXd::Random(N);

    // check first JacobianAction
    const auto jacActionBoost0 = graphmod->JacobianAction(0, 0, vec, l);
    const Eigen::VectorXd& jacAction0 = boost::any_cast<const Eigen::VectorXd&>(jacActionBoost0);

    const Eigen::VectorXd model_jacAction0 = model_jac0*vec;

    EXPECT_EQ(jacAction0.size(), 4);
    EXPECT_EQ(model_jacAction0.size(), 4);
    for( unsigned int i=0; i<4; ++i ) {
      EXPECT_NEAR(jacAction0(i), model_jacAction0(i), 1.0e-12);
    }

    // check second JacobianAction
    const auto jacActionBoost1 = graphmod->JacobianAction(0, 1, vec, l);
    const double jacAction1 = boost::any_cast<const double>(jacActionBoost1);

    const Eigen::VectorXd model_jacAction1 = model_jac1*vec;

    EXPECT_EQ(model_jacAction1.size(), 1);
    EXPECT_NEAR(jacAction1, model_jacAction1(0), 1.0e-10);
  }

    { // test JacobianAction
    // apply the Jacobian tranpose to this vector
    const Eigen::VectorXd vec0 = Eigen::VectorXd::Random(4);

    // check first JacobianTransposeAction
    const auto jacTransActionBoost0 = graphmod->JacobianTransposeAction(0, 0, vec0, l);
    const Eigen::VectorXd& jacTransAction0 = boost::any_cast<const Eigen::VectorXd&>(jacTransActionBoost0);

    const Eigen::VectorXd model_jacTransAction0 = model_jac0.transpose()*vec0;

    EXPECT_EQ(jacTransAction0.size(), N);
    EXPECT_EQ(model_jacTransAction0.size(), N);
    for( unsigned int i=0; i<N; ++i ) {
      EXPECT_NEAR(jacTransAction0(i), model_jacTransAction0(i), 1.0e-12);
    }

    // apply the Jacobian tranpose to this vector
    const Eigen::VectorXd vec1 = Eigen::VectorXd::Random(1);

    // check second JacobianTransposeAction
    const auto jacTransActionBoost1 = graphmod->JacobianTransposeAction(0, 1, vec1(0), l);
    const Eigen::VectorXd& jacTransAction1 = boost::any_cast<const Eigen::VectorXd&>(jacTransActionBoost1);

    const Eigen::VectorXd model_jacTransAction1 = model_jac1.transpose()*vec1;

    EXPECT_EQ(jacTransAction1.size(), N);
    EXPECT_EQ(model_jacTransAction1.size(), N);
    for( unsigned int i=0; i<N; ++i ) {
      EXPECT_NEAR(jacTransAction1(i), model_jacTransAction1(i), 1.0e-12);
    }
  }
}
