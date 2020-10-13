#include <gtest/gtest.h>

#include "MUQ/Approximation/Regression/Regression.h"

namespace pt = boost::property_tree;
using namespace muq::Approximation;

class RegressionTest : public::testing::Test {
public:
  inline RegressionTest() {
    const unsigned int Npts = 14;

    // generate the input points
    ins.resize(Npts, Eigen::VectorXd::Constant(2, std::numeric_limits<double>::quiet_NaN()));
    for( auto it=ins.begin(); it!=ins.end(); ++it ) { *it = Eigen::Vector2d::Random(); }

    // generate the output points
    outs.resize(ins.size());
    for( std::vector<Eigen::VectorXd>::size_type i=0; i<ins.size(); ++i ) { outs[i] = f(ins[i]); }

    pt.put<unsigned int>("MyRegression.Order", order);
    pt.put<unsigned int>("MyRegression.InputSize", 2);
  }

  inline virtual ~RegressionTest() {}

  inline Eigen::VectorXd f(Eigen::VectorXd const& x) const {
    Eigen::VectorXd out = Eigen::VectorXd::Constant(2, std::numeric_limits<double>::quiet_NaN());
    out(0) = x(1)*x(1)*x(1)+x(0)*x(1)-x(0);
    out(1) = x(0)*x(1)*x(1)+x(0)*x(0)+1.5;

    return out;
  }

  inline virtual void TearDown() override {
    // create the regression
    auto reg = std::make_shared<Regression>(pt.get_child("MyRegression"));

    // fit the polynomial coefficients
    reg->Fit(ins, outs);

    // points to test the evaluate
    const Eigen::VectorXd x = 0.1*Eigen::VectorXd::Ones(2);//Eigen::Vector2d::Random();
    const Eigen::VectorXd y = 0.2*Eigen::VectorXd::Ones(2);//Eigen::Vector2d::Random();
    const Eigen::VectorXd z = 0.5*Eigen::VectorXd::Ones(2);//Eigen::Vector2d::Random();

    // evaluate the polynomial
    const std::vector<boost::any>& output = reg->Evaluate(x, y, z);
    const Eigen::MatrixXd& result = boost::any_cast<Eigen::MatrixXd const&>(output[0]);
    EXPECT_EQ(result.rows(), 2);
    EXPECT_EQ(result.cols(), 3);

    // compute the true function values---should be exact (we are using 3rd order to estimate a degree 3 polynomial)
    const Eigen::Vector2d x_true = f(x);
    const Eigen::Vector2d y_true = f(y);
    const Eigen::Vector2d z_true = f(z);

    EXPECT_NEAR((x_true-result.col(0)).norm(), 0.0, 1.0e-9);
    EXPECT_NEAR((y_true-result.col(1)).norm(), 0.0, 1.0e-9);
    EXPECT_NEAR((z_true-result.col(2)).norm(), 0.0, 1.0e-9);

    std::pair<Eigen::VectorXd, double> lambda = reg->PoisednessConstant(ins, x);
    unsigned int count = 0;
    while( lambda.second>25.0 && count++<1000 ) { // adding the computed point should improve bad poisedness
      ins.push_back(lambda.first);
      lambda = reg->PoisednessConstant(ins, x);
    }
    EXPECT_TRUE(lambda.second<25.0);
    EXPECT_TRUE(count<1000);
  }

  /// A matrix holding the input points.
  std::vector<Eigen::VectorXd> ins;

  /// A matrix holding the output points.
  std::vector<Eigen::VectorXd> outs;

  /// The order of the polynomial regression
  const unsigned int order = 3;

  /// The options for the regression
  pt::ptree pt;
};

TEST_F(RegressionTest, LegendreBasis) {} // should use Legendre basis by default

TEST_F(RegressionTest, MonomialBasis) {
  pt.put<std::string>("MyRegression.PolynomialBasis", "Monomial");
}

TEST_F(RegressionTest, PhysicistHermiteBasis) {
  pt.put<std::string>("MyRegression.PolynomialBasis", "PhysicistHermite");
}

TEST_F(RegressionTest, ProbabilistHermiteBasis) {
  pt.put<std::string>("MyRegression.PolynomialBasis", "ProbabilistHermite");
}

TEST_F(RegressionTest, LaguerreBasis) {
  pt.put<std::string>("MyRegression.PolynomialBasis", "Laguerre");
}

TEST_F(RegressionTest, JacobiBasis) {
  pt.put<std::string>("MyRegression.PolynomialBasis", "Jacobi");
}
