#include "MUQ/Optimization/CostFunction.h"

/// A muq::Modeling::ModPiece version
class RosenbrockModPiece : public muq::Modeling::ModPiece {
public:
  inline RosenbrockModPiece() : muq::Modeling::ModPiece(Eigen::Vector2i(2, 1), Eigen::VectorXi::Constant(1, 1)) {}

  virtual inline ~RosenbrockModPiece() {}

 private:

  inline virtual void EvaluateImpl(muq::Modeling::ref_vector<Eigen::VectorXd> const& input) override {
    const Eigen::VectorXd& xc = input[0];
    const Eigen::VectorXd& a = input[1];

    outputs.resize(1);
    outputs[0] = (Eigen::VectorXd)Eigen::VectorXd::Constant(1, (1.0-xc(0))*(1.0-xc(0))+a(0)*(xc(1)-xc(0)*xc(0))*(xc(1)-xc(0)*xc(0)));
  }

  inline virtual void GradientImpl(unsigned int const outputDimWrt, unsigned int const inputDimWrt, muq::Modeling::ref_vector<Eigen::VectorXd> const& input, Eigen::VectorXd const& sensitivity) override {
    assert(outputDimWrt==0);
    assert(inputDimWrt==0);

    const Eigen::VectorXd& xc = input[0];
    const Eigen::VectorXd& a = input[1];

    gradient = Eigen::Vector2d::Constant(2, std::numeric_limits<double>::quiet_NaN());
    gradient(0) = -4.0*a(0)*(xc(1)-xc(0)*xc(0))*xc(0)-2.0*(1.0-xc(0));
    gradient(1) = 2.0*a(0)*(xc(1)-xc(0)*xc(0));

    gradient *= sensitivity(0);
  }
};

/// A muq::Optimization::CostFunction version
class RosenbrockFunction : public muq::Optimization::CostFunction {
public:
  inline RosenbrockFunction() : CostFunction(2*Eigen::VectorXi::Ones(1)) {}

  virtual inline ~RosenbrockFunction() {}

 private:

  inline virtual double CostImpl(muq::Modeling::ref_vector<Eigen::VectorXd> const& input) override {

    const Eigen::VectorXd& xc = input[0];
    const Eigen::VectorXd a = Eigen::VectorXd::Constant(1, 5.0);

    return (1.0-xc(0))*(1.0-xc(0))+a(0)*(xc(1)-xc(0)*xc(0))*(xc(1)-xc(0)*xc(0));

  }

  inline virtual
   void GradientImpl(unsigned int const inputDimWrt,
                     muq::Modeling::ref_vector<Eigen::VectorXd> const& input,
                     Eigen::VectorXd const& sensitivity) override {

    assert(inputDimWrt==0);

    const Eigen::VectorXd& xc = input[0];
    const Eigen::VectorXd a = Eigen::VectorXd::Constant(1, 5.0);

    gradient = Eigen::Vector2d::Constant(2, std::numeric_limits<double>::quiet_NaN());
    gradient(0) = -4.0*a(0)*(xc(1)-xc(0)*xc(0))*xc(0)-2.0*(1.0-xc(0));
    gradient(1) = 2.0*a(0)*(xc(1)-xc(0)*xc(0));

    gradient *= sensitivity(0);
  }
};
