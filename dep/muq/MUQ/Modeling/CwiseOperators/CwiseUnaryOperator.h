#ifndef CWISEUNARYOPERATOR_H
#define CWISEUNARYOPERATOR_H

#include <Eigen/Core>
#include <stan/math/fwd/scal.hpp>

#include "MUQ/Modeling/ModPiece.h"

namespace muq{
  namespace Modeling{

    template<double (*T1)(double), stan::math::fvar<double> (*T2)(stan::math::fvar<double> const&)>
    class CwiseUnaryOperator : public ModPiece{

    public:
      CwiseUnaryOperator(unsigned int dim) : ModPiece(dim*Eigen::VectorXi::Ones(1), dim*Eigen::VectorXi::Ones(1)){};

    private:
      virtual void EvaluateImpl(ref_vector<Eigen::VectorXd> const& in) override
      {
        outputs.resize(1);
        outputs.at(0).resize(inputSizes(0));
        for(int i=0; i<inputSizes(0); ++i)
          outputs.at(0)(i) = T1(in.at(0).get()(i));
      };

      virtual void JacobianImpl(unsigned int outWrt,
                                unsigned int inWrt,
                                ref_vector<Eigen::VectorXd> const& in) override
      {
        jacobian = Eigen::MatrixXd::Zero(inputSizes(0), inputSizes(0));
        for(int i=0; i<inputSizes(0); ++i){
          stan::math::fvar<double> x(in.at(0).get()(i), 1.0);
          jacobian(i,i) = T2(x).tangent();
        }
      };

      virtual void GradientImpl(unsigned int outWrt,
                        unsigned int inWrt,
                        ref_vector<Eigen::VectorXd> const& in,
                        Eigen::VectorXd const& sens) override
      {
        gradient.resize(inputSizes(0));
        for(int i=0; i<inputSizes(0); ++i){
          stan::math::fvar<double> x(in.at(0).get()(i), 1.0);
          gradient(i) = sens(i)*T2(x).tangent();
        }
      };

      virtual void ApplyJacobianImpl(unsigned int outWrt,
                        unsigned int inWrt,
                        ref_vector<Eigen::VectorXd> const& in,
                        Eigen::VectorXd const& vec) override
      {
        jacobianAction.resize(inputSizes(0));
        for(int i=0; i<inputSizes(0); ++i){
          stan::math::fvar<double> x(in.at(0).get()(i), 1.0);
          jacobianAction(i) = vec(i)*T2(x).tangent();
        }
      };


    };

    // Define specific operators
    typedef CwiseUnaryOperator<std::exp, stan::math::exp> ExpOperator;
    typedef CwiseUnaryOperator<std::cos, stan::math::cos> CosOperator;
    typedef CwiseUnaryOperator<std::sin, stan::math::sin> SinOperator;
    typedef CwiseUnaryOperator<std::abs, stan::math::abs> AbsOperator;
    typedef CwiseUnaryOperator<std::acos, stan::math::acos> AcosOperator;
    typedef CwiseUnaryOperator<std::asin, stan::math::asin> AsinOperator;
    typedef CwiseUnaryOperator<std::atan, stan::math::atan> AtanOperator;
    typedef CwiseUnaryOperator<std::atanh, stan::math::atanh> AtanhOperator;
    typedef CwiseUnaryOperator<std::cbrt, stan::math::cbrt> CbrtOperator;
    typedef CwiseUnaryOperator<std::ceil, stan::math::ceil> CeilOperator;
    typedef CwiseUnaryOperator<std::cosh, stan::math::cosh> CoshOperator;
    typedef CwiseUnaryOperator<stan::math::digamma, stan::math::digamma> DigammaOperator;
    typedef CwiseUnaryOperator<std::erf, stan::math::erf> ErfOperator;
    typedef CwiseUnaryOperator<std::erfc, stan::math::erfc> ErfcOperator;
    typedef CwiseUnaryOperator<std::floor, stan::math::floor> FloorOperator;
    typedef CwiseUnaryOperator<stan::math::inv_logit, stan::math::inv_logit> InvLogitOperator;
    typedef CwiseUnaryOperator<stan::math::inv_Phi, stan::math::inv_Phi> InvPhiOperator;
    typedef CwiseUnaryOperator<stan::math::inv_sqrt, stan::math::inv_sqrt> InvSqrtOperator;
    typedef CwiseUnaryOperator<stan::math::inv_square, stan::math::inv_square> InvSquareOperator;
    typedef CwiseUnaryOperator<stan::math::inv, stan::math::inv> InvOperator;
    typedef CwiseUnaryOperator<stan::math::lgamma, stan::math::lgamma> LogGammaOperator;
    typedef CwiseUnaryOperator<stan::math::log_inv_logit, stan::math::log_inv_logit> LogInvLogitOperator;
    typedef CwiseUnaryOperator<std::log, stan::math::log> LogOperator;
    typedef CwiseUnaryOperator<std::log2, stan::math::log2> Log2Operator;
    typedef CwiseUnaryOperator<std::log10, stan::math::log10> Log10Operator;
    typedef CwiseUnaryOperator<stan::math::logit, stan::math::logit> LogitOperator;
    typedef CwiseUnaryOperator<stan::math::Phi, stan::math::Phi> PhiOperator;
    typedef CwiseUnaryOperator<std::round, stan::math::round> RoundOperator;
    typedef CwiseUnaryOperator<std::sinh, stan::math::sinh> SinhOperator;
    typedef CwiseUnaryOperator<std::sqrt, stan::math::sqrt> SqrtOperator;
    typedef CwiseUnaryOperator<stan::math::square, stan::math::square> SquareOperator;
    typedef CwiseUnaryOperator<std::tan, stan::math::tan> TanOperator;
    typedef CwiseUnaryOperator<std::tanh, stan::math::tanh> TanhOperator;
    typedef CwiseUnaryOperator<std::tgamma, stan::math::tgamma> TgammaOperator;
    typedef CwiseUnaryOperator<stan::math::trigamma, stan::math::trigamma> TrigammaOperator;

  }
}

#endif
