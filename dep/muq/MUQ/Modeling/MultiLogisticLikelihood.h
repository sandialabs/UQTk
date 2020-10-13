#ifndef MULTILOGISTICLIKELIHOOD_H
#define MULTILOGISTICLIKELIHOOD_H

#include "MUQ/Modeling/ModPiece.h"

namespace muq{
  namespace Modeling{

    /** @class MultiLogisticLikelihood
        @brief Class that defines the likelihood function for multinomial logistic regression
        @details Consider a discrete random variable \f$y\in\{1,2,\ldots, M\}\f$
        that can take one of \f$M\f$ values.  Assume this random depends on another
        indpendent variable \f$x\in\mathbb{R}^N\f$ and we are interested in modeling
        the probabilities \f$p_i(x) = \mathbb{P}[y=i]\f$, which are functions of
        the independent varible \f$x\f$.  Logistic regression is the classic
        approach when \f$M=2\f$ and Multinomial Logistic Regression is a natural
        way to extend logistic regression to the case where \f$M>2\f$.

        The basic idea is to model the log probabilities \f$\log(p_i)\f$ through
        an expression of the form
        \f[
        \log(p_i) = f_i(x) - \log(z(x)),
        \f]
        where \f$\log(z(x))\f$ is a normalization term that ensure the sum of
        the probabilities is equal to \f$1\f$.  More precisely,
        \f[
        z = \sum_{i=1}^M \exp(f_i(x)).
        \f]
        Typically, the functions \f$f_i\f$ would also be represent through an
        expansion of the form
        \f[
        \log(p_i) = -\log(z) + \sum_{j=1}^P \theta_j \phi_j(x),
        \f]
        for some basis functions \f$\phi_j\f$.

        This class evaluates \f$\log(p_i)\f$ for each \f$i\f$ given the values
        of \f$f_i(x)\f$ and then computes the likelihood of \f$K\f$ observations
        \f$y_1,\ldots, y_K\f$.  The data is passed to the constructor.

        This piece assumes that each observation \f$y_k\f$ corresponds to a condition
        \f$x_k\f$.  The input to this function is thus an unrolled vector
        containing \f$f_i(x_k)\f$ for all \f$i\in\left{1,\ldots,M\right}\f$ and
        \f$k\in\{1,\ldots, K\}\f$.  It is assumed that \f$i\f$ is the faster index,
        so that the input vector would look like \f$[f_1(x_1),\ldots, f_M(x_1),f_1(x_2),\ldots]\f$.

        This piece returns a scalar value, with the value of the loglikelihood
        \f[
        \log \pi(y| x) = \sum_{k=1}^K \log(p_{y_k}(x_k))
        \f]
    */
    class MultiLogisticLikelihood : public ModPiece {
    public:

      MultiLogisticLikelihood(unsigned int numClasses, Eigen::VectorXi const& data);

    private:

      virtual void EvaluateImpl(muq::Modeling::ref_vector<Eigen::VectorXd> const& inputs) override;

      virtual void GradientImpl(unsigned int                const  outputDimWrt,
                                unsigned int                const  inputDimWrt,
                                ref_vector<Eigen::VectorXd> const& input,
                                Eigen::VectorXd             const& sensitivity) override;

      virtual void JacobianImpl(unsigned int                const  outputDimWrt,
                                unsigned int                const  inputDimWrt,
                                ref_vector<Eigen::VectorXd> const& input) override;

      virtual void ApplyJacobianImpl(unsigned int                const  outputDimWrt,
                                     unsigned int                const  inputDimWrt,
                                     ref_vector<Eigen::VectorXd> const& input,
                                     Eigen::VectorXd             const& vec) override;

      const unsigned int numClasses;
      Eigen::VectorXi data;
    };

  } // namespace Modeling
} // namespace muq


#endif
