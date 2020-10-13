#ifndef AFFINEOPERATOR_H
#define AFFINEOPERATOR_H

#include "MUQ/Modeling/ModPiece.h"
#include "MUQ/Modeling/LinearAlgebra/LinearOperator.h"

#include <Eigen/Core>

#include <memory>
#include <iostream>
#include <exception>

namespace muq
{
namespace Modeling
{

/** @class AffineOperator
 *  @ingroup LinearOperators
 *  @brief Generic affine operator which adds an offset to a linear operator.
    @seealso LinearOperator
 *  @details Builds on the LinearOperator class to add offsets and define affine operators.
 @code{.cpp}

#include "MUQ/Modeling/LinearAlgebra/AffineOperator.h"

// ...

Eigen::MatrixXd A = Eigen::MatrixXd::Random(10,10);
Eigen::VectorXd B = Eigen::VectorXd::Random(10);
std::shared_ptr<muq::Modeling::AffineOperator> Aop = muq::Modeling::AffineOperator::Create(A,b);

@endcode
 */
class AffineOperator : public muq::Modeling::ModPiece{
public:

  template<typename T>
  AffineOperator(T const& Ain, Eigen::VectorXd const& bIn) : AffineOperator(LinearOperator::Create(Ain), bIn) {}

  AffineOperator(std::shared_ptr<LinearOperator> const& Ain, Eigen::VectorXd const& bIn);

  virtual ~AffineOperator(){};

  /** Returns the linear operator part of this affine operator. */
  std::shared_ptr<LinearOperator> Linear() const{return A;};

  /** Returns the offset part of this affine operator. */
  Eigen::VectorXd const& Offset() const{return b;};

  /** The output dimension of the linear operator. */
  int rows() const { return A->rows(); }

  /** The input dimension of the linear operator. */
  int cols() const { return A->cols(); }

  template<typename OtherType>
  static std::shared_ptr<AffineOperator> Create(OtherType const& A, Eigen::VectorXd const& b)
  {
      return std::make_shared<AffineOperator>(LinearOperatorFactory<OtherType>::Create(A),b);
  }

protected:

  std::shared_ptr<LinearOperator> A;
  Eigen::VectorXd b;

  // ModPiece overrides
  virtual void EvaluateImpl(muq::Modeling::ref_vector<Eigen::VectorXd> const& input) override;

  virtual void GradientImpl(unsigned int                const  outputDimWrt,
                            unsigned int                const  inputDimWrt,
                            muq::Modeling::ref_vector<Eigen::VectorXd> const& input,
                            Eigen::VectorXd             const& sensitivity) override;

  virtual void JacobianImpl(unsigned int                const  outputDimWrt,
                            unsigned int                const  inputDimWrt,
                            muq::Modeling::ref_vector<Eigen::VectorXd> const& input) override;

  virtual void ApplyJacobianImpl(unsigned int                const  outputDimWrt,
                                 unsigned int                const  inputDimWrt,
                                 muq::Modeling::ref_vector<Eigen::VectorXd> const& input,
                                 Eigen::VectorXd             const& vec) override;
};




} // namespace Modeling
} // namespace MUQ

#endif
