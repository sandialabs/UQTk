#ifndef EIGENLINEAROPERATOR_H
#define EIGENLINEAROPERATOR_H

#include "MUQ/Modeling/LinearAlgebra/LinearOperator.h"

#include <Eigen/SparseCore>

#include <memory>

namespace muq
{
namespace Modeling
{

    
/** @class DenseLinearOperator
 *  @ingroup Utilities
 *  @brief Wraps a general Eigen::MatrixXd into a linear operator
 */
template<typename EigenType>
class EigenLinearOperator : public LinearOperator {
public:

EigenLinearOperator(EigenType const& Ain) : LinearOperator(Ain.rows(), Ain.cols()), A(Ain){}

  virtual ~EigenLinearOperator(){};
  
  /** Apply the linear operator to a vector */
  virtual Eigen::MatrixXd Apply(Eigen::Ref<const Eigen::MatrixXd> const& x) override {return A*x;};

  /** Apply the transpose of the linear operator to a vector. */
  virtual Eigen::MatrixXd ApplyTranspose(Eigen::Ref<const Eigen::MatrixXd> const& x) override {return A.transpose()*x;};

  virtual Eigen::MatrixXd GetMatrix() override{ return Eigen::MatrixXd(A);};
  
protected:
  EigenType A;
  
};

//template<typename Derived>
//std::shared_ptr<LinearOperator> LinearOperator::Create<Eigen::MatrixBase<Derived>>(Eigen::MatrixBase<Derived> const& A)
//{
//   return std::make_shared<EigenLinearOperator>(A);
//}

template<typename ScalarType>
struct LinearOperatorFactory<Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic>> {
    
  static std::shared_ptr<LinearOperator> Create(Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic> const& A)
  {
      return std::make_shared<muq::Modeling::EigenLinearOperator<Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic>>>(A);
  }
};

template<typename ScalarType>
struct LinearOperatorFactory<Eigen::SparseMatrix<ScalarType>> {
    
  static std::shared_ptr<LinearOperator> Create(Eigen::SparseMatrix<ScalarType> const& A)
  {
    return std::make_shared<muq::Modeling::EigenLinearOperator<Eigen::SparseMatrix<ScalarType>>>(A);
  }
};

    
//template<typename Derived>
//std::shared_ptr<LinearOperator> muq::Modeling::LinearOperatorFactory<Eigen::MatrixBase<Derived>>::Create(Eigen::MatrixBase<Derived> const& A)
//{
//    return std::make_shared<muq::Modeling::EigenLinearOperator<Eigen::MatrixBase<Derived>>>(A);
//}

} // namespace Modeling
} // namespace MUQ



#endif
