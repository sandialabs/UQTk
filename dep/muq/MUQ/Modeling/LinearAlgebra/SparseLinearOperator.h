#ifndef SPARSELINEAROPERATOR_H
#define SPARSELINEAROPERATOR_H

#include "MUQ/Modeling/LinearAlgebra/LinearOperator.h"
#include <Eigen/SparseCore>

namespace muq
{
namespace Modeling
{

    
/** @class SparseLinearOperator
 *  @ingroup Utilities
 *  @brief Wraps a sparse Eigen matrix into a linear operator
 */
class SparseLinearOperator : public LinearOperator {
public:

  SparseLinearOperator(Eigen::SparseMatrix const& Ain) : LinearOperator(A.rows(), A.cols()), A(Ain){}

  virtual ~SparseLinearOperator(){};
  
  /** Apply the linear operator to a vector */
  virtual Eigen::MatrixXd Apply(Eigen::Ref<Eigen::MatrixXd> const& x) override {return A*x;};

  /** Apply the transpose of the linear operator to a vector. */
  virtual Eigen::MatrixXd ApplyTranspose(Eigen::Ref<Eigen::MatrixXd> const& x) override {return A.transpose()*x};

  /** Fills in the reference \f$y\f$ with \f$y=Ax\f$ */
  virtual void Apply(Eigen::Ref<Eigen::MatrixXd> const& x, Eigen::Ref<Eigen::MatrixXd> y) override { y = A*x;};

  /** Fill in the reference \f$y\f$ with \f$y = A^Txf$ */
  virtual void ApplyTranspose(Eigen::Ref<Eigen::MatrixXd> const& x, Eigen::Ref<Eigen::MatrixXd> y) override {y = A.transpose()*x;};

protected:
  Eigen::MatrixXd A;
  
};

} // namespace Modeling
} // namespace MUQ



#endif
