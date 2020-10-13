#ifndef ZEROOPERATOR_H
#define ZEROOPERATOR_H

#include "MUQ/Modeling/LinearAlgebra/LinearOperator.h"

namespace muq{
namespace Modeling{

  class ZeroOperator : public LinearOperator
  {
  public:
      ZeroOperator(int rowsIn, int colsIn) : LinearOperator(rowsIn, colsIn){};

      /** Apply the linear operator to a vector */
      virtual Eigen::MatrixXd Apply(Eigen::Ref<const Eigen::MatrixXd> const& x) override
      {
          return Eigen::MatrixXd::Zero(nrows, x.cols());
      };
      
      /** Apply the transpose of the linear operator to a vector. */
      virtual Eigen::MatrixXd ApplyTranspose(Eigen::Ref<const Eigen::MatrixXd> const& x) override
      {
          return Eigen::MatrixXd::Zero(ncols, x.cols());
      };
  };

} // namespace muq
} // namespace Modeling



#endif // ZEROOPERATOR
