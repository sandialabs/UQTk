#ifndef COMPANIONMATRIX_H
#define COMPANIONMATRIX_H

#include "MUQ/Modeling/LinearAlgebra/LinearOperator.h"
#include <Eigen/Core>

namespace muq
{
namespace Modeling
{

    /** @class CompanionMatrix
        @ingroup Utilities
        @brief Implments a companion matrix as a linear operator.
        @details The companion matrix of a monomial of order \f$m-1\f$ takes the form
\f[
A = \left[\begin{array}{cccc} 0 & 1 & & \\ & \ddots & \ddots & \\ & & 0 & 1\\ -a_0 & \cdots & -a_{m-2} & -a_{m-1} \end{array}\right],
\f]
where \f$a_i\f$ are the coefficients of the polynomial.
    */
    class CompanionMatrix : public LinearOperator
    {
    public:
        CompanionMatrix(Eigen::VectorXd const& lastRowIn) : LinearOperator(lastRowIn.size(), lastRowIn.size()), lastRow(lastRowIn){};
        
        virtual Eigen::MatrixXd Apply(Eigen::Ref<const Eigen::MatrixXd> const& x) override;

        virtual Eigen::MatrixXd ApplyTranspose(Eigen::Ref<const Eigen::MatrixXd> const& x) override;

        virtual Eigen::MatrixXd GetMatrix() override;
        
    private:
        Eigen::VectorXd lastRow;
    };

} // namespace Modeling
} // namespace muq

#endif // #ifndef COMPANIONMATRIX_H
