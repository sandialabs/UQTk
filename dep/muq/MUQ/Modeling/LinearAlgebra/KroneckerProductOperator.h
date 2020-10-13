#ifndef KRONECKERPRODUCTOPERATOR_H
#define KRONECKERPRODUCTOPERATOR_H

#include "MUQ/Modeling/LinearAlgebra/LinearOperator.h"
#include "MUQ/Modeling/LinearAlgebra/SumOperator.h"

namespace muq
{
namespace Modeling
{

    /** @class KroneckerProductOperator
        @ingroup LinearOperators
        @brief Defines the Kronecker product of two other linear operators.
        @details Let \f$A\in\mathbb{R}^M\times\mathbb{R}^N\f$ and \f$B\f$ be two rectangular matrices.  The Kronecker product of \f$A\f$ and \f$B\f$ is given in block form by
\f[
A\otimes B = \left[ \begin{array}{ccc} A_{11} B & \cdots & A_{1N} B\\ \vdots & \ddots & \vdots \\ A_{M1} & \cdots & A_{MN} B \end{array}\right],
\f]
    */
    class KroneckerProductOperator : public LinearOperator
    {

    public:
        KroneckerProductOperator(std::shared_ptr<LinearOperator> Ain,
                                 std::shared_ptr<LinearOperator> Bin);


        /** Fills in the reference \f$y\f$ with \f$y=Ax\f$ */
        virtual Eigen::MatrixXd Apply(Eigen::Ref<const Eigen::MatrixXd> const& x) override;
  
        /** Fill in the reference \f$y\f$ with \f$y = A^Txf$ */
        virtual Eigen::MatrixXd ApplyTranspose(Eigen::Ref<const Eigen::MatrixXd> const& x) override;

    private:
        std::shared_ptr<LinearOperator> A;
        std::shared_ptr<LinearOperator> B;

    }; // class KroneckerProductOperator



    /** Returns the Kronecker sum of \f$A\in\mathbb{R}^{M_A}\times\mathbb{R}^{N_A}\f$ and \f$B\in\mathbb{R}^{M_B}\times\mathbb{R}^{N_B}\f$.  The Kronecker sum is given by
        \f[
           A\otimes I_{N_B} + I_{M_A} \otimes B,
         \f]
where \f$I_N\f$ is the identity matrix of size \f$N\f$.
    */
    std::shared_ptr<LinearOperator> KroneckerSum(std::shared_ptr<LinearOperator> A,
                                                 std::shared_ptr<LinearOperator> B);

    Eigen::MatrixXd KroneckerProduct(Eigen::Ref<const Eigen::MatrixXd> const& A,
                                     Eigen::Ref<const Eigen::MatrixXd> const& B);

} // namespace Modeling
} // namespace muq


#endif // #ifndef KRONECKERPRODUCTOPERATOR_H
