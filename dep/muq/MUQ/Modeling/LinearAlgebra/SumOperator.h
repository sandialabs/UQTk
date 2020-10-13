#ifndef SUMOPERATOR_H
#define SUMOPERATOR_H

#include "MUQ/Modeling/LinearAlgebra/LinearOperator.h"


namespace muq
{
namespace Modeling
{

    class SumOperator : public LinearOperator
    {
    public:
        SumOperator(std::shared_ptr<LinearOperator> Ain,
                    std::shared_ptr<LinearOperator> Bin);

        virtual ~SumOperator(){};
        
        /** Apply the linear operator to a vector */
        virtual Eigen::MatrixXd Apply(Eigen::Ref<const Eigen::MatrixXd> const& x) override;
  
        /** Apply the transpose of the linear operator to a vector. */
        virtual Eigen::MatrixXd ApplyTranspose(Eigen::Ref<const Eigen::MatrixXd> const& x) override;

        virtual Eigen::MatrixXd GetMatrix() override{return A->GetMatrix() + B->GetMatrix();};
        
    private:
        std::shared_ptr<LinearOperator> A, B;

    }; // class SumOperator


} // namespace Modeling
} // namespace muq

#endif // #ifndef SUMOPERATOR_H
