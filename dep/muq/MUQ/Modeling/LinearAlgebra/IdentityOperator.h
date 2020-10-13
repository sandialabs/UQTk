#ifndef IDENTITYOPERATOR_H
#define IDENTITYOPERATOR_H

#include "MUQ/Modeling/LinearAlgebra/LinearOperator.h"

namespace muq{
namespace Modeling{


    class IdentityOperator : public LinearOperator
    {
    public:

        IdentityOperator(const int dimIn) : LinearOperator(dimIn, dimIn){};

        virtual Eigen::MatrixXd Apply(Eigen::Ref<const Eigen::MatrixXd> const& x){return x; }

        virtual Eigen::MatrixXd  ApplyTranspose(Eigen::Ref<const Eigen::MatrixXd> const& x){return x;}

        virtual Eigen::MatrixXd GetMatrix(){return Eigen::MatrixXd::Identity(nrows, ncols);};

    }; // class IdentityOperator

} // namespace Modeling
} // namespace muq

#endif // #ifndef IDENTITYOPERATOR_H
