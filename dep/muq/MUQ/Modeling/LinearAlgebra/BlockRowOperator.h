#ifndef BLOCKROWOPERATOR_H
#define BLOCKROWOPERATOR_H

#include "MUQ/Modeling/LinearAlgebra/LinearOperator.h"

#include <vector>

namespace muq{
namespace Modeling{

    /** @class BlockRowOperator
        @ingroup LinearOperators
        @brief Defines a block row matrix in terms of other linear operators
        @details This class defines a \f$1\times N\f$ block matrix of the form
\f[
A = \left[\begin{array}{cccc} A_{1,1} & A_{1,2} & \cdots & A_{1,N} \end{array}\right]
\f]
where each block \f$A_{1,i}\f$ is defined through another linear operator.   Note that the blocks do not need to have the same number of columns, but must have the same number of rows.
    */
    class BlockRowOperator : public LinearOperator
    {

    public:

        BlockRowOperator(std::vector<std::shared_ptr<LinearOperator>> const& blocksIn);

        virtual Eigen::MatrixXd Apply(Eigen::Ref<const Eigen::MatrixXd> const& x) override;
        
        virtual Eigen::MatrixXd ApplyTranspose(Eigen::Ref<const Eigen::MatrixXd> const& x) override;
        
        virtual Eigen::MatrixXd GetMatrix() override;

        std::shared_ptr<LinearOperator> GetBlock(int i) const{return blocks.at(i);};
        std::vector<std::shared_ptr<LinearOperator>> const& GetBlocks() const{return blocks;};
        
    private:
        std::vector<std::shared_ptr<LinearOperator>> blocks;

        static int SumCols(std::vector<std::shared_ptr<LinearOperator>> const& blocksIn);
        
    };
}
}



        



#endif // BLOCKROWOPERATOR_H
