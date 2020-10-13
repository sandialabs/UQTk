#ifndef BLOCKDIAGONALOPERATOR_H
#define BLOCKDIAGONALOPERATOR_H

#include "MUQ/Modeling/LinearAlgebra/LinearOperator.h"

#include <vector>

namespace muq{
namespace Modeling{

    /** @class BlockDiagonalOperator
        @ingroup LinearOperators
        @brief Defines a block diagonal linear operator in terms of other linear operators
        @details This class defines a \f$N\times N\f$ matrix of the form
\f[
A = \left[\begin{array}{cccc} A_{1,1} & 0 & & \\ 0 & A_{2,2} & 0 & \\ \ddots & \ddots & \ddots & \\ & & 0 & A_{N,N} \end{array}\right],
\f]
where each block \f$A_{i,i}\f$ is defined through another linear operator.   Note that the blocks do not need to be the same size or square.
    */
    class BlockDiagonalOperator : public LinearOperator
    {

    public:

        BlockDiagonalOperator(std::vector<std::shared_ptr<LinearOperator>> const& blocksIn);

        virtual Eigen::MatrixXd Apply(Eigen::Ref<const Eigen::MatrixXd> const& x) override;
        
        virtual Eigen::MatrixXd ApplyTranspose(Eigen::Ref<const Eigen::MatrixXd> const& x) override;
        
        virtual Eigen::MatrixXd GetMatrix() override;

        std::shared_ptr<LinearOperator> GetBlock(int i) const{return blocks.at(i);};
        std::vector<std::shared_ptr<LinearOperator>> const& GetBlocks() const{return blocks;};
        
    private:
        std::vector<std::shared_ptr<LinearOperator>> blocks;

        static int SumRows(std::vector<std::shared_ptr<LinearOperator>> const& blocksIn);

        static int SumCols(std::vector<std::shared_ptr<LinearOperator>> const& blocksIn);
        
    };
}
}



        



#endif // BLOCKDIAGONALOPERATOR_H
