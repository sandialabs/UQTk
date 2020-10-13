#ifndef CONCATENATEOPERATOR_H
#define CONCATENATEOPERATOR_H

#include "MUQ/Modeling/LinearAlgebra/LinearOperator.h"

#include <vector>

namespace muq
{
namespace Modeling
{

    /** @class ConcatenateOperator
        @ingroup LinearOperators
        @brief Vertical or horizontal concatenation of other linear operators.
        @details Given two linear operators \f$A\f$ and \f$B\f$, this operator defines either the vertically stacked operator
\f[
C = \left[\begin{array}{c} A\\ B\end{array}\right],
\f]
or the horizontally stacked operator
\f[
C = \left[A, B\right].
\f]
    */
    class ConcatenateOperator : public LinearOperator
    {
    public:
        ConcatenateOperator(std::vector<std::shared_ptr<LinearOperator>> const& opsIn,
                            const int                                           rowColIn);

        virtual ~ConcatenateOperator(){};
        
        /** Apply the linear operator to a vector */
        virtual Eigen::MatrixXd Apply(Eigen::Ref<const Eigen::MatrixXd> const& x) override;
  
        /** Apply the transpose of the linear operator to a vector. */
        virtual Eigen::MatrixXd ApplyTranspose(Eigen::Ref<const Eigen::MatrixXd> const& x) override;


        virtual Eigen::MatrixXd GetMatrix() override;
        
        static std::shared_ptr<ConcatenateOperator> VStack(std::shared_ptr<LinearOperator> Ain,
                                                           std::shared_ptr<LinearOperator> Bin);

        static std::shared_ptr<ConcatenateOperator> HStack(std::shared_ptr<LinearOperator> Ain,
                                                           std::shared_ptr<LinearOperator> Bin);


        
    private:
        std::vector<std::shared_ptr<LinearOperator>> ops;
        
        const int rowCol; // zero if stacked vertically along rows, 1 if stacked horizontally

        static int GetRows(std::vector<std::shared_ptr<LinearOperator>> const& opsIn,
                           const int                                           rowColIn);
        
        static int GetCols(std::vector<std::shared_ptr<LinearOperator>> const& opsIn,
                           const int                                           rowColIn);

        void CheckSizes();
        
    }; // class ConcatenateOperator

    

} // namespace Modeling
} // namespace muq


#endif // #ifndef CONCATENATEOPERATOR_H
