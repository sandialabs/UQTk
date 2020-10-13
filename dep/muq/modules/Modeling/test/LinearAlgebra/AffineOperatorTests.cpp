#include "MUQ/Modeling/LinearAlgebra/AffineOperator.h"
#include "MUQ/Modeling/LinearAlgebra/EigenLinearOperator.h"

#include <random>

#include <gtest/gtest.h>

#include <Eigen/Core>
#include <Eigen/SparseCore>

#include "MUQ/Utilities/Exceptions.h"

using namespace muq::Modeling;

TEST(Utilties_AffineOperator, EigenDense)
{
    const int rows = 5;
    const int cols = 2;

    Eigen::MatrixXd A = Eigen::MatrixXd::Random(rows,cols);
    Eigen::VectorXd b = Eigen::VectorXd::Random(rows);

    auto denseOp = AffineOperator::Create(A,b);

    Eigen::VectorXd x = Eigen::VectorXd::Ones(cols);
    Eigen::VectorXd ytrue, ytest;

    ytrue = b+A*x;
    ytest = denseOp->Evaluate(x).at(0);

    for(int i=0; i<rows; ++i)
        EXPECT_DOUBLE_EQ(ytrue(i), ytest(i));
}
