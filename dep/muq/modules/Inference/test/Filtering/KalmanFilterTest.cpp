#include "gtest/gtest.h"

#include "MUQ/Inference/Filtering/KalmanFilter.h"
#include "MUQ/Modeling/LinearAlgebra/LinearOperator.h"
#include "MUQ/Modeling/LinearAlgebra/EigenLinearOperator.h"
#include "MUQ/Utilities/Exceptions.h"

#include <Eigen/Dense>

using namespace muq::Inference;
using namespace muq::Modeling;


class Inference_KalmanFilter : public ::testing::Test
{

protected:

    virtual void SetUp()
    {
        prior.first = Eigen::VectorXd::Ones(stateDim);
        prior.second = Eigen::MatrixXd::Identity(stateDim,stateDim);

        Eigen::MatrixXd Hmat = Eigen::MatrixXd::Zero(obsDim, stateDim);
        Hmat(0,0) = 1.0;
        Hmat(0,2) = 1.0;
        Hmat(1,1) = 1.0;
        Hmat(1,3) = 1.0;

        H = LinearOperator::Create(Hmat);

        data = Hmat*Eigen::VectorXd::Random(stateDim);

        // Now compute the true posterior
        Eigen::VectorXd y = data - Hmat*prior.first;
        Eigen::MatrixXd S = obsVar*Eigen::MatrixXd::Identity(obsDim, obsDim) + Hmat*prior.second*Hmat.transpose();

        Eigen::MatrixXd K = S.llt().solve(Hmat*prior.second).transpose();
        truePost.first = prior.first + K*y;
        truePost.second = prior.second - K*Hmat*prior.second;
        truePost.second = 0.5*(truePost.second + truePost.second.transpose()).eval();

    };

    virtual void TearDown(){
        if(post.first.size()>0){
            for(int j=0; j<stateDim; ++j){
                EXPECT_DOUBLE_EQ(truePost.first(j), post.first(j));
                for(int i=0; i<stateDim; ++i)
                    EXPECT_DOUBLE_EQ(truePost.second(i,j), post.second(i,j));
            }
        }
    }

    const int stateDim = 4;
    const int obsDim = 2;
    const double obsVar = 1e-4;

    std::pair<Eigen::VectorXd, Eigen::MatrixXd> prior;
    std::pair<Eigen::VectorXd, Eigen::MatrixXd> post;
    std::shared_ptr<LinearOperator> H;
    Eigen::VectorXd data;


    std::pair<Eigen::VectorXd, Eigen::MatrixXd> truePost;
};


TEST_F(Inference_KalmanFilter, SizeFail)
{
    Eigen::VectorXd obsCov = obsVar*Eigen::VectorXd::Ones(obsDim);
    Eigen::MatrixXd Hmat = Eigen::MatrixXd::Ones(stateDim, obsDim-1);
    auto H2 = LinearOperator::Create(Hmat);

    EXPECT_THROW(KalmanFilter::Analyze(prior, H2, data, obsCov), muq::WrongSizeError);

    std::pair<Eigen::VectorXd, Eigen::MatrixXd> prior2(prior);
    prior2.first = Eigen::VectorXd::Ones(stateDim-1);
    EXPECT_THROW(KalmanFilter::Analyze(prior2, H, data, obsCov), muq::WrongSizeError);
}


TEST_F(Inference_KalmanFilter, VectorNoise)
{
    Eigen::VectorXd obsCov = obsVar*Eigen::VectorXd::Ones(obsDim);
    post = KalmanFilter::Analyze(prior, H, data, obsCov);
}

TEST_F(Inference_KalmanFilter, FullNoise)
{
    Eigen::MatrixXd obsCov = 1e-4*Eigen::MatrixXd::Identity(obsDim, obsDim);
    post = KalmanFilter::Analyze(prior, H, data, obsCov);
}
