#include "MUQ/Inference/Filtering/KalmanSmoother.h"

#include <Eigen/Dense>

using namespace muq::Modeling;
using namespace muq::Inference;



Eigen::MatrixXd KalmanSmoother::ComputeC(Eigen::MatrixXd                          const& currDist_t_cov,
                                         Eigen::MatrixXd                          const& nextDist_t_cov,
                                         std::shared_ptr<muq::Modeling::LinearOperator> F)
{
    return  nextDist_t_cov.llt().solve( F->Apply(currDist_t_cov) ).transpose();
}

Eigen::MatrixXd KalmanSmoother::ComputeC( Eigen::MatrixXd const& currDist_t_cov,
                                          Eigen::MatrixXd const& nextDist_t_cov,
                                          Eigen::MatrixXd const& F)
{
    return nextDist_t_cov.llt().solve( F*currDist_t_cov ).transpose();
}
