#include "MUQ/Inference/Filtering/KalmanFilter.h"

#include "MUQ/Utilities/Exceptions.h"

#include <Eigen/Dense>

using namespace muq::Modeling;
using namespace muq::Inference;


Eigen::MatrixXd KalmanFilter::ComputeGain(Eigen::MatrixXd                           const& HP,
                                          std::shared_ptr<muq::Modeling::LinearOperator>  H,
                                          Eigen::Ref<const Eigen::MatrixXd> const&         obsCov)
{

    Eigen::MatrixXd S = H->Apply( HP.transpose() );

    if((obsCov.rows()==1)||(obsCov.cols()==1)){
        S += obsCov.asDiagonal();
    }else{
        S += obsCov;
    }

    Eigen::LLT<Eigen::MatrixXd> solver(S); // <- In place Cholesky decomposition
    Eigen::MatrixXd K = solver.solve( HP ).transpose();

    return K;
}



std::pair<Eigen::VectorXd, Eigen::MatrixXd> KalmanFilter::Analyze(std::pair<Eigen::VectorXd, Eigen::MatrixXd> const& dist,
                                                                  std::shared_ptr<muq::Modeling::LinearOperator>    H,
                                                                  Eigen::Ref<const Eigen::VectorXd> const&           obsMean,
                                                                  Eigen::Ref<const Eigen::MatrixXd> const&           obsCov)
{


    const int obsDim = std::max(obsCov.rows(), obsCov.cols());
    if(H->rows() != obsDim)
        throw muq::WrongSizeError("In KalmanFilter::Analyze: The size of the observation noise covariance does not match the size of the observation operator. The observation operator returns a vector with  " + std::to_string(H->rows()) + " components, but the noise covariance is a " + std::to_string(obsDim) + "x" + std::to_string(obsDim) + " matrix.");

    if(H->cols() != dist.first.rows())
        throw muq::WrongSizeError("In KalmanFilter::Analyze: The size of the observation operator (" + std::to_string(H->rows()) + "x" + std::to_string(H->cols()) + ") does not match the size of the prior mean (" + std::to_string(dist.first.rows()) + ").");

    Eigen::VectorXd y = obsMean - H->Apply(dist.first);

    Eigen::MatrixXd HP = H->Apply(dist.second);

    // Compute the Kalman Gain
    Eigen::MatrixXd K = ComputeGain(HP, H, obsCov);

    // Posterior Mean and Covariance
    std::pair<Eigen::VectorXd, Eigen::MatrixXd> output;

    output.first = dist.first + K*y;

    output.second = dist.second - K*HP;

    return output;
}
