#ifndef KALMANFILTER_H
#define KALMANFILTER_H

#include "MUQ/Modeling/LinearAlgebra/LinearOperator.h"

#include <Eigen/Core>

#include <memory>


namespace muq
{
namespace Inference
{

    class KalmanFilter
    {

    public:
        
        static std::pair<Eigen::VectorXd, Eigen::MatrixXd> Analyze(std::pair<Eigen::VectorXd, Eigen::MatrixXd> const& dist,
                                                                   std::shared_ptr<muq::Modeling::LinearOperator>    H,
                                                                   Eigen::Ref<const Eigen::VectorXd> const&           obsMean,
                                                                   Eigen::Ref<const Eigen::MatrixXd> const&           obsCov);

    private:
        
        static Eigen::MatrixXd ComputeGain(Eigen::MatrixXd                           const& HP,
                                           std::shared_ptr<muq::Modeling::LinearOperator>  H,
                                           Eigen::Ref<const Eigen::MatrixXd> const&         obsCov);
        
        
    }; // class KalmanFilter

} // namespace Inference
} // namespace muq

#endif
