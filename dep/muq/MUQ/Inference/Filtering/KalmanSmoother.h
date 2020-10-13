#ifndef KALMANSMOOTHER_H
#define KALMANSMOOTHER_H

#include "MUQ/Inference/Filtering/KalmanFilter.h"

namespace muq
{
namespace Inference
{

    /** @class KalmanSmoother
        Implements the Rauch–Tung–Striebel smoother.
    */
    class KalmanSmoother
    {

    public:

        /** @param[in] currDist_t The distribution at time t after the forward Kalman filtering step (i.e., using all data up to and including time t).
            @param[in] nextDist_t The distribution at time t+1 from the prediction phase of the Kalman filtering step (i.e., using all data up to time t).
            @param[in] nextDist_n The distribution at time t+1 after the RTS smoothing step (i.e., using all data).
            @param[in] F The linear operator acting on the state at time t, to produce the state at time t+1
            @returns A distribution (mean and covariance) at time t after accounting for all data.
        */
        template<typename FType>
        static std::pair<Eigen::VectorXd, Eigen::MatrixXd> Analyze(std::pair<Eigen::VectorXd, Eigen::MatrixXd> const& currDist_t,
                                                                   std::pair<Eigen::VectorXd, Eigen::MatrixXd> const& nextDist_t,
                                                                   std::pair<Eigen::VectorXd, Eigen::MatrixXd> const& nextDist_n,
                                                                   FType                                       const& F)
        {
            std::pair<Eigen::VectorXd, Eigen::MatrixXd> output;
            
            Eigen::MatrixXd C = ComputeC(currDist_t.second, nextDist_t.second, F);
            
            output.first = currDist_t.first + C*(nextDist_n.first - nextDist_t.first);
            output.second = currDist_t.second + C*(nextDist_n.second - nextDist_t.second).selfadjointView<Eigen::Lower>()*C.transpose();
            
            return output;

        }

        
        
    private:                                                           

        static Eigen::MatrixXd ComputeC( Eigen::MatrixXd                          const& currDist_t_cov,
                                         Eigen::MatrixXd                          const& nextDist_t_cov,
                                         std::shared_ptr<muq::Modeling::LinearOperator> F);

        
        static Eigen::MatrixXd ComputeC( Eigen::MatrixXd const& currDist_t_cov,
                                         Eigen::MatrixXd const& nextDist_t_cov,
                                         Eigen::MatrixXd const& F);
        

    }; // class KalmanSmoother 
    
}// namespace Inference
}// namespace muq


#endif 
