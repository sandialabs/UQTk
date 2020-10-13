#ifndef LINEARSDE_H
#define LINEARSDE_H


#include "MUQ/Modeling/LinearAlgebra/LinearOperator.h"
#include "MUQ/Modeling/LinearAlgebra/EigenLinearOperator.h"

#include <boost/property_tree/ptree.hpp>
#include <random>

namespace muq
{
namespace Modeling
{

    /** @brief Defines a linear time invariant stochastic differential equation.
        @details This class defines a LTI SDE of the form
\f[
\frac{\partial f(t)}{\partial t} = F f(t) + L w(t),
\f]
where \f$f(t)\f$ is the solution in \f$\mathbb{R}^M\f$, \f$F\f$ is an \f$M\timesM\f$ matrix, \f$L\f$ is an \f$M\times N\f$ matrix, and \f$w(t)\f$ is a white noise process with an \f$N\timesN\f$ covariance matrix \f$Q\f$.
    */
    class LinearSDE
    {

    public:

        template<typename Derived1, typename Derived2>
        LinearSDE(Eigen::Matrix<Derived1,Eigen::Dynamic, Eigen::Dynamic>   const& Fin,
                  Eigen::Matrix<Derived2,Eigen::Dynamic, Eigen::Dynamic>   const& Lin,
                  Eigen::MatrixXd                                          const& Qin,
                  boost::property_tree::ptree options) : LinearSDE(muq::Modeling::LinearOperator::Create(Fin),
                                                                   muq::Modeling::LinearOperator::Create(Lin),
                                                                   Qin,
                                                                   options)
        {};
        
        
        LinearSDE(std::shared_ptr<muq::Modeling::LinearOperator>    Fin,
                  std::shared_ptr<muq::Modeling::LinearOperator>    Lin,
                  Eigen::MatrixXd                             const& Qin,
                  boost::property_tree::ptree                        options);


        /** Given \f$f(t)\f$, the state of the system at time \f$t\f$, return a random realization of the state at time \f$t+\delta t\f$.
         */
        Eigen::VectorXd EvolveState(Eigen::VectorXd const& f0,
                                    double                 T) const;

        /** Given the mean and covariance of the solution at time \f$t\f$, compute the mean and covariance of the solution at time \f$t+T\f$.
         */
        std::pair<Eigen::VectorXd, Eigen::MatrixXd> EvolveDistribution(Eigen::VectorXd const& muIn,
                                                                       Eigen::MatrixXd const& gammaIn,
                                                                       double                 T) const; 

        /** Evolve the mean and covariance of the system using a std::pair to hold the distribution.
         */
        std::pair<Eigen::VectorXd, Eigen::MatrixXd> EvolveDistribution(std::pair<Eigen::VectorXd,Eigen::MatrixXd> const& muCov,
                                                                       double                                            T) const{
            return EvolveDistribution(muCov.first, muCov.second, T);
        }; 


        /** 
           Compute a matrix A and covariance Q such that \f$x(t+\delta t) = A x(t) + q\f$ where \f$q\f$ is a normal random variable with covariance \f$Q\f$.
         */
        std::pair<Eigen::MatrixXd, Eigen::MatrixXd> Discretize(double deltaT);
        
        /** @brief Combines the states of multiple SDEs into a single monolitch SDE.
            @details Consider \f$N\f$ different stochastic differential equations defined by matrices \f$F_i\f$, \f$L_i\f$, and process covariances \f$Q_i\f$.   This function creates a new SDE defined by block diagonal matrices \f$F\f$, \f$L\f$, and \f$Q\f$:
\f[
F = \left[\begin{array}{cccc} F_1 & 0 & \cdots & \\ 0 & F_2 & 0 \\ \vdots & & \ddots & \\ 0 & \cdots & & F_N \end{array}\right]
\f]
\f[
L = \left[\begin{array}{cccc} L_1 & 0 & \cdots & \\ 0 & L_2 & 0 \\ \vdots & & \ddots & \\ 0 & \cdots & & L_N \end{array}\right]
\f]
\f[
Q = \left[\begin{array}{cccc} Q_1 & 0 & \cdots & \\ 0 & Q_2 & 0 \\ \vdots & & \ddots & \\ 0 & \cdots & & Q_N \end{array}\right]
\f]
         */
        static std::shared_ptr<LinearSDE> Concatenate(std::vector<std::shared_ptr<LinearSDE>> const& sdes,
                                                      boost::property_tree::ptree                    options = boost::property_tree::ptree());
        

        /// The dimension of the state variable \f$f(t)\f$.
        const int stateDim;


        std::shared_ptr<muq::Modeling::LinearOperator> GetF() const{return F;};
        std::shared_ptr<muq::Modeling::LinearOperator> GetL() const{return L;};
        Eigen::MatrixXd const& GetQ() const{return Q;};
        
        
    protected:

        void ExtractOptions(boost::property_tree::ptree options);
        
        std::shared_ptr<muq::Modeling::LinearOperator> F;
        std::shared_ptr<muq::Modeling::LinearOperator> L;

        Eigen::MatrixXd Q;
        Eigen::MatrixXd sqrtQ;

        double dt; // time step used in SDE integration

    };


}// namespace Modeling
}// namespace muq




#endif
