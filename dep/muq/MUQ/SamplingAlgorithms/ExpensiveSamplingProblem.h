#ifndef EXPENSIVESAMPLINGPROBLEM_H_
#define EXPENSIVESAMPLINGPROBLEM_H_

#include "MUQ/config.h"

#if MUQ_HAS_PARCER
#include <parcer/Communicator.h>
#endif

#include "MUQ/Approximation/Regression/LocalRegression.h"

#include "MUQ/SamplingAlgorithms/SamplingProblem.h"

namespace muq {
  namespace SamplingAlgorithms {
    class ExpensiveSamplingProblem : public SamplingProblem {
    public:

      /**
	 @param[in] target The target distribution
	 @param[in] pt Options for the sampling problem
       */
      ExpensiveSamplingProblem(std::shared_ptr<muq::Modeling::ModPiece> const& target, boost::property_tree::ptree pt);

      /**
   @param[in] target The target distribution
   @param[in] centroid The centroid of the distribution
   @param[in] pt Options for the sampling problem
       */
      ExpensiveSamplingProblem(std::shared_ptr<muq::Modeling::ModPiece> const& target, Eigen::VectorXd const& centroid, boost::property_tree::ptree pt);

#if MUQ_HAS_PARCER
      /**
    @param[in] target The target distribution
    @param[in] pt Options for the sampling problem
    @param[in] comm Parallel communicator for regression
      */
      ExpensiveSamplingProblem(std::shared_ptr<muq::Modeling::ModPiece> const& target, boost::property_tree::ptree pt, std::shared_ptr<parcer::Communicator> comm);

      /**
	 @param[in] target The target distribution
   @param[in] centroid The centroid of the distribution
	 @param[in] pt Options for the sampling problem
	 @param[in] comm Parallel communicator for regression
       */
      ExpensiveSamplingProblem(std::shared_ptr<muq::Modeling::ModPiece> const& target, Eigen::VectorXd const& centroid, boost::property_tree::ptree pt, std::shared_ptr<parcer::Communicator> comm);
#endif

      virtual ~ExpensiveSamplingProblem() = default;

      virtual double LogDensity(unsigned int const t, std::shared_ptr<SamplingState> const& state, AbstractSamplingProblem::SampleType type) override;

      unsigned int CacheSize() const;


      virtual void AddOptions(boost::property_tree::ptree & pt) const override;

    private:

      /// Set up the sampling problem
      /**
	 @param[in] pt Options for the sampling problem
       */
      void SetUp(boost::property_tree::ptree& pt);

      /**
	 @param[in] step The current MCMC step
	 @param[in] state The state where we are evalauting the log target
	 @param[out] neighbors The nearest neighbors
	 @param[out] results The log-target at the nearest neighbors
   \return The error threshold
       */
      double RefineSurrogate(
        unsigned int const step,
        std::shared_ptr<SamplingState> const& state,
        std::vector<Eigen::VectorXd>& neighbors,
        std::vector<Eigen::VectorXd>& results);

        /**
     @param[in] state The state where we are evalauting the log target
     @param[in] radius The radus of the nearest neighbor ball
     @param[out] neighbors The nearest neighbors
     @param[out] results The log-target at the nearest neighbors
     \return The poisedness constant \f$\Lambda\f$
         */
        double RefineSurrogate(
          std::shared_ptr<SamplingState> const& state,
          double const radius,
          std::vector<Eigen::VectorXd>& neighbors,
          std::vector<Eigen::VectorXd>& results);

      /**
	 Replace the point that is farthest from the new point
	 @param[in] state The point where we are evalauting the log target
	 @param[in] index The index of the point that we are replacing
	 @param[out] neighbors The nearest neighbors
	 @param[out] results The log-target at the nearest neighbors
       */
      void RefineSurrogate(
        Eigen::VectorXd const& point,
        unsigned int const index,
        std::vector<Eigen::VectorXd>& neighbors,
        std::vector<Eigen::VectorXd>& results);

      /// Check to make sure we have enough model evaluations
      /**
      Check to see if the state has the nearest neighbors already stored.  If yes, return them, if not, find and return them.
        @param[in] state Find the nearest neighbors closest to this point
        @param[out] neighbors The nearest neighbors
        @param[out] results The log-target at the nearest neighbors
      */
      void CheckNeighbors(
        std::shared_ptr<SamplingState> const& state,
        std::vector<Eigen::VectorXd>& neighbors,
        std::vector<Eigen::VectorXd>& results) const;

      /// Check to make sure we have enough model evaluations
      /**
        @param[in] state If we do not have enough points, sample from a standard Gaussian centered at this point
      */
      void CheckNumNeighbors(std::shared_ptr<SamplingState> const& state);

      /// Compute the error threshold
      /**
        @param[in] x The state of the MCMC chain
        \return The error threshild
      */
      double ErrorThreshold(Eigen::VectorXd const& x) const;

      /// Estimate the Lyapunov function
      /**
        Defaults to \f$V(x) = \exp{(\|x-\bar{x}\|_2)}\f$.
        @param[in] x the input point \f$x\f$
        \return \f$V(x)\f$
      */
      double LogLyapunovFunction(Eigen::VectorXd const& x) const;

      std::shared_ptr<muq::Approximation::LocalRegression> reg;

      /// Parameters for random refinement
      /**
	 Refine with probability \f$\beta = \beta_0 t^{-\beta_1}\f$.  \f$\beta_0\f$ is "BetaScale" and it defaults to \f$0\f$.  \f$\beta_1\f$ is "BetaExponent" and it defaults to \f$RAND_MAX\f$.
       */
      std::pair<double, double> beta;

      /// The length of the first level
      /**
        \f$\tau_0\f$ is "FirstLevelLength" and defaults to \f$1.0\f$.
      */
      double tau0;

      /// Parameters for structural refinement
      /**
	 Refine if the error threshold exceeds \f$\gamma = \gamma_0 l^{-\gamma_1}\f$, where \f$l\f$ is the current error threshold level (ExpensiveSamplingProblem::level).  \f$\gamma_0\f$ is "GammaScale" and it defaults to \f$1\f$.  \f$\gamma_1\f$ is "GammaExponent" and it defaults to \f$1.0\f$.
       */
      std::pair<double, double> gamma;

      /// The scaling and exponent for the Lyapunov function
      std::pair<double, double> lyapunovPara;

      /// The centroid of the distribution (input parameter)
      const Eigen::VectorXd centroid;

      /// Tail correction parameter
      double eta;

      /// The current error threshold level
      unsigned int level = 1;

      /// Cumulative beta refinements
      unsigned int cumbeta = 0;

      /// Cumulative gamma refinements
      unsigned int cumgamma = 0;

      /// Cumulative kappa refinements
      unsigned int cumkappa = 0;

      /// The Lyapunov function at the most recently accepted point
      double lastLyapunov = 0.0;

      /// The error threshold at the most recently accepted point
      double lastThreshold = 0.0;
    };
  } // namespace SamplingAlgorithms
} // namespace muq

#endif
