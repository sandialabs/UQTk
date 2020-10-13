#ifndef SAMPLECOLLECTION_H
#define SAMPLECOLLECTION_H

#include <memory>
#include <vector>

#include <Eigen/Core>

#include "MUQ/Modeling/ModPiece.h"

#include "MUQ/SamplingAlgorithms/SamplingState.h"

#include "MUQ/Utilities/HDF5/HDF5File.h"

namespace muq{
  namespace SamplingAlgorithms{

    class SamplingStateIdentity{
    public:
      SamplingStateIdentity(int blockIndIn) : blockInd(blockIndIn){};

      virtual ~SamplingStateIdentity() = default;

      Eigen::VectorXd const& operator()(SamplingState const& a);

      const int blockInd;

    private:
      Eigen::VectorXd output;
    };

    class ExpectedModPieceValue {
    public:
      ExpectedModPieceValue(std::shared_ptr<muq::Modeling::ModPiece> const& f, std::vector<std::string> const& metains);

      virtual ~ExpectedModPieceValue() = default;

      Eigen::VectorXd const& operator()(SamplingState const& a);

    private:
      std::shared_ptr<muq::Modeling::ModPiece> f;

      const std::vector<std::string> metains;
    };

    class SamplingStatePartialMoment{
    public:
      SamplingStatePartialMoment(int                    blockIndIn,
                                 int                    momentPowerIn,
                                 Eigen::VectorXd const& muIn) : blockInd(blockIndIn), momentPower(momentPowerIn), mu(muIn){};

      virtual ~SamplingStatePartialMoment() = default;

      Eigen::VectorXd const& operator()(SamplingState const& a);

      const int blockInd;
      const int momentPower;
      const Eigen::VectorXd& mu;

    private:
      Eigen::VectorXd output;
    };

    /** @ingroup SamplingAlgorithms
        @class SampleCollection
        @brief A class to hold and analyze a collection of SamplingState objects
    */
    class SampleCollection{
    public:
      SampleCollection() = default;

      virtual ~SampleCollection() = default;

      virtual void Add(std::shared_ptr<SamplingState> newSamp);

      virtual std::shared_ptr<SamplingState> at(unsigned i);
      virtual const std::shared_ptr<SamplingState> at(unsigned i) const;

      virtual unsigned size() const;

      ///  Computes the componentwise central moments (e.g., variance, skewness, kurtosis, etc..) of a specific order
      virtual Eigen::VectorXd CentralMoment(unsigned order, int blockNum=-1) const;

      ///  Computes the componentwise central moments (e.g., variance, skewness, kurtosis, etc..) of a specific order given that we already know the mean
      virtual Eigen::VectorXd CentralMoment(unsigned order, Eigen::VectorXd const& mean, int blockNum=-1) const;

      virtual Eigen::VectorXd Mean(int blockInd=-1) const;
      virtual Eigen::VectorXd Variance(int blockInd=-1) const;
      virtual Eigen::MatrixXd Covariance(int blockInd=-1) const;
      virtual Eigen::MatrixXd Covariance(Eigen::VectorXd const& mean, int blockInd=-1) const;

      virtual std::vector<Eigen::MatrixXd> RunningCovariance(int blockInd=-1) const;
      virtual std::vector<Eigen::MatrixXd> RunningCovariance(Eigen::VectorXd const& mean, int blockInd=-1) const;

      /** @brief Returns the effective sample size of the samples
          @details For almost all random variables of interest, the central limit
                   theorem states that the variance of a Monte Carlo estimator
                   \f[\hat{\mu}_N = \frac{1}{N}\sum_{i=1}^N \theta^{(i)},\f]
                   satisfies
                   \f[ \sqrt{N}\mathbb{V}\left[\hat{\mu}_N-\mu\right] \rightarrow N(0,\sigma^2),\f]
                   as the number of samples \f$N\f$ increases.  Here, \f$mu\f$ is
                   the true mean of the random variable, \f$\sigma^2\f$ is the true
                   variance of the random variable.  This assumes that each \f$\theta^{(i)}\f$
                   is an independent sample with unity weight.

                   When the samples are not independent or the weights are not unity,
                   the variance of the estimator will generally be larger.   The
                   "effective sample size" (ESS) describes how many independent
                   samples would have been needed to obtain the same estimator
                   variance.  In particular, let \f$\tilde{\mu}_N\f$ be a Monte
                   Carlo estimator based on \f$N\f$ correlated samples.  The ESS
                   is then given by the ratio of the esimator variances:
                   \f[ ESS = N \frac{\mathbb{V}[\hat{\mu}-\mu]}{\mathbb{V}[\tilde{\mu}-\mu]}. \f]

                   In SampleCollection, the samples are assumed independent, but
                   not equally weighted, which is typically the case with importance
                   sampling.  In this setting, the ESS is estimated using
                   \f[
                      ESS = \frac{\left(\sum_{i=1}^N w_i}\right)^2}{\sum_{i=1}^N w_i^2}.
                   \f]

                   Note that children of this class may compute the ESS with different
                   approaches.  For example, the MarkovChain class computes the
                   ESS using the approach of "Monte Carlo errors with less error" by
                   Ulli Wolff.
      */
      virtual Eigen::VectorXd ESS(int blockDim=-1) const;

      virtual Eigen::MatrixXd AsMatrix(int blockDim=-1) const;

      virtual Eigen::VectorXd Weights() const;

      /**
	 @param[in] filename The name of the file
	 @param[in] dataset The name of the group within the file
       */
      virtual void WriteToFile(std::string const& filename, std::string const& dataset = "/") const;

      /**
      @param[in] name Need the this piece of meta data for each sample
      \return A matrix of meta data associated with in the input string
      */
      Eigen::MatrixXd GetMeta(std::string const& name) const;

      virtual Eigen::VectorXd ExpectedValue(std::shared_ptr<muq::Modeling::ModPiece> const& f, std::vector<std::string> const& metains = std::vector<std::string>()) const;

      std::vector<Eigen::VectorXd> RunningExpectedValue(std::shared_ptr<muq::Modeling::ModPiece> const& f, std::vector<std::string> const& metains = std::vector<std::string>()) const;

    protected:

      std::vector<std::shared_ptr<SamplingState>> samples;

      /** Returns the sum of the weights and the sum of the squared weights. */
      static std::pair<double,double> RecursiveWeightSum(std::vector<std::shared_ptr<SamplingState>>::const_iterator start,
                                                         std::vector<std::shared_ptr<SamplingState>>::const_iterator end);


      template<typename FuncType>
      static std::pair<double,Eigen::VectorXd> RecursiveSum(std::vector<std::shared_ptr<SamplingState>>::const_iterator                         start,
                                                            std::vector<std::shared_ptr<SamplingState>>::const_iterator                         end,
                                                            FuncType& f)
      {
        int numSamps = std::distance(start,end);
        const int maxSamps = 20;

        // If the number of samples is small enough, we can safely add them up directly
        if(numSamps<maxSamps){

          Eigen::VectorXd sum = (*start)->weight * f(**start);
          double weightSum = (*start)->weight;

          for(auto it=start+1; it!=end; ++it){
              sum += (*it)->weight * f(**it);
              weightSum += (*it)->weight;
          }
          return std::make_pair(weightSum, sum);

        // Otherwise, it's more numerically stable to add things pairwise
        }else{
          int halfDist = std::floor(0.5*numSamps);
          double weight1, weight2;
          Eigen::VectorXd sum1, sum2;
          std::tie(weight1,sum1) = RecursiveSum(start, start+halfDist, f);
          std::tie(weight2,sum2) = RecursiveSum(start+halfDist, end, f);

          return std::make_pair(weight1+weight2, (sum1+sum2).eval());
        }
      }

      /**
	 \return A map from meta data name to a matrix where each column corresponds to a sample
       */
      std::unordered_map<std::string, Eigen::MatrixXd> GetMeta() const;

    private:

      /**
	 @param[in] hdf5file The hdf5 file where the data will be written
	 @param[in] dataname The name of the data set we (may) want to create
	 @param[in] dataSize The number of rows (the size of the data in one sample)
	 @param[in] totSamps The total number of samples we need to write to the file (max. number of samples)
	 \return true: the data set exists and is the right size, false: the data set does not exist or is the wrong size
       */
      bool CreateDataset(std::shared_ptr<muq::Utilities::HDF5File> hdf5file, std::string const& dataname, int const dataSize, int const totSamps) const;
    };
  }
}

#endif // SAMPLECOLLECTION_H
