#ifndef DISTRIBUTEDCOLLECTION_H_
#define DISTRIBUTEDCOLLECTION_H_

#include "MUQ/config.h"

#if MUQ_HAS_MPI

#if !MUQ_HAS_PARCER
#error
#endif

#include <parcer/Communicator.h>

#include "MUQ/SamplingAlgorithms/SampleCollection.h"

namespace muq {
  namespace SamplingAlgorithms {
    class DistributedCollection : public SampleCollection {
    public:

      DistributedCollection(std::shared_ptr<SampleCollection> collection, std::shared_ptr<parcer::Communicator> comm);

      virtual ~DistributedCollection() = default;

      /// Add a sample that is stored on this processor
      /**
	 @param[in] newSamp The sample to be added
       */
      virtual void Add(std::shared_ptr<SamplingState> newSamp) override;

      /// Get the local state at the \f$i^{th}\f$ index
      /**
	 @param[in] i The local index
       */
      std::shared_ptr<SamplingState> LocalAt(unsigned i);

      /// Get the local state at the \f$i^{th}\f$ index
      /**
	 @param[in] i The local index
       */
      const std::shared_ptr<SamplingState> LocalAt(unsigned i) const;

      /// Get the global state at the \f$i^{th}\f$ index
      /**
	 @param[in] i The global index
       */
      std::shared_ptr<SamplingState> GlobalAt(unsigned i);

      /// Get the global state at the \f$i^{th}\f$ index
      /**
	 @param[in] i The global index
       */
      const std::shared_ptr<SamplingState> GlobalAt(unsigned i) const;

      /// Get the local state at the \f$i^{th}\f$ index
      /**
	 @param[in] i The local index
       */
      virtual std::shared_ptr<SamplingState> at(unsigned i) override;

      /// Get the local state at the \f$i^{th}\f$ index
      /**
	 @param[in] i The local index
       */
      virtual const std::shared_ptr<SamplingState> at(unsigned i) const override;

      /// The number of samples stored locally
      /**
	 \return The number of samples stored on this processor
       */
      unsigned int LocalSize() const;

      /// The total number of samples
      /**
	 \return The sum of the number of samples stored on each processor
       */
      unsigned int GlobalSize() const;

      /// The total number of samples
      /**
	 By default the size returns the global size.
	 \return The sum of the number of samples stored on each processor
       */
      virtual unsigned int size() const override;

      /// Computes the componentwise central moments (e.g., variance, skewness, kurtosis, etc..) of a specific order
      /**
	 @param[in] order The order of the central moment
	 @param[in] blockDim Compute the central moment of this block
	 \return The central moment using only local samples
       */
      Eigen::VectorXd LocalCentralMoment(unsigned order, int blockDim=-1) const;

      /// Computes the componentwise central moments (e.g., variance, skewness, kurtosis, etc..) of a specific order
      /**
	 @param[in] order The order of the central moment
	 @param[in] blockDim Compute the central moment of this block
	 \return The central moment using all of the samples
       */
      Eigen::VectorXd GlobalCentralMoment(unsigned order, int blockDim=-1) const;

      /// Computes the componentwise central moments (e.g., variance, skewness, kurtosis, etc..) of a specific order
      /**
	 Defaults to using all of the samples.
	 @param[in] order The order of the central moment
	 @param[in] blockDim Compute the central moment of this block
	 \return The central moment using all of the samples
       */
      virtual Eigen::VectorXd CentralMoment(unsigned order, int blockDim=-1) const override;

      /// Compute the mean using only local samples
      /**
	 @param[in] blockDim Compute the mean of this block
	 \return The local mean
      */
      Eigen::VectorXd LocalMean(int blockDim=-1) const;

      /// Compute the mean using all of the samples
      /**
	 @param[in] blockDim Compute the mean of this block
	 \return The global mean
      */
      Eigen::VectorXd GlobalMean(int blockDim=-1) const;

      /// Compute the mean using all of the samples
      /**
	 Defaults to the global mean
	 @param[in] blockDim Compute the mean of this block
	 \return The global mean
      */
      virtual Eigen::VectorXd Mean(int blockDim=-1) const override;

      /// Compute the variance using only local samples
      /**
	 @param[in] blockDim Compute the variance of this block
	 \return The local variance
      */
      Eigen::VectorXd LocalVariance(int blockDim=-1) const;

      /// Compute the variance using all of the samples
      /**
	 @param[in] blockDim Compute the variance of this block
	 \return The global variance
      */
      Eigen::VectorXd GlobalVariance(int blockDim=-1) const;

      /// Compute the variance using all of the samples
      /**
	 Defaults to using the global variance.
	 @param[in] blockDim Compute the variance of this block
	 \return The global variance
      */
      virtual Eigen::VectorXd Variance(int blockDim=-1) const override;

      /// Compute the covariance using only local samples
      /**
	 @param[in] blockDim Compute the covariance of this block
	 \return The local covariance
      */
      Eigen::MatrixXd LocalCovariance(int blockDim=-1) const;

      /// Compute the covariance using all of the samples
      /**
	 @param[in] blockDim Compute the covariance of this block
	 \return The global covariance
      */
      Eigen::MatrixXd GlobalCovariance(int blockDim=-1) const;

      /// Compute the covariance using all of the samples
      /**
	 Defaults to using the global covariance.
	 @param[in] blockDim Compute the covariance of this block
	 \return The global covariance
      */
      virtual Eigen::MatrixXd Covariance(int blockDim=-1) const override;

      /// Compute the ESS using only local samples
      /**
	 @param[in] blockDim Compute the ESS of this block
	 \return The local ESS
      */
      Eigen::VectorXd LocalESS(int blockDim=-1) const;

      /// Compute the ESS using all of the samples
      /**
	 @param[in] blockDim Compute the ESS of this block
	 \return The global ESS
      */
      Eigen::VectorXd GlobalESS(int blockDim=-1) const;

      /// Compute the ESS using all of the samples
      /**
	 Defaults to using the global ESS.
	 @param[in] blockDim Compute the ESS of this block
	 \return The global ESS
      */
      virtual Eigen::VectorXd ESS(int blockDim=-1) const override;

      /// Return all of the samples on the local processor as a matrix
      /**
	 @param[in] blockDim Return the samples of this block
	 \return Each column is a sample
       */
      Eigen::MatrixXd AsLocalMatrix(int blockDim=-1) const;

      /// Return all of the samples as a matrix
      /**
	 @param[in] blockDim Return the samples of this block
	 \return Each column is a sample
       */
      Eigen::MatrixXd AsGlobalMatrix(int blockDim=-1) const;

      /// Return all of the samples as a matrix
      /**
	 @param[in] blockDim Return the samples of this block
	 \return Each column is a sample
       */
      virtual Eigen::MatrixXd AsMatrix(int blockDim=-1) const override;

      /// Return all of the weights on this processor as a vector
      /**
	 @param[in] blockDim Return the weights of this block
	 \return The weights
       */
      Eigen::VectorXd LocalWeights() const;

      /// Return all of the weights as a vector
      /**
	 @param[in] blockDim Return the weights of this block
	 \return The weights
       */
      Eigen::VectorXd GlobalWeights() const;

      /// Return all of the weights as a vector
      /**
	 @param[in] blockDim Return the weights of this block
	 \return The weights
       */
      virtual Eigen::VectorXd Weights() const override;

      /**
      Writes the local samples to file
	 @param[in] filename The name of the file
	 @param[in] dataset The name of the group within the file
      */
      virtual void WriteToFile(std::string const& filename, std::string const& dataset = "/") const override;

      Eigen::VectorXd GlobalExpectedValue(std::shared_ptr<muq::Modeling::ModPiece> const& f, std::vector<std::string> const& metains = std::vector<std::string>()) const;

      Eigen::VectorXd LocalExpectedValue(std::shared_ptr<muq::Modeling::ModPiece> const& f, std::vector<std::string> const& metains = std::vector<std::string>()) const;

      virtual Eigen::VectorXd ExpectedValue(std::shared_ptr<muq::Modeling::ModPiece> const& f, std::vector<std::string> const& metains = std::vector<std::string>()) const override;

    private:

      template <class scalar, int rows, int cols, int options, int maxRows, int maxCols>
    	Eigen::Matrix<scalar, rows, cols, options, maxRows, maxCols> GlobalEigenMean(Eigen::Matrix<scalar, rows, cols, options, maxRows, maxCols> const& local) const {

        	typedef Eigen::Matrix<scalar, rows, cols, options, maxRows, maxCols> MatType;
        	MatType global = MatType::Zero(local.rows(), local.cols());

        	for( unsigned int i=0; i<comm->GetSize(); ++i ) {
        	  MatType l(local.rows(), local.cols());
        	  if( comm->GetRank()==i ) { l = local; }
        	  comm->Bcast(l, i);

        	  global += l/(double)comm->GetSize();
        	}

        	return global;
      }

      /// The local sample collection (stored on this processor)
      std::shared_ptr<SampleCollection> collection;

      /// The communicator for this collection
      std::shared_ptr<parcer::Communicator> comm;
    };
  } // namespace SamplingAlgorithms
} // namespace muq

#endif // end MUQ_HAS_MPI
#endif
