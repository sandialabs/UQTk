#ifndef GAUSSIANBASE_H_
#define GAUSSIANBASE_H_

#include "MUQ/Modeling/Distributions/Distribution.h"

namespace muq {
  namespace Modeling {

    /**
    @brief Defines an abstract Gaussian class.
    @class GaussianBase
    @ingroup Distributions
    @seealso Gaussian
    */
    class GaussianBase : public Distribution {
    public:

      /** Basic constructor for Gaussians with not additional hyperparameters.
          @param[in] dim The dimension of the Gaussian distribution.
      */
      GaussianBase(unsigned int dim);

      /** Basic constructor for Gaussians with hyperparameter inputs.
          @param[in] dim The dimension of the Gaussian distribution
          @param[in] hyperSizesIn A vector of integers specify the size of any
                                  additional inputs (e.g., mean, standard deviation).
      */
      GaussianBase(unsigned int dim,
                   Eigen::VectorXi const& hyperSizesIn);


      /** Construct a Gaussian with no hyperparameter inputs and a specified
          mean vector.
          @param[in] meanIn A vector containing the mean of the Gaussian distribution.
      */
      GaussianBase(Eigen::VectorXd const& meanIn);

      /** Construct a Gaussian with hyperparameters and a specified
          mean vector.
          @param[in] meanIn A vector containing the mean of the Gaussian distribution.
          @param[in] hyperSizesIn A vector of integers specify the size of any
                                  additional inputs (e.g., mean, standard deviation).
      */
      GaussianBase(Eigen::VectorXd const& meanIn,
                   Eigen::VectorXi const& hyperSizesIn);


      virtual ~GaussianBase() = default;

      /**
      @return The dimension of the distribution.
      */
      virtual unsigned int Dimension() const;

      /**
        Applies the covariance matrix to a vector \f$x\f$
        @param[in] x A reference to the input vector
        @return \f$\Sigma x\f$, where \f$\Sigma\f$ is the covariance matrix.
      */
      virtual Eigen::MatrixXd ApplyCovariance(Eigen::Ref<const Eigen::MatrixXd> const& x) const = 0;

      /**
        Applies the precision matrix (inverse covariance) to a vector \f$x\f$
        @param[in] x A reference to the input vector
        @return \f$\Sigma^{-1} x\f$, where \f$\Sigma\f$ is the covariance matrix.
      */
      virtual Eigen::MatrixXd ApplyPrecision(Eigen::Ref<const Eigen::MatrixXd> const& x) const = 0;

      /**
        Applies a matrix square root of the covariance to a vector x.
        @param[in] x A reference to the input vector
        @return \f$\Sigma^{1/2} x\f$, where \f$\Sigma^{1/2}\f$ is a square root of the covariance matrix (e.g., Cholesky factor).
      */
      virtual Eigen::MatrixXd ApplyCovSqrt(Eigen::Ref<const Eigen::MatrixXd> const& x) const = 0;

      /**
        Applies a matrix square root of the precision to a vector x.  If \f$\Sigma= L L^T\f$
        for a covariance matrix \f$\Sigma\f$ and a square root \f$L\f$, then
        \f$\Sigma^{-1} = L^{-T}L^{-1}\f$ is a decomposition of the precision matrix
        and \f$L^{-T}\f$ is a square root of the precision matrix.  This function
        returns the action of \f$L^{-T}\f$ on a vector \f$x\f$.
        @param[in] x A reference to the input vector
        @return \f$L^{-T} x\f$
      */
      virtual Eigen::MatrixXd ApplyPrecSqrt(Eigen::Ref<const Eigen::MatrixXd> const& x) const = 0;

      /**
        Returns a vector containing the mean of the Gaussian.
        @return A vector containing the mean.
      */
      virtual Eigen::VectorXd const& GetMean() const{return mean;};

      /**
        Set the mean of this Gaussian.  The new mean must be the same size as the old mean.
        @param[in] newMu The new mean.
      */
      virtual void SetMean(Eigen::VectorXd const& newMu);

      /**
        Return the log determinant of the covariance matrix.  Defaults to 0.0, so
        this should be overridden by child classes that are used in applications
        where the density needs to known completely, not only up to a normalizing
        constant.
      */
      virtual double LogDeterminant() const{return 0.0;};

      /** Process the hyperparameter inputs.  This should be overridden by
          children that have extra inputs.
          @param[in] params A vector of hyperparameter vectors.
      */
      virtual void ResetHyperparameters(ref_vector<Eigen::VectorXd> const& params){};

      /** Compute the gradient of the log density with respect to either the
          distribution input or the hyperparameters.  This should be overridden
          by child classes that have hyperparameters.  If not overridden, this
          function will assert false if wrt>0.
          @param[in] wrt Specifies the index of the variable we wish to take the gradient wrt.  If wrt==0, then the gradient should be taken wrt the input variable.
          @return The gradient of the log density.
      */
      virtual Eigen::VectorXd GradLogDensity(unsigned int wrt, ref_vector<Eigen::VectorXd> const& inputs) override;

    protected:

      /**
      Compute the log density.
      @param[inputs] A vector of extra hyperparameter vectors.
      @return A double containing the log density.
      */
      virtual double LogDensityImpl(ref_vector<Eigen::VectorXd> const& inputs) override;

      /// Sample the distribution
      virtual Eigen::VectorXd SampleImpl(ref_vector<Eigen::VectorXd> const& inputs) override;

      // Space to store the mean of the distribution
      Eigen::VectorXd mean;

    };
  } // namespace Modeling
} // namespace muq

#endif
