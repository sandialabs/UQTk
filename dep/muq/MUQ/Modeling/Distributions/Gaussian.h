#ifndef GAUSSIAN_H_
#define GAUSSIAN_H_

#include "MUQ/Modeling/Distributions/GaussianBase.h"

namespace muq {
  namespace Modeling {

    /**
    @class Gaussian
    @ingroup Distributions
    @seealso GaussianBase
    */
    class Gaussian : public GaussianBase {
    public:

      /// Are we specifying the mean, covariance matrix, or precision matrix
      enum Mode {

        /// We are specifying the covariance
        Covariance,

        /// We are specifying the precision
	      Precision
      };

      enum ExtraInputs {
        None           = 1 << 0,
        Mean           = 1 << 1,
        DiagCovariance = 1 << 2,
        DiagPrecision  = 1 << 3,
        FullCovariance = 1 << 4,
        FullPrecision  = 1 << 5,
      };
      typedef uint8_t InputMask;

      Gaussian(unsigned int dim,
               InputMask    extraInputs = ExtraInputs::None);

      /// Construct a Gaussian with scaled identity covariance/precision
      /**
	      @param[in] mu The mean
      */
      Gaussian(Eigen::VectorXd const& mu,
               InputMask              extraInputs = ExtraInputs::None);

      /// Construct a Gaussian by specifying both the mean and covariance or precision matrix
      /**
      @param[in] mu The mean vector
      @param[in] obj A matrix holding either the covariance or precision
      @param[in] mode A flag indicating whether the matrix should be treated as the covariance or precision
      */
      Gaussian(Eigen::VectorXd const& mu,
               Eigen::MatrixXd const& obj,
               Gaussian::Mode         mode = Gaussian::Mode::Covariance,
               InputMask              extraInputs = ExtraInputs::None);


      virtual ~Gaussian() = default;

      Mode GetMode() const{return mode;};

      virtual Eigen::MatrixXd ApplyCovariance(Eigen::Ref<const Eigen::MatrixXd> const& x) const override;
      virtual Eigen::MatrixXd ApplyPrecision(Eigen::Ref<const Eigen::MatrixXd> const& x) const override;

      virtual Eigen::MatrixXd ApplyCovSqrt(Eigen::Ref<const Eigen::MatrixXd> const& x) const override;
      virtual Eigen::MatrixXd ApplyPrecSqrt(Eigen::Ref<const Eigen::MatrixXd> const& x) const override;

      /// Get the covariance
      /**
	     @return The covariance
      */
      Eigen::MatrixXd GetCovariance() const;
      Eigen::MatrixXd GetPrecision() const;

      /// Set the covariance matrix
      /**
        @param[in] newcov The new covariance
      */
      void SetCovariance(Eigen::MatrixXd const& newCov);

      /// Set the precision matrix
      /**
        @param[in] newprec The new precision
      */
      void SetPrecision(Eigen::MatrixXd const& newPrec);

      Gaussian::InputMask GetInputTypes() const{return inputTypes;};

      /// Returns a new Gaussian distribution conditioned on a linear observation
      std::shared_ptr<Gaussian> Condition(Eigen::MatrixXd const& obsMat,
                                          Eigen::VectorXd const& data,
                                          Eigen::MatrixXd const& obsCov) const;

      void ResetHyperparameters(ref_vector<Eigen::VectorXd> const& params) override;

      virtual double LogDeterminant() const override{return logDet;};

    protected:

      /// Compute the distribution's scaling constant
      /**
	      @return Scaling constant
      */
      void ComputeNormalization();

      static Eigen::VectorXi GetExtraSizes(unsigned dim, InputMask extraInputs);

      static Gaussian::Mode ModeFromExtras(InputMask extraInputs);
      static void CheckInputTypes(InputMask extraInputs, Mode mode);

      /// Have we specified the covariance or the precision
      Gaussian::Mode mode;

      /// What form do the extra inputs take? Just the mean, or the mean and covariance?
      Gaussian::InputMask inputTypes;

      // Space to store either the covariance or precision matrix (depending on mode)
      Eigen::MatrixXd covPrec;

      // Space to tore the matrix square root of the covariance or precision
      Eigen::LLT<Eigen::MatrixXd> sqrtCovPrec;

      /// The log determinant of the covariance matrix
      double logDet;

    };
  } // namespace Modeling
} // namespace muq

#endif
