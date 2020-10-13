#ifndef DISTRIBUTION_H_
#define DISTRIBUTION_H_

#include "MUQ/Modeling/WorkPiece.h"
#include "MUQ/Utilities/VariadicMacros.h"

#include "MUQ/Utilities/Exceptions.h"

#include <Eigen/Core>
#include <vector>

namespace muq {
  namespace Modeling {

    class Density;
    class RandomVariable;

    class Distribution : public std::enable_shared_from_this<Distribution> {
    public:
      friend class Density;
      friend class RandomVariable;

      Distribution(unsigned int           varSizeIn,
                   Eigen::VectorXi const& hyperSizesIn = Eigen::VectorXi()) : varSize(varSizeIn),
                                                                              hyperSizes(hyperSizesIn){};

      virtual ~Distribution() = default;

      /// Evaluate the log-density
      /**
	      If known, the log-density should be implemented by a child in the LogDensityImpl class.
	      @param[in] inputs the vector of inputs to the log-density
	      @return The log density
      */
      virtual double LogDensity(ref_vector<Eigen::VectorXd> const& inputs);
      virtual double LogDensity(std::vector<Eigen::VectorXd> const& inputs){return LogDensity(ToRefVector(inputs));};;
      VARIADIC_TO_REFVECTOR(LogDensity, Eigen::VectorXd, double);


      virtual Eigen::VectorXd GradLogDensity(unsigned int wrt, std::vector<Eigen::VectorXd> const& inputs){return GradLogDensity(wrt, ToRefVector(inputs));};;
      virtual Eigen::VectorXd GradLogDensity(unsigned int wrt, ref_vector<Eigen::VectorXd> const& inputs);

      template<typename... Args>
	    Eigen::VectorXd GradLogDensity(unsigned int wrt, Args... args) {
	      ref_vector<Eigen::VectorXd> inputs;
	      return GradLogDensity(wrt, inputs, args...);
      }

      /// Sample the distribution
      /**
	      Calls SampleImpl, the default behavior is to return boost::none
	      @param[in] inputs the vector of inputs to the log-density
	      @return A sample
       */
      Eigen::VectorXd Sample(ref_vector<Eigen::VectorXd> const& inputs);
      Eigen::VectorXd Sample(std::vector<Eigen::VectorXd> const& inputs){return Sample(ToRefVector(inputs));};

      /// Sample the distribution with no inputs
      /**
	      Allows the user to call Sample without any inputs
	      @return A sample
       */
      Eigen::VectorXd Sample();

      VARIADIC_TO_REFVECTOR(Sample, Eigen::VectorXd, Eigen::VectorXd);


      /** @brief Returns a density built from this distribution.
          @details The distribution class allows users to both evaluate the density
          corresponding to a probability distribution and draw a sample of the
          corresponding random variable.  Both of these actions are supported
          through the Evaluate function, where the first input to Evaluate specifies
          what action to perform: evaluate the density or sample the RV.
          However, sometimes we only want to focus on the density part of the
          distribution.  This function returns a Density object that only supports
          evaluting the log density and does not require the extra input specifying
          what type of action to perform.

          For example,
@code
std::shared_ptr<Distribution> dist = std::make_shared<Gaussian>(mu,cov);

Eigen::VectorXd x;
// ... initialize the point x

// To obtain the density through Evaluate, we need to call
boost::any densVal = dist->Evaluate(Distribution::Mode::EvaluateLogDensity, x);

// With Density however, we don't need to add the additional flag
boost::any densVal2 = dist->AsDensity()->Evaluate(x);

// or, analogously
std::shared_ptr<Density> dens = dist->AsDensity();
boost::any densVal3 = dens->Evaluate(x);
@endcode
      */
      std::shared_ptr<Density> AsDensity();

      /** @brief Returns a random variable built from this distribution.
          @details The distribution class allows users to both evaluate the density
          corresponding to a probability distribution and draw a sample of the
          corresponding random variable.  Both of these actions are supported
          through the Evaluate function, where the first input to Evaluate specifies
          what action to perform: evaluate the density or sample the RV.
          However, sometimes we only want to focus on the random variable part
          (i.e., the "Sample" function) of the distribution.  This function returns
           a RandomVariable object that only supports drawing realizations of the
           distribution, not evaluting the density. Because the  RandomVariable
           only does one thing, it does not require the extra input specifying
          what type of action to perform.

          For example,
@code
std::shared_ptr<Distribution> dist = std::make_shared<Gaussian>(mu,cov);

Eigen::VectorXd x;
// ... initialize the point x

// To obtain a sample, we need to call
boost::any sample = dist->Evaluate(Distribution::Mode::SampleDistribution, x);

// With the RandomVariable class however, we don't need to add the additional flag
boost::any sample2 = dist->AsVariable()->Evaluate(x);

// or, analogously
std::shared_ptr<RandomVariable> rv = dist->AsVariable();
boost::any sample3 = rv->Evaluate(x);
@endcode
      */
      std::shared_ptr<RandomVariable> AsVariable();

      const unsigned int    varSize;
      const Eigen::VectorXi hyperSizes;

    protected:

      /// Implement the log-density
      /**
	      If known, the log-density should be implemented by a child.  If it is not overridden then the default behavior is to return negative infinity (-1.0*std::numeric_limits<double>::infinity()).
	      @param[in] inputs the vector of inputs to the log-density
	      @return The log density
      */
      virtual double LogDensityImpl(ref_vector<Eigen::VectorXd> const& inputs){
        throw  muq::NotImplementedError(std::string("LogDensityImpl function has not been implemented for class ") + typeid(*this).name());
        return 0;
      };


      virtual Eigen::VectorXd GradLogDensityImpl(unsigned int wrt, ref_vector<Eigen::VectorXd> const& inputs);

      /// Sample the distribution
      /**
	      Should be overwritten by a child.  The default behavior is to return boost::none
      */
      virtual Eigen::VectorXd SampleImpl(ref_vector<Eigen::VectorXd> const& inputs){
        throw  muq::NotImplementedError(std::string("SampleImpl function has not been implemented for class ") + typeid(*this).name());
        return Eigen::VectorXd();
      };

      ref_vector<const Eigen::VectorXd> ToRefVector(std::vector<Eigen::VectorXd> const& anyVec) const;

    private:

      template<typename... Args>
	    Eigen::VectorXd GradLogDensity(unsigned int wrt, ref_vector<Eigen::VectorXd>& inputs, Eigen::VectorXd const& in, Args... args) {
      	inputs.push_back(std::cref(in));
      	return GradLogDensity(wrt, inputs, args...);
      }
    };
  } // namespace Modeling
} // namespace muq

#endif
