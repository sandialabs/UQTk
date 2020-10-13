#ifndef SMOLYAKESTIMATOR_H
#define SMOLYAKESTIMATOR_H

#include "MUQ/Modeling/ModPiece.h"
#include "MUQ/Modeling/Flann/FlannCache.h"

#include "MUQ/Utilities/MultiIndices/MultiIndexSet.h"

#include <vector>
#include <boost/property_tree/ptree.hpp>

namespace muq {
namespace Approximation {

  template<typename EstimateType>
  class SmolyakEstimator {
  public:

    SmolyakEstimator(std::shared_ptr<muq::Modeling::ModPiece> const& modelIn);

    virtual ~SmolyakEstimator() = default;

    /** This is the main function to constructing static or adaptive Smolyak
        estimates.
    */
    virtual EstimateType Compute(std::shared_ptr<muq::Utilities::MultiIndexSet> const& fixedSet,
                                 boost::property_tree::ptree                           options = boost::property_tree::ptree());


    /** To be called after Compute, this will continue to refine the estimate
        until the stopping criteria in options is met.
    */
    virtual EstimateType Adapt(boost::property_tree::ptree options);

    /** Returns the current estimate of the global error in the Smolyak
        approximation.
    */
    virtual double Error() const{return globalError;};

    /** Returns the number of model evaluations this factory has performed. */
    virtual unsigned int NumEvals() const{return numEvals;};

    /** Returns the history of the global error for each adaptation iteration.*/
    virtual std::vector<double> ErrorHistory() const{return errorHistory;};

    /** Returns the cumulative number of evaluations after each adaptation iteration.
    Note that when Compute() is called, the number of evaluations is reset but the
    evaluation cache is not.  If Compute() is called more than once, this may lead
    to counterintuitive evaluation histories because the points that were evaluated
    in the first call to Compute() will not be added to the history.
    */
    virtual std::vector<int> EvalHistory() const{return evalHistory;};

    /** Returns the cumulative runtime (in seconds) for each adaptation iteration. */
    virtual std::vector<double> TimeHistory() const{return timeHistory;};

    /** Return the points that were evaluated during each adaptation step. */
    virtual std::vector<std::vector<Eigen::VectorXd>> PointHistory() const;

    /** Returns the terms used in the estimate as the a function of adaptation iterations. */
    virtual std::vector<std::vector<std::shared_ptr<muq::Utilities::MultiIndex>>> TermHistory() const{return termHistory;};

  protected:

    virtual void Reset();

    virtual void AddTerms(std::shared_ptr<muq::Utilities::MultiIndexSet> const& fixedSet);

    virtual void AddTerms(std::vector<std::shared_ptr<muq::Utilities::MultiIndex>> const& fixedSet);

    /** Updates the local error indicators for terms on the leading edge. */
    virtual void UpdateErrors();

    /** Refine the approximation by adding to the computed terms.  Returns true
        if any new terms were added.
    */
    virtual bool Refine();

    /** Evaluates the model at specified points in the cache and saves the results
        to the evalCache vector.
        @param[in] ptsToEval A set of indices into the cache with points that need evaluation.
    */
    virtual void EvaluatePoints(std::set<unsigned int> const& ptsToEval);

    /** Computes the locations where the model will need to be evaluated in order
        to construct a single tensor-product estimate.  For example, in the
        Smolyak quadrature setting, this function will return the quadrature points
        coming from the tensor product quadrature rule defined by the multiindex.

        This function works in tandem with the ComputeOneTerm function, which takes
        model evaluations and actually returns the tensor product estimate.  These
        two functions are split to allow the SmolyakEstimator to handle any job
        scheduling or caching that may be needed for parallel model evaluations.
    */
    virtual std::vector<Eigen::VectorXd> OneTermPoints(std::shared_ptr<muq::Utilities::MultiIndex> const& multi) = 0;

    /**
      This function works in tandem with the OneTermPoints function.  After the model
      has been evaluated at the points returned by OneTermPoints, this function will
      compute an estimate with the new model evaluations.  In the quadrature setting,
      the estimate will be computed with a weighted sum of the model evaluations stored
      in the modEvals vector.
    */
    virtual EstimateType ComputeOneTerm(std::shared_ptr<muq::Utilities::MultiIndex>                const& multi,
                                        std::vector<std::reference_wrapper<const Eigen::VectorXd>> const& modEvals) = 0;

    /** Should compute sum(smolyVals[i] * smolyWeights[i]) and return the result*/
    virtual EstimateType ComputeWeightedSum(Eigen::VectorXd const& weights) const;
    virtual EstimateType ComputeWeightedSum() const;


    virtual EstimateType AddEstimates(double w1, EstimateType const& part1, double w2, EstimateType const& part2) const = 0;

    virtual double ComputeMagnitude(EstimateType const& estimate) const = 0;

    /// The model used to construct the approximations
    std::shared_ptr<muq::Modeling::ModPiece> model;

    /// Multiindices defining each tensor product term in the Smolyak approximation
    std::shared_ptr<muq::Utilities::MultiIndexSet> termMultis;

    /// Holds the history of the error.  Each component corresponds to an iteration
    std::vector<double> errorHistory;

    /// Holds the history of how many model evaluations have occured.  Each component corresponds to an adaptation iteration
    std::vector<int> evalHistory;

    /// Holds the history of how.  Each component corresponds to an adaptation iteration.
    std::vector<double> timeHistory;

    /// Indices in the cache for the points that were evaluated during each adaptation iteration.
    std::vector<std::set<unsigned int>> pointHistory;

    /// The terms that are added during each refinement
    std::vector<std::vector<std::shared_ptr<muq::Utilities::MultiIndex>>> termHistory;

    /// A cache of model evaluations
    muq::Modeling::DynamicKDTreeAdaptor<> pointCache;
    std::vector<Eigen::VectorXd> evalCache;

    int InCache(Eigen::VectorXd const& input) const;
    Eigen::VectorXd const& GetFromCache(unsigned int index) const{return pointCache.m_data.at(index);};
    int AddToCache(Eigen::VectorXd const& newPt);
    int CacheSize() const{return pointCache.m_data.size();};

    const double cacheTol = 10.0*std::numeric_limits<double>::epsilon(); // <- points are considered equal if they are closer than this

    /// Tolerance on the time (in seconds) allowed to continue adapting
    double timeTol = std::numeric_limits<double>::infinity();

    /// Tolerance on the global error indicator to continue adapting
    double errorTol = 0.0;

    /// Tolerance on the maximum number of evaluations
    unsigned int maxNumEvals = std::numeric_limits<unsigned int>::max();

    /// The number of model evaluations that have been performed
    unsigned int numEvals=0;

    struct SmolyTerm {

      // Value of the tensor product approximation for one term in the Smolyak expansion
      EstimateType val;

      /* Weight on this tensor product approximation -- will change as terms are
         added to the Smolyak rule.
      */
      double weight = 0.0;

      /* Has val been computed for this term yet? */
      bool isComputed = false;

      /* Whether or not this term is part of the "old" set as defined by G&G */
      bool isOld = false;

      /* Is this term needed for the estimate or an error indicator? */
      bool isNeeded = false;

      // A local error indicator.  See eq (5.1) of Conrad and Marzouk
      double localError = -1.0;

      /* A vector containing the indices of points in the evaluation cache
        (see evalCache and pointCache variables) that are needed to compute the
        value for this term.  This is used for lazy evaluation and helps avoid
        reevaluation of the model.
      */
      std::vector<unsigned int> evalInds;

      Eigen::VectorXd diffWeights;

    };

    std::vector<SmolyTerm> terms;

    double globalError;

    const double nzTol = 10.0*std::numeric_limits<double>::epsilon();

  }; // class SmolyakEstimator

} // namespace muq
} // namespace Approximation


#endif
