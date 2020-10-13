#ifndef NLOPTOPTIMIZATION_H_
#define NLOPTOPTIMIZATION_H_

#include <boost/property_tree/ptree.hpp>

#include <nlopt.h>

#include "MUQ/Optimization/Optimizer.h"
#include "MUQ/Optimization/CostFunction.h"

namespace muq {
  namespace Optimization {

    class NLoptOptimizer : public Optimizer {
    public:

      NLoptOptimizer(std::shared_ptr<CostFunction> cost,
                     boost::property_tree::ptree const& pt);

      virtual ~NLoptOptimizer();

      virtual std::pair<Eigen::VectorXd, double>
        Solve(std::vector<Eigen::VectorXd> const& inputs) override;

    private:


      /// Evaluate either the cost function or a constraint
      /**
         @param[in] n The size of the input
         @param[in] x The current point
         @param[out] grad The gradient of the cost/constraint
         @param[in] f_data A CostHelper
         \return The cost/constraint value
      */
      static double Cost(unsigned int n,
                         const double* x,
                         double* grad,
                         void* f_data);


      static void Constraint(unsigned int m,
                             double* result,
                             unsigned int n,
                             const double* x,
                             double* grad,
                             void* f_data);


      /// Override the evaluate impl method (solve the optimization problem)
      /**
	 @param[in] args The first input is the variable we are optimizing over
       */
      virtual void
      EvaluateImpl(muq::Modeling::ref_vector<boost::any> const& inputs) override;

      /// Get the NLOPT algorithm we are using
      /**
	 @param[in] alg User-input algorithm
	 \return The NLOPT algorithm
       */
      nlopt_algorithm NLOptAlgorithm(std::string const& alg) const;

      /// The algorithm used to solve the problem
      const nlopt_algorithm algorithm;

      /// True: minimize the cost function, False: maximize the cost function
      const bool minimize = true;


    }; // class NLoptOptimizer

  } // namespace Optimization
} // namespace muq

#endif
