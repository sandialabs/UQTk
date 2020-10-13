#ifndef OPTIMIZER_H_
#define OPTIMIZER_H_

#include <boost/property_tree/ptree.hpp>

#include "MUQ/Modeling/ModPiece.h"
#include "MUQ/Optimization/CostFunction.h"

namespace muq {
namespace Optimization {
  /// Solve an optimization problem
  /**
     \f{eqnarray}{
     c &=& \min{J(x; \theta_1, ..., \theta_1)} \        \
     f_i(x) &\leq& 0 \                                  \
     g_i(x) &=& 0
     \f}
  */
  class Optimizer : public muq::Modeling::WorkPiece {
  public:

    Optimizer(std::shared_ptr<CostFunction> cost,
                 boost::property_tree::ptree const& pt);

    virtual ~Optimizer();

    /// Add an inequality constraint to the optimization
    /**
       @param[in] ineq The constraint
    */
    virtual void AddInequalityConstraint(std::vector<std::shared_ptr<muq::Modeling::ModPiece>> const& ineq);

    /// Add an inequality constraint to the optimization
    /**
       @param[in] ineq The constraint
    */
    virtual void AddInequalityConstraint(std::shared_ptr<muq::Modeling::ModPiece> const& ineq);


    /// Clear all inequality constraints
    void ClearInequalityConstraint();

    /// Add an equality constraint to the optimization
    /**
       NOTE: the NLOPT algorithm used must be able to handle equality constraints
       @param[in] ineq The constraint
    */
    virtual void AddEqualityConstraint(std::vector<std::shared_ptr<muq::Modeling::ModPiece>> const& eq);

    /// Add an equality constraint to the optimization
    /**
       NOTE: the NLOPT algorithm used must be able to handle equality constraints
       @param[in] ineq The constraint
    */
    virtual void AddEqualityConstraint(std::shared_ptr<muq::Modeling::ModPiece> const& eq);

    /// Clear all equality constraints
    void ClearEqualityConstraint();
    
    /// Solve the optimization problem
    /**
       @param[in] inputs The first input is the variable we are optimizing over, second input are the
                  cost function parameters, and the third input are the constraint parameters 
       \return First: the argmin, second: the minimum cost
    */
    virtual std::pair<Eigen::VectorXd, double> Solve(std::vector<Eigen::VectorXd> const& inputs)=0;

  protected:

    /// The cost function that we are trying to minimize
    std::shared_ptr<CostFunction> opt;
    
    /// Inequality constraints
    std::vector<std::shared_ptr<muq::Modeling::ModPiece>> ineqConstraints;
    
    /// Equality constraints
    /**
       NOTE: the solver muq::Optimization::Optimization::algorithm must be able to handle equality constraints
    */
    std::vector<std::shared_ptr<muq::Modeling::ModPiece>> eqConstraints;
    
    /// Relative and absolute tolerances on the cost function value and on the difference between successive values of the state
    const double ftol_rel, ftol_abs, xtol_rel, xtol_abs;
    
    /// Tolerance on the constraints
    const double constraint_tol;
    
    /// Maximum number of cost function evaluations
    const unsigned int maxEvals;
  };

} // namespace Optimization
} // namespace muq

#endif
