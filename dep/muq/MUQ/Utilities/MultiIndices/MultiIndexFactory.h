#ifndef MULTIINDEXFACTORY_H_
#define MULTIINDEXFACTORY_H_

#include "MUQ/config.h"

#include "MUQ/Utilities/MultiIndices/MultiIndexSet.h"
#include "MUQ/Utilities/MultiIndices/MultiIndexLimiter.h"

namespace muq {
  namespace Utilities{

    /** @class MultiIndexFactory
     @ingroup MultiIndices
     @brief A factory class with static methods for generating MultiIndexSets
     @details Functions in this class are static and used to create MultiIndex sets with different structure. If the default behavior of these methods is insufficient to create a multiindex for your application, a custom MultiIndexLimiter can be used to ``filter" the results of one of these methods.
     @seealso muq::Utilities::MultiIndexLimiter
     */
    class  MultiIndexFactory{

    public:

      /** @brief Construct a total order limited MultiIndex
          @details The total order of a multiindex \f$\mbox{j}\f$ is given by the \f$\ell_1\f$ norm \f$\|\mathbf{j}\|_1\f$. This function creates the set of all MultiIndices such that \f$p_L \leq \|\mathbf{j}\|_1 \leq p_U\f$ for two nonnegative integers \f$p_L\f$ and \f$p_U\f$.  For example, in two dimensions, setting \f$p_L=0\f$ and \f$p_U=3\f$ would result in a set containing \f$ [0, 0], [0, 1], [0, 2], [0, 3], [1, 0], [1, 1], [1, 2], [2,0], [2,1], [3,0]\f$.
          @param[in] length The dimension of the multiindex \f$\mathbf{j}\f$.
          @param[in] maxOrder The upper bound on the total order \f$p_U\f$.
          @param[in] minOrder The lower bound on the total order \f$p_L\f$.  The default value is 0.
          @param[in] limiter A shared_ptr to a child of the abstract MultiIndexLimiter class.  This input can be used to further restrict what multiIndices are put in the set.  However, the default limiter is a NoLimiter, which has no effect.
          @return A shared_ptr to a MultiIndexSet containing the total order limited set.  Note that the MultiIndexSet does NOT use the same limiter passed to this function.
       */
      static std::shared_ptr<MultiIndexSet> CreateTotalOrder(unsigned int const length,
                                                             unsigned int const maxOrder,
                                                             unsigned int const minOrder = 0,
                                                             std::shared_ptr<MultiIndexLimiter> limiter = std::make_shared<NoLimiter>());


      static std::vector<std::shared_ptr<MultiIndexSet>> CreateTriTotalOrder(unsigned int const length,
                                                                             unsigned int const maxOrder,
                                                                             unsigned int const minOrder = 0,
                                                                             std::shared_ptr<MultiIndexLimiter> limiter = std::make_shared<NoLimiter>());

      /** @brief Construct a general hyperbolic MultiIndex
          @details Following <a href="http://dl.acm.org/citation.cfm?id=1931121">Blatman and Sudret 2011</a>, this function creates a hyperbolic index set.  Define the \f$q\f$-norm of the multindex as \f[ \|\mathbf{j}\|_q = \left(\sum_{d=1}^Dj_d^q\right)^{1/q} \f] where \f$D\f$ is the dimension of the multiindex and \f$q\f$ is a positive scalar.  This function creates a set such that \f$ \|\mathbf{j}\|_q \leq p_U\f$ for some nonnegative integer \f$p_U\f$.  When \f$q=1\f$, we obtain a total order limited set and as \f$q\rightarrow\infty\f$, we would obtain the full tensor set.  Notice that in general \f$q<1\f$ to ensure that the resulting set is smaller than the total order set.
          @param[in] length The dimension of the multiindex \f$\mathbf{j}\f$.
          @param[in] maxOrder The upper bound \f$p_U\f$.  Notice that this is also an upper bound on the maximum value in the multiindex.
          @param[in] q A scalar defining the norm used.  The default value is 0.5.
          @param[in] limiter A shared_ptr to a child of the abstract MultiIndexLimiter class.  This input can be used to further restrict what multiIndices are put in the set.  However, the default limiter is a NoLimiter, which has no effect.
       */
      static std::shared_ptr<MultiIndexSet> CreateHyperbolic(unsigned int const length,
                                                             unsigned int const maxOrder,
                                                             const double q = 0.5,
                                                             std::shared_ptr<MultiIndexLimiter> limiter = std::make_shared<NoLimiter>());


      static std::vector<std::shared_ptr<MultiIndexSet>> CreateTriHyperbolic(unsigned int const length,
                                                                             unsigned int const maxOrder,
                                                                             const double q = 0.5,
                                                                             std::shared_ptr<MultiIndexLimiter> limiter = std::make_shared<NoLimiter>());

      /** @brief Construct a full tensor product multiindex set.
          @details The output of this function is a multiindex set containing a full tensor product of the one dimensional orders.  For example, a two dimensional set with order 2 would include \f$ [0, 0], [0, 1], [0, 2], [1, 0], [1, 1], [1, 2], [2, 0], [2,1], [2,2]\f$.  Notice that full tensor sets are significantly larger than total order limited sets.  In this function, the upper bound on the order is the same for each dimension.  For constructing tensor product sets with a different upper bound for each dimension, @see CreateFullTensor(const Eigen::RowVectorXi& orders).
          @param[in] length The dimension of the multiindex
          @param[in] order The upper bound on the order for all dimensions.
          @return A shared_ptr to a MultiIndexSet containing the tensor product set.  Note that the MultiIndexSet does NOT use the same limiter passed to this function.
       */
      static std::shared_ptr<MultiIndexSet> CreateFullTensor(unsigned int const length,
                                                             unsigned int const order,
                                                             std::shared_ptr<MultiIndexLimiter> limiter = std::make_shared<NoLimiter>());

      /** @brief Construct a full tensor product multiindex set with the potential for different orders in each dimension.
       @details Like the CreateFullTensor(int,int) function, the output of the this function is a tensor product multiindex set.  However, this function allows each dimension of the multiindex to have a different maximum.  For example, if [1,2] was passed as the input to this function, the returned multiindex set would include
           \f$[0, 0], [0, 1], [0, 2], [1, 0], [1, 1], [1, 2] \f$.
       @param[in] orders A row vector of nonnegative unsigned integers dictating the maximum order in each dimension.  Note that the maximum order is inclusive.
       @return A shared_ptr to a MultiIndexSet containing the tensor product set.  Note that the MultiIndexSet does NOT use the same limiter passed to this function.
       */
      static std::shared_ptr<MultiIndexSet> CreateFullTensor(const Eigen::RowVectorXi& orders,
                                                             std::shared_ptr<MultiIndexLimiter> limiter = std::make_shared<NoLimiter>());


      /** @brief Creates a single multiindex with one nonzero term.
          @details Constructs a MultiIndex of length totalDim that contains a
                   single nonzero element.  Component nonzeroDim is set to value
                   of order.
          @param[in] totalDim The length of the MultiIndex (i.e., number of components)
          @param[in] nonzeroDim The index of the single nonzero component.
          @param[in] order The value of the single nonzero component.
          @return A MultiIndex of the form [0,...,0,order,0,...,0]
      */                                              
      static std::shared_ptr<MultiIndex> CreateSingleTerm(int totalDim, int nonzeroDim, int order);


    private:

      static void RecursiveHyperbolicFill(const double maxOrderPow,
                                          std::shared_ptr<MultiIndexSet> output,
                                          unsigned int const currDim,
                                          Eigen::RowVectorXi &base,
                                          const double q,
                                          std::shared_ptr<MultiIndexLimiter> limiter);


      static void RecursiveTotalOrderFill(unsigned int const maxOrder,
                                          unsigned int const minOrder,
                                          std::shared_ptr<MultiIndexSet> output,
                                          unsigned int const currDim,
                                          Eigen::RowVectorXi &base,
                                          std::shared_ptr<MultiIndexLimiter> limiter);


      static void RecursiveTensor(const Eigen::RowVectorXi& orders,
                                  std::shared_ptr<MultiIndexSet> output,
                                  unsigned int const currDim,
                                  Eigen::RowVectorXi &base,
                                  std::shared_ptr<MultiIndexLimiter> limiter,
                                  bool allInactive);

    };// class MultiIndexFactory
  }// namespace Utilities
}// namespace muq


#endif
