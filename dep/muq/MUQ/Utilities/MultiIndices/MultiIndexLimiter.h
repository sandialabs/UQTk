#ifndef MULTIINDEXLIMITER_H_
#define MULTIINDEXLIMITER_H_

#include <memory>

#include "MUQ/Utilities/MultiIndices/MultiIndex.h"

namespace muq{
namespace Utilities{
  
  /** @class MultiIndexLimiter
      @ingroup MultiIndices
      @brief An abstract base class for multi index limiters
      @details When constructing or adapting MultiIndex sets, some applications may need to impose additional requirements on what terms should be included.  This class provides a mechanism for adding that information to the MultiIndexSet.  The pure virtual function IsFeasible is overloaded by children to test individual multiindices.  If IsFeasible returns false for a given MultiIndex, the multiindex will not be included as an active member of the MultiIndexSet.
   @see muq::Utilities::MultiIndexSet muq::Utilities::MultiIndex
   */
  class MultiIndexLimiter{
  
  public:
    virtual ~MultiIndexLimiter() = default;
    
    /** This function is overloaded by children to define what terms are included. */
    virtual bool IsFeasible(std::shared_ptr<MultiIndex> multi) const = 0;
    
  }; // class MultiIndexLimiter
  
  /** @class TotalOrderLimiter
      @ingroup MultiIndices
      @brief Provides a cap on the total-order allowed
      @details This limter only allows terms that satisfy \f$\|\mathbf{j}\|_1\leq p_U\f$, where \f$\mathbf{j}\f$ is the multiindex, and \f$p_U\f$ is a nonnegative integer passed to the constructor of this class.
   */
  class TotalOrderLimiter : public MultiIndexLimiter{

  public:
    
    TotalOrderLimiter(unsigned int totalOrderIn) : totalOrder(totalOrderIn){};
    virtual ~TotalOrderLimiter() = default;
    
    virtual bool IsFeasible(std::shared_ptr<MultiIndex> multi) const override {return (multi->Sum() <= totalOrder);};
    
  private:
    TotalOrderLimiter(){};
    
    unsigned int totalOrder;
    
  }; // class TotalOrderLimiter


  /** @class DimensionLimiter
   @ingroup MultiIndices
   @brief Provides bounds on what dimensions are allowed to have nonzero values.
   @details This limiter only allows terms that satisfy \f$\mathbf{j}_d = 0 \f$ for \f$d<D_L\f$ or \f$d>=D_L+M\f$ for a lower bound \f$D_L\f$ and length \f$M\f$.
   */
  class DimensionLimiter : public MultiIndexLimiter{
      
  public:
    
    DimensionLimiter(unsigned int lowerDimIn, unsigned int lengthIn) : lowerDim(lowerDimIn), length(lengthIn){};
    virtual ~DimensionLimiter() = default;
    
    virtual bool IsFeasible(std::shared_ptr<MultiIndex> multi) const override;
    
  private:
    DimensionLimiter(){};
    
    unsigned int lowerDim;
    unsigned int length;
    
  }; // class DimensionLimiter
  
  
  /** @class GeneralLimiter
   @ingroup MultiIndices
   @brief Checks if a multiindex is in another set
   @details This limter only allows terms that are present in another MultiIndexSet.  The limiting set is passed as an argument to the constructor.
   */
  /*
  class GeneralLimiter : public MultiIndexLimiter{
    friend class boost::serialization::access;
    
  public:
    GeneralLimiter(std::shared_ptr<MultiIndexSet> limitingSetIn) : limitingSet(limitingSetIn){};
    virtual ~GeneralLimiter() = default;
    
    virtual bool IsFeasible(std::shared_ptr<MultiIndex> multi) const override;
    
  private:
    GeneralLimiter(){};
    
    std::shared_ptr<MultiIndexSet> limitingSet;
    
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version){
      ar & boost::serialization::base_object<MultiIndexLimiter>(*this);
      ar & limitingSet;
    };
    
    };*/ // class GeneralLimiter
  
  /** @class MaxOrderLimiter
   @ingroup MultiIndices
   @brief Provides a cap on the maximum value of each component the multiindex
   @details This limter only allows terms that satisfy \f$\mathbf{j}_i\leq p_i\f$ for \f$i\in \{1,2,\ldots,D\}\f$, where \f$p\f$ is a vector of upper bounds.
   */
  class MaxOrderLimiter : public MultiIndexLimiter{
    
  public:
    MaxOrderLimiter(unsigned int maxOrderIn) : maxOrder(maxOrderIn){};
    MaxOrderLimiter(Eigen::VectorXi const& maxOrdersIn) : maxOrder(-1), maxOrders(maxOrdersIn), vectorMin(maxOrdersIn.minCoeff()){};
    
    virtual ~MaxOrderLimiter() = default;
    
    virtual bool IsFeasible(std::shared_ptr<MultiIndex> multi) const override;
    
  private:
    MaxOrderLimiter(){};
    
    int maxOrder;
    Eigen::VectorXi maxOrders;
    int vectorMin;
    
  }; // class MaxOrderLimiter
 
  /** @class NoLimiter
   @ingroup MultiIndices
   @brief Returns true for an multiindex
   @details This class is used as a default in many places where a limiter is not always needed.  IsFeasible will return true for any multiindex.
   */
 class NoLimiter : public MultiIndexLimiter{
   
  public:
    virtual ~NoLimiter() = default;
    virtual bool IsFeasible(std::shared_ptr<MultiIndex> multi) const override {return true;};
    
  }; // class NoLimiter
  
  /** @class AndLimiter
   @ingroup MultiIndices
   @brief Combines two limiters through an AND operation
   @details This class will return true if both limiters given to the constructor return true.
   */
 class AndLimiter : public MultiIndexLimiter{
   
  public:
    AndLimiter(std::shared_ptr<MultiIndexLimiter> limitA, std::shared_ptr<MultiIndexLimiter> limitB) : a(limitA), b(limitB){};
    virtual ~AndLimiter() = default;
    virtual bool IsFeasible(std::shared_ptr<MultiIndex> multi) const override {return (a->IsFeasible(multi)&&b->IsFeasible(multi));};
    
  private:
    AndLimiter(){};
    std::shared_ptr<MultiIndexLimiter> a, b;
    
  }; // class AndLimiter
  
  /** @class OrLimiter
   @ingroup MultiIndices
   @brief Combines two limiters through an OR operation
   @details This class will return true if either of the limiters given to the constructor return true.
   */
  class OrLimiter : public MultiIndexLimiter{
    
  public:
    OrLimiter(std::shared_ptr<MultiIndexLimiter> limitA, std::shared_ptr<MultiIndexLimiter> limitB) : a(limitA), b(limitB){};
    virtual ~OrLimiter() = default;
    virtual bool IsFeasible(std::shared_ptr<MultiIndex> multi) const override {return (a->IsFeasible(multi)||b->IsFeasible(multi));};
    
  private:
    OrLimiter(){};
    std::shared_ptr<MultiIndexLimiter> a, b;
    
  }; // class OrLimiter
  
  
  /** @class XorLimiter
   @ingroup MultiIndices
   @brief Combines two limiters through an XOR operation
   @details This class will return true if exactly one of the limiters given to the constructor returns true.
   */
  class XorLimiter : public MultiIndexLimiter{
    
  public:
    XorLimiter(std::shared_ptr<MultiIndexLimiter> limitA, std::shared_ptr<MultiIndexLimiter> limitB) : a(limitA), b(limitB){};
    virtual ~XorLimiter() = default;
    virtual bool IsFeasible(std::shared_ptr<MultiIndex> multi) const override {return (a->IsFeasible(multi)^b->IsFeasible(multi));};
    
  private:
    XorLimiter(){};
    std::shared_ptr<MultiIndexLimiter> a, b;
    
  }; // class XorLimiter
  
} // namespace muq
} // namespace Utilities

#endif
