#ifndef MULTIINDEX_H_
#define MULTIINDEX_H_

#include <unordered_map>
#include <memory>
#include <initializer_list>

#include <Eigen/Core>

namespace muq {
  namespace Utilities {

    class MultiIndexSet;

    ///
    /**
    @class MultiIndex
    @ingroup MultiIndices
    @brief A multi-index determines the powers of a multi-dimensional polynomial

    <p> In its simplest form, a multiindex is simply a vector of nonnegative
     integers, say \f$ \mathbf{j}=[j_1,j_2,\dots,j_D]\f$, where \f$D\f$ is a
      user-specified dimension.  In MUQ, these multiindices are used to define
       multivariate polynomial expansions.

      <p> The reason we use this class, instead of just storing the vectors
        directly, is that this class provides sparse storage scheme for the
        multiindex; only elements of \f$\mathbf{j}\f$ that are nonzero are
        actually stored.  This type of sparse storage is particularly
        advantageous for polynomial expansions that do not have a large number
        of cross terms, e.g., diagonal transport maps or highly anisotropic
        polynomial chaos expansions.
      </p>

      <p> This class is mostly used behind the scenes.  However, the GetVector()
          function may come in useful for users that need to extract the
          multiindex vector. </p>

     */
    class MultiIndex {
      friend class MultiIndexSet;

    public:

      /** Constructor that creates a multiindex of all zeros.
          @param[in] lengthIn The length (i.e., number of components) in the
           multiindex.
       */
      MultiIndex(unsigned lengthIn);

      /** Constructor that creates a multiindex with some default value.
       *          @param[in] lengthIn The length (i.e., number of components) in the
       *           multiindex.
       *          @param[in] val The value to be set for all entries.
       */
      MultiIndex(unsigned lengthIn, unsigned val);

      /** Takes a dense vector description of the multiindex and extracts the
          nonzero components.
          @param[in] indIn Row vector of unsigned integers defining the
                           multiindex.
       */
      MultiIndex(Eigen::RowVectorXi const& indIn);

      /** Allows users to intiailize the multiindex with curly braces.  For
          example, @code MultiIndex temp{1,0,2,3} #endcode would create a
          multiindex with length four and values 1, 0, 2, and 3.
      */
      MultiIndex(std::initializer_list<unsigned> const& indIn);

      /** Create a deep copy of MultiIndex pointed to by the input.
          @param[in] indIn A shared_ptr to a MultiIndex instance
          @return A shared_ptr to a new MultiIndex instance containing the same information as the input.
       */
      static std::shared_ptr<MultiIndex> Copy(std::shared_ptr<MultiIndex> const& indIn){return std::make_shared<MultiIndex>(*indIn);};

      /** Default irtual destructor.  */
      virtual ~MultiIndex() = default;

      /** Get the dense representation of this multiindex.
          @return A row vector of unsigned integers containing the multiindex.
       */
      Eigen::RowVectorXi GetVector() const;

      /** Get the total order of this multiindex: the \f$\ell_1\f$ norm.
       @return The sum of the nonzero components: the total order of this multiindex.
       */
      unsigned Sum() const{return totalOrder;};

      /** This function returns the maximum of this multiindex: the \f$\ell_\infty\f$ norm.
       @return The maximum value of the multiindex.
       */
      unsigned Max() const{return maxValue;};

      /** Use this function to set the value of the an entry in the multiindex.
       @param[in] ind The component of the multiindex to set (starting with 0).
       @param[in] val A non-negative value for the dim component of the multiindex.
       @return True if this function updated an already nonzero component, or false if this function added a new nonzero entry.
       */
      bool SetValue(unsigned ind, unsigned val);

      /** Obtain the a particular component of the multiindex.  Notice that this function can be slow for multiindices with many nonzero components.  The worst case performance requires \f$O(|\mathbf{j}|_0)\f$ integer comparisons, where \f$|\mathbf{j}|_0\f$ denotes the number of nonzero entries in the multiindex.
       @param[in] ind The component to return.
       @return The integer stored in component dim of the multiindex.
       */
      unsigned GetValue(unsigned ind) const;

      /** Change the length of the multiindex.  If the new length is shorter
          than the current length, any nonzero components after newLength will
          be removed.
          @param[in] newLength The new length of the multiindex
      */
      void SetLength(unsigned newLength);

      /** Returns the number of nonzero components in the multiindex.
      */
      unsigned int NumNz() const;

      std::string ToString() const;

      /** Get the number of components in the index.  When used to define a
          multivariate polynomial, this will return the dimension of the
          polynomial.
          @return The length of the multiindex.
      */
      unsigned GetLength() const{return length;};

      bool operator==(const MultiIndex &b) const;
      bool operator!=(const MultiIndex &b) const;
      bool operator<(const MultiIndex &b) const;
      bool operator>(const MultiIndex &b) const;

      MultiIndex& operator+=(const MultiIndex &b);
      MultiIndex& operator++();
      MultiIndex operator+(const MultiIndex &b) const;
      MultiIndex& operator-=(const MultiIndex &b);
      MultiIndex& operator--();
      MultiIndex operator-(const MultiIndex &b) const;

      std::unordered_map<unsigned, unsigned>::const_iterator GetNzBegin() const{return nzInds.begin();};
      std::unordered_map<unsigned, unsigned>::const_iterator GetNzEnd()   const{return nzInds.end();};

    private:

      unsigned int length;

      /// a vector holding pairs of (dimension,index) for nonzero values of index.
      std::unordered_map<unsigned, unsigned> nzInds;

      /// The maximum index over all nzInds pairs.
      unsigned maxValue;

      // the total order of the multiindex (i.e. the sum of the indices)
      unsigned totalOrder;

    }; // class MultiIndex


    struct MultiPtrComp{
      bool operator()(std::shared_ptr<MultiIndex> const& a, std::shared_ptr<MultiIndex> const& b) const{return (*a)<(*b);};
    };

    std::ostream& operator<< (std::ostream &out, const muq::Utilities::MultiIndex &ind);

  } // namespace Utilities
} // namespace muq



#endif
