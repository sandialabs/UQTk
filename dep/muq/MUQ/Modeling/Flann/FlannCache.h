#ifndef FLANNCACHE_H_
#define FLANNCACHE_H_

#include <deque>

#include <nanoflann.hpp>

#include "MUQ/Modeling/LinearAlgebra/AnyAlgebra.h"

#include "MUQ/Modeling/ModPiece.h"

namespace muq {
  namespace Modeling {

    class FlannCache;

    template <class Distance = nanoflann::metric_L2, typename IndexType = size_t>
    struct DynamicKDTreeAdaptor {

      friend class FlannCache;

      typedef DynamicKDTreeAdaptor<Distance,IndexType> self_t;
      typedef typename Distance::template traits<double,self_t>::distance_t metric_t;
      typedef nanoflann::KDTreeSingleIndexDynamicAdaptor< metric_t,self_t,-1,IndexType>  index_t;

      std::shared_ptr<index_t> index; //! The kd-tree index for the user to call its methods as usual with any other FLANN index.
      std::deque<Eigen::VectorXd> m_data;

      /// Constructor: takes a const ref to the vector of vectors object with the data points
	    inline DynamicKDTreeAdaptor(const int dim, const int leaf_max_size = 10) {
        index = std::make_shared<index_t>(dim, *this /* adaptor */, nanoflann::KDTreeSingleIndexAdaptorParams(leaf_max_size ) );
      }

	    inline void UpdateIndex(const int leaf_max_size = 10) {
	      assert(m_data.size()>0);
        index = std::make_shared<index_t>(m_data.at(0).size(), *this /* adaptor */, nanoflann::KDTreeSingleIndexAdaptorParams(leaf_max_size ) );
        index->addPoints(0, m_data.size()-1);
      }

      inline virtual ~DynamicKDTreeAdaptor() {}

      inline void add(Eigen::VectorXd const& newPt) {
        m_data.push_back(newPt);
	      index->addPoints(m_data.size()-1, m_data.size()-1);
      }

      /** Query for the \a num_closest closest points to a given point (entered as query_point[0:dim-1]).
    	 *  Note that this is a short-cut method for index->findNeighbors().
    	 *  The user can also call index->... methods as desired.
    	 * \note nChecks_IGNORED is ignored but kept for compatibility with the original FLANN interface.
    	 */
      inline std::pair<std::vector<IndexType>, std::vector<double>> query(Eigen::VectorXd const& query_point, const size_t num_closest, const int nChecks_IGNORED = 10) const {
        std::vector<IndexType> out_indices(num_closest);
        std::vector<double> out_distances_sq(num_closest);

        nanoflann::KNNResultSet<double,IndexType> resultSet(num_closest);
        resultSet.init(&out_indices[0], &out_distances_sq[0]);
        index->findNeighbors(resultSet, query_point.data(), nanoflann::SearchParams());

        return std::make_pair(out_indices, out_distances_sq);
      }

      inline const self_t & derived() const {
        return *this;
      }

      inline self_t& derived() {
        return *this;
      }

      // Must return the number of data points
      inline size_t kdtree_get_point_count() const {
        return m_data.size();
      }

      // Returns the dim'th component of the idx'th point in the class:
      inline double kdtree_get_pt(const size_t idx, int dim) const {
	      assert(idx<m_data.size());
	      assert(dim<m_data[idx].size());

	      return m_data[idx][dim];
      }

      // Optional bounding-box computation: return false to default to a standard bbox computation loop.
      //   Return true if the BBOX was already computed by the class and returned in "bb" so it can be avoided to redo it again.
      //   Look at bb.size() to find out the expected dimensionality (e.g. 2 or 3 for point clouds)
      template <class BBOX>
      inline bool kdtree_get_bbox(BBOX & /*bb*/) const {
        return false;
      }
    }; // end of DynamicKDTreeAdaptor


    /// Create a cache of model evaluations (input/output pairs)
    /**
       Caches the input/output pairs for a muq::Modeling::WorkPiece that has one input and one output.
     */
    class FlannCache : public ModPiece {
    public:

      /**
	 If we pass a nullptr (default), then only construct a cache of points (no outputs)
	 @param[in] function The function whose input/output pairs we want to cache
       */
      FlannCache(std::shared_ptr<ModPiece> function/* = nullptr*/);

      ~FlannCache();

      /// Determine if an entry is in the cache
      /**
	     @param[in] input Check if this input is in the cache
	      \return -1: the input vector does not have an entry in the cache, otherise: return the index of the input in the cache
       */
      int InCache(Eigen::VectorXd const& input) const;

      /// Add new points to the cache
      /**
	 @param[in] inputs The points to add
	 \return The function evaluation at each input point
       */
      std::vector<Eigen::VectorXd> Add(std::vector<Eigen::VectorXd> const& inputs);

      /// Add new points to the cache given the result
      /**
	 @param[in] inputs The points to add
	 @param[in] results The function evaluation at each input point
       */
      void Add(std::vector<Eigen::VectorXd> const& inputs, std::vector<Eigen::VectorXd> const& results);

      /// Get an input point from the cache
      /**
	 @param[in] index The point to returnn
	 \return The input point associated with the index
       */
      const Eigen::VectorXd at(unsigned int const index) const;

      /// Get an input point from the cache
      /**
	 @param[in] index The point to returnn
	 \return The input point associated with the index
       */
      Eigen::VectorXd at(unsigned int const index);

      /// Returns the model for a specific cache index
      Eigen::VectorXd const& OutputValue(unsigned int index) const;

      /// Add a new point to the cache
      /**
	 @param[in] input The entry we would like to add to the cache (if it is not there already)
	 \return The function result at that point
       */
      Eigen::VectorXd Add(Eigen::VectorXd const& input);

      /// Add a new point to the cache given the result
      /**
	 @param[in] input The entry we would like to add to the cache (if it is not there already)
	 @param[in] result The result of the function
   @return Returns the index of the added point
      */
      unsigned int Add(Eigen::VectorXd const& input, Eigen::VectorXd const& result);

      /// Remove point from the cache
      /**
	 @param[in] input The entry we would like to remove from the cache
       */
      void Remove(Eigen::VectorXd const& input);

      /// The index of the nearest neighbor
      /**
	 @param[in] point We want the index of the nearest neighbor that is closest to this point.
	 \return The index of the nearest neighbor
       */
      size_t NearestNeighborIndex(Eigen::VectorXd const& point) const;

      /// Find the \f$k\f$ nearest neighbors
      /**
	 @param[in] point The point whose nearest neighbors we want to find
	 @param[in] k We want to find this many nearest neighbors
	 @param[out] neighbors A vector of the \fk\f$ nearest neighbors
	 @param[out] result The output corresponding to the \f$k\f$ nearest neighbors
       */
      void NearestNeighbors(Eigen::VectorXd const& point,
                            unsigned int const k,
                            std::vector<Eigen::VectorXd>& neighbors,
                            std::vector<Eigen::VectorXd>& result) const;

      /// Find the \f$k\f$ nearest neighbors (don't bother getting the result too)
      /**
	 @param[in] point The point whose nearest neighbors we want to find
	 @param[in] k We want to find this many nearest neighbors
	 @param[out] neighbors A vector of the \fk\f$ nearest neighbors
       */
      void NearestNeighbors(Eigen::VectorXd const& point,
                            unsigned int const k,
                            std::vector<Eigen::VectorXd>& neighbors) const;

      /// Get the size of the cache
      /**
	     \return The size of the cache
       */
      unsigned int Size() const;

      /// Get the centroid of the cache
      /**
      \return The cache centroid
      */
      Eigen::VectorXd Centroid() const;

      /// Get the underlying function
      /**
      \return The function whose output we are storing
      */
      std::shared_ptr<ModPiece> Function() const;

    private:

      /// Update the centroid when a new point is added to the cache
      /**
        @param[in] point The new point
      */
      void UpdateCentroid(Eigen::VectorXd const& point);

      // The vector of previous results
      std::vector<Eigen::VectorXd> outputCache;

      virtual void EvaluateImpl(ref_vector<Eigen::VectorXd> const& inputs) override;

      /// The function whose input/outputs we are caching
      std::shared_ptr<ModPiece> function;

      /// The nearest neighbor index, used to perform searches
      std::shared_ptr<DynamicKDTreeAdaptor<>> kdTree;

      /// The centroid of the input addPoints
      /**
        The centroid is
        \f{equation}{
          \bar{x} = \frac{1}{N} \sum_{i=1}^N x_i,
        \f}
        where \f$N\f$ is the cache size.

        Note: the centroid is not necessarily in the cache.
      */
      Eigen::VectorXd centroid;
    };
  } // namespace Utilities
} // namespace muq

#endif
